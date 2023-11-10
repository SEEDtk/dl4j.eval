/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeMultiDirectory;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.GroupSpec;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.RoleMatrix;
import org.theseed.proteins.kmers.reps.RepGenome;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.FloatList;

/**
 * This class processes a master genome directory and a genome quality report file to create a
 * rep-genome based completeness engine.
 *
 * The basic strategy is to extract the singly-occurring roles from each good genome, then
 * assign that genome to each applicable representative group.  Groups with 100 or more genomes
 * will be processed to find completeness roles.  The list of representative genomes in the
 * largest group will be processed to find completeness roles for the root.
 *
 * The positional parameters are the name of the evaluation directory, the name of the master genome
 * directory, and the names of the representative-genome database files.  The standard input should
 * contain the summary output from a full evaluation of the master directory in TEXT mode.  We only
 * care about the "Genome" and "Good" columns, as this yields us the list of IDs for the good genomes.
 *
 * The output will be to "comp.tbl" in the evaluation directory.
 *
 * The command-line options are as follows.
 *
 * -i	specifies a file name to use for the standard input
 * -m	minimum number of acceptable roles for a completeness set
 *
 * --target		target completeness fractions for computing marker roles, coded as a comma-delimited
 * 				list; if a particular fraction doesn't yield enough roles, the next one will be used;
 * 				if none work, then common roles will be sought
 * --common		minimum commonality fraction for computing common roles (default 0.95)
 * --minSize	minimum number of genomes for a group to be useful
 *
 * @author Bruce Parrello
 *
 */
public class CompletenessRolesProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CompletenessRolesProcessor.class);
    /** master genome directory */
    private GenomeMultiDirectory master;
    /** list of target completeness fractions */
    private FloatList targets;
    /** input stream reader */
    private TabbedLineReader inStream;
    /** input column index for genome IDs */
    private int idColIdx;
    /** input column index for good/bad flag */
    private int goodColIdx;
    /** role definition map */
    private RoleMap rolesInSubsystems;
    /** ordered set of genome groups */
    private SortedSet<GroupSpec> groups;
    /** root group */
    private GroupSpec rootGroup;
    /** root genome list */
    private Set<String> rootGenomes;
    /** ID for the seed protein role */
    private static final String SEED_ROLE = "PhenTrnaSyntAlph";
    /** seed protein of current genome */
    private ProteinKmers seedProteinKmers;
    /** singleton roles in current genome */
    private List<String> genomeRoles;

    // COMMAND-LINE OPTIONS

    /** input file (if not using STDIN) */
    @Option(name = "--input", aliases = { "-i" }, usage = "input file (if not STDIN)")
    private File inFile;

    /** target completeness fractions */
    @Option(name = "--target", usage="comma-delimited list of target completeness fractions for computing marker roles")
    private void setTargets(String targetString) {
        this.targets = new FloatList(targetString);
    }

    /** minimum commonality fraction */
    @Option(name="--common", usage="target commonality fraction for computing common roles")
    private double commonOccurrence;

    /** minimum number of roles for an acceptable role set */
    @Option(name="-m", aliases = { "--minRoles", "--min" }, usage="minimum number of roles for a set")
    private int minRoles;

    /** minimum number of genomes for an acceptable genome group */
    @Option(name="--minSize", aliases = { "--minRoles", "--min" }, usage="minimum number of genomes for an acceptable completeness group")
    private int minSize;

    /** evaluation directory name */
    @Argument(index = 0, metaVar = "Eval.XX", usage = "evaluation directory name", required = true)
    private File evalDir;

    /** master genome directory name */
    @Argument(index = 1, metaVar = "P3Master", usage = "name of PATRIC master genome directory", required = true)
    private File masterDir;

    /** repgen database file names */
    @Argument(index = 2, metaVar = "rep1db.ser rep2db.ser ...", required = true)
    private List<File> repDbFiles;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.targets = new FloatList(new double[] { 0.90, 0.80 });
        this.commonOccurrence = 0.95;
        this.minRoles = 20;
        this.minSize = 100;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (this.inFile == null) {
            this.inStream = new TabbedLineReader(System.in);
            log.info("Reading from standard input.");
        } else {
            this.inStream = new TabbedLineReader(this.inFile);
            log.info("Reading from {}.", this.inFile);
        }
        this.goodColIdx = this.inStream.findField("Good");
        this.idColIdx = this.inStream.findField("Genome");
        // Verify the tuning parameters.
        if (this.commonOccurrence < 0.0 || this.commonOccurrence > 1.0)
            throw new ParseFailureException("Target commonality fraction must be between 0 and 1.");
        if (this.minRoles <= 0)
            throw new ParseFailureException("Minimum role count must be positive.");
        if (this.minSize <= 0)
            throw new ParseFailureException("Minimum genome group size must be positive.");
        // Verify the evaluation directory.
        if (! this.evalDir.isDirectory())
            throw new FileNotFoundException("Evaluation directory " + this.evalDir + " not found or invalid.");
        File roleFile = new File(this.evalDir, "roles.in.subsystems");
        if (! roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + roleFile + " not found in evaluation directory.");
        this.rolesInSubsystems = RoleMap.load(roleFile);
        // Verify the genome input directory.
        if (! this.masterDir.isDirectory())
            throw new FileNotFoundException("Master genome directory " + this.masterDir + " not found or invalid.");
        this.master = new GenomeMultiDirectory(this.masterDir);
        // Verify the representative-genome files.
        for (File repDbFile : repDbFiles) {
            if (! repDbFile.canRead())
                throw new FileNotFoundException("Representative-genome database file " + repDbFile + " not found or unreadable.");
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        log.info("Processing input.");
        // First, we create a group for each representative-genome set.
        this.createGroups();
        // Loop through the input, processing each genome.
        for (TabbedLineReader.Line line : this.inStream) {
            String genomeId = line.get(this.idColIdx);
            // If the genome is good, process it to get the seed protein and singleton roles.  Only
            // proceed if successful.
            if (line.getFlag(this.goodColIdx) && this.processGenome(genomeId)) {
                // Here we have a genome we want to process.
                int included = 0;
                // If it is in the root group, add it there.
                if (this.rootGenomes.contains(genomeId)) {
                    this.rootGroup.recordGenome(genomeId, this.genomeRoles);
                    included++;
                }
                // Find all the other groups containing this genome.
                for (GroupSpec group : this.groups) {
                    if (group.isClose(this.seedProteinKmers)) {
                        group.recordGenome(genomeId, this.genomeRoles);
                        included++;
                    }
                }
                log.info("Genome {} added to {} groups.", genomeId, included);
            }
        }
        // Add the root group to the group list.
        this.groups.add(this.rootGroup);
        // Loop through the groups.  For each group, if it is big enough, produce output for it.
        File outFile = new File(this.evalDir, "comp.tbl");
        int outCount = 0;
        try (PrintWriter outStream = new PrintWriter(outFile)) {
            for (GroupSpec group : this.groups) {
                if (group.gCount() >= this.minSize) {
                    this.processGroup(group, outStream);
                    outCount++;
                }
            }
        }
        // Clean up the files.
        this.inStream.close();
        log.info("All done. {} completeness groups output.", outCount);
    }

    /**
     * Process the specified genome to determine its seed protein and singleton roles.
     *
     * @param genomeId	ID of the genome to process
     *
     * @return TRUE if successful, FALSE if the genome is unsuitable or not found
     */
    private boolean processGenome(String genomeId) {
        boolean retVal = false;
        Genome genome = this.master.get(genomeId);
        if (genome == null)
            log.warn("WARNING:  Genome {} not found in master database.", genomeId);
        else {
            log.info("Processing {}.", genome);
            String seedProtein = "";
            CountMap<String> roleCounts = new CountMap<String>();
            // Loop through the pegs, counting roles.  Save the sequence for the PheS.
            for (Feature feat : genome.getPegs()) {
                Collection<String> roles = feat.getUsefulRoles(this.rolesInSubsystems).stream()
                        .map(x -> x.getId()).collect(Collectors.toList());
                for (String role : roles) {
                    if (role.contentEquals(SEED_ROLE))
                        seedProtein = feat.getProteinTranslation();
                    roleCounts.count(role);
                }
            }
            // Make sure we found a seed protein.
            if (seedProtein == null || seedProtein.isEmpty())
                log.warn("WARNING: No seed protein found in {}:  genome skipped.", genome.getId());
            else {
                // Prepare to process this genome.
                this.genomeRoles = roleCounts.counts().stream().filter(x -> x.getCount() == 1)
                        .map(x -> x.getKey()).collect(Collectors.toList());
                this.seedProteinKmers = new ProteinKmers(seedProtein);
                retVal = true;
            }
        }
        return retVal;
    }

    /**
     * Loop through the representative-genome databases, creating groups.
     *
     * @throws IOException
     */
    private void createGroups() throws IOException {
        this.groups = new TreeSet<GroupSpec>();
        // The first representative-genome database is used to create the root group.
        this.rootGroup = new GroupSpec(GroupSpec.ROOT_GENOME, 0);
        this.rootGenomes = new HashSet<String>(8000);
        Iterator<File> iter = this.repDbFiles.iterator();
        File repFile = iter.next();
        RepGenomeDb repDb = RepGenomeDb.load(repFile);
        // Save the kmer size.  We will verify that it is the same in the other databases.
        ProteinKmers.setKmerSize(repDb.getKmerSize());
        log.info("Primary repGen DB is {} with kmer size {} and {} genomes.", repFile, repDb.getKmerSize(), repDb.size());
        // Now process this database.
        for (RepGenome rep : repDb) {
            this.rootGenomes.add(rep.getGenomeId());
            GroupSpec newGroup = new GroupSpec(rep, repDb.getThreshold());
            this.groups.add(newGroup);
        }
        // Loop through the remaining databases.
        while (iter.hasNext()) {
            repFile = iter.next();
            repDb = RepGenomeDb.load(repFile);
            log.info("Processing repGen DB {} with {} genomes.", repFile, repDb.size());
            if (repDb.getKmerSize() != ProteinKmers.kmerSize())
                throw new IllegalArgumentException("Invalid kmer size " + Integer.toString(repDb.getKmerSize()) + " in " +
                        repFile.toString() + ". All databases must use the same size.");
            // Process this database's groups.
            for (RepGenome rep : repDb)
                this.groups.add(new GroupSpec(rep, repDb.getThreshold()));
        }
    }

    /**
     * Output the data for a completeness group.
     *
     * @param group			completeness group to process
     * @param outStream		output writer to contain the completeness information
     */
    private void processGroup(GroupSpec group, PrintWriter outStream) {
        log.info("Creating completeness data for group {}.", group);
        // We are searching for a good role set.  We try marker roles at each of the
        // completeness levels in the float list, and if none of them work, we fall back
        // to common roles.  The role set found is put in here.
        Collection<String> roleSet = new ArrayList<String>();
        // Get the role matrix.
        RoleMatrix roleMtx = group.getMatrix();
        // Start with the first target value.
        this.targets.reset();
        // Loop through the targets.
        double target = 0.0;
        while (this.targets.hasNext() && roleSet.size() < this.minRoles) {
            target = targets.next();
            log.debug("Computing marker role set for target {}.");
            roleSet = roleMtx.getMarkerRoles(target);
        }
        if (roleSet.size() < this.minRoles) {
            // Here we failed to find a role set, so we look for common roles.
            log.info("Computing common role set at threshold {}.", this.commonOccurrence);
            roleSet = roleMtx.getCommonRoles(this.commonOccurrence);
            // Count the genomes that fail the lowest target.
            int fCount = 0;
            for (String genome : roleMtx.genomes()) {
                if (roleMtx.completeness(roleSet, genome) < target) fCount++;
            }
            log.info("{} genomes failed the lowest target.", fCount);
        }
        log.info("{} roles found for {}.", roleSet.size(), group.getId());
        // Write out the group header.
        outStream.println(group.getHeader());
        // Write out the roles.
        for (String role : roleSet) {
            outStream.println("   " + role);
        }
        // Write out the trailer.
        outStream.println("//");
    }

}
