/**
 *
 */
package org.theseed.dl4j.eval;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.counters.Shuffler;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.MarkerFile;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;
import org.theseed.subsystems.SubsystemRoleFactory;

/**
 * This command reads a PATRIC CoreSEED dump and builds training files.  The PATRIC CoreSEED dump consists
 * of an annotation directory and a compressed subsystem directory.  Additional information is taken from
 * the current CoreSEED evaluator directory.
 *
 * The CoreSEED evaluator directory contains "roles.init.tbl", which contains the fixed role IDs used to
 * initialize the role table, "parms.prm", which is the RandomForest configuration file, and "training.tbl",
 * which contains the IDs of the acceptable genomes, in order.
 *
 * The compressed subsystem directory is processed by the SubsystemRoleFactory object to create a role map.
 * This map is saved as the master "roles.in.subsystems" file.
 *
 * The annotation directory contains one file per genome.  Each genome file has the genome ID as its name.
 * The genomes are processed in order as they appear in "training.tbl".  Each genome file is tab-delimited,
 * without headers, and contains in each record a feature ID followed by a role description.  Only the roles
 * identified by the SubsystemRoleFactory are considered.  For each genome, we need a count of the number of
 * times each role occurs.  We also keep a global count of the number of genomes containing each role.
 *
 * The output training file contains the roles in roles.to.use across the top as column headers and has one row
 * per genome.  In each column of a genome's record is the number of times the role occurs in that genome.
 *
 * The positional parameters are the name of the PATRIC CoreSEED dump annotation directory, the name of the compressed
 * subsystem-directory file, and the name of the output directory.
 *
 * Use RolesProcessor in kmer.reps to build the completeness file comp.tbl.
 *
 * The command-line options are as follows:
 *
 * -h	display usage information
 * -v	display more detailed progress messages on the log
 *
 * --core		the name of the CoreSEED evaluator directory; the default is "CoreSEED/Eval.New" in the current
 * 				directory
 * --minOccur	the minimum number of genomes containing a role required for a role to be considered useful; the
 * 				default is 100
 * --maxMulti	the maxmimum number of times a role can be contained in a genome for it to be considered useful; the
 * 				default is 5
 * --clear		erase the output directory before processing
 * --test		name of a directory of GTOs; if specified, the GTOs will be used to create a testing set
 *
 * @author Bruce Parrello
 *
 */
public class BuildProcessor extends BaseProcessor {

    // FIELDS
    /** role counts per genome */
    private Map<String, CountMap<String>> genomeMap;
    /** ordered list of genome IDs */
    private List<String> genomeIDs;
    /** subsystem role map */
    private RoleMap subsystemRoles;
    /** number of genomes containing each role */
    private CountMap<String> roleGenomeCounts;
    /** set of bad roles (occur too often in a genome) */
    private Set<String> badRoles;
    /** role initialization file */
    private File roleInitFile;
    /** list of roles to use, in order */
    private List<String> roleList;

    // COMMAND-LINE OPTIONS

    /** coreSEED evaluation directory */
    @Option(name = "--core", metaVar = "Core/Eval", usage = "coreSEED evaluator directory")
    private File coreEvalDir;

    /** minimum number of genomes that must contain a role for it to be useful */
    @Option(name = "--minOccur", metaVar = "50", usage = "minimum number of genomes that must contain a role for it to be useful")
    private int minOccur;

    /** maximum number of role occurrences per genome for a role to be useful */
    @Option(name = "--maxMulti", metaVar = "3", usage = "maximum number of times a useful role can occur in a genome")
    private int maxMulti;

    /** if specified, the output directory will be erased before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** testing-set GTO directory */
    @Option(name = "--test", metaVar = "GTOrefs", usage = "directory of GTOs to use to create a testing set")
    private File testDir;

    /** input annotation directory */
    @Argument(index = 0, metaVar = "Annotations/0", usage = "input annotation directory", required = true)
    private File annoDir;

    /** input subsystem-directory tgz file */
    @Argument(index = 1, metaVar = "subsystems.tgz", usage = "gzipped tar file containing subsystem directory", required = true)
    private File subFile;

    /** output directory */
    @Argument(index = 2, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    @Override
    protected void setDefaults() {
        this.coreEvalDir = new File(System.getProperty("user.dir"), "CoreSEED/Eval.New");
        this.maxMulti = 5;
        this.minOccur = 100;
        this.clearFlag = false;
        this.testDir = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the coreSEED directory.
        if (! this.coreEvalDir.isDirectory())
            throw new FileNotFoundException(this.coreEvalDir + " is not found or is not a directory.");
        File trainFile = new File(this.coreEvalDir, "training.tbl");
        if (! trainFile.canRead())
            throw new FileNotFoundException(this.coreEvalDir + " does not have a readable training file.");
        this.genomeIDs = TabbedLineReader.readColumn(trainFile, "genome_id");
        log.info("{} good genome IDs found in {}.", this.genomeIDs.size(), trainFile);
        this.roleInitFile = new File(this.coreEvalDir, "roles.init.tbl");
        if (! this.roleInitFile.canRead())
            throw new FileNotFoundException(this.coreEvalDir + " does not have a readable role initialization file.");
        File parmFile = new File(this.coreEvalDir, "parms.prm");
        if (! parmFile.canRead())
            throw new FileNotFoundException(this.coreEvalDir + " does not have a readable parms.prm file.");
        // verify the testing directory.
        if (this.testDir != null && ! this.testDir.isDirectory())
            throw new FileNotFoundException(this.testDir + " is not found or invalid.");
        // Verify the input files.
        if (! this.annoDir.isDirectory())
            throw new FileNotFoundException(this.annoDir + " is not found or is notoes not a directory.");
        if (! this.subFile.canRead())
            throw new FileNotFoundException(this.subFile + " is not found or unreadable.");
        // Verify the tuning parameters.
        if (this.maxMulti < 1 || this.maxMulti > 10)
            throw new ParseFailureException("Maxmimum multiplicity must be between 1 and 10 inclusive.");
        if (this.minOccur < 1)
            throw new ParseFailureException("Minimum occurrences must be at least 1.");
        // Set up the output directory.
        if (this.outDir.isDirectory()) {
            if (! this.clearFlag)
                log.info("Evaluation files will be created in {}.", this.outDir);
            else {
                log.info("Erasing output directory {}.", this.outDir);
                FileUtils.cleanDirectory(this.outDir);
            }
        } else {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        }
        // Copy the parm file to the output directory.
        FileUtils.copyFileToDirectory(parmFile, this.outDir);
        // Create the version stamp.
        String versionStamp = this.outDir.getName() + "-" + DateTimeFormatter.ISO_DATE.format(LocalDate.now());
        log.info("Evaluator version stamp is {}.", versionStamp);
        MarkerFile.write(new File(this.outDir, "VERSION"), versionStamp);
        // Create the roles directory.
        File roleDir = new File(this.outDir, "Roles");
        if (! roleDir.isDirectory()) {
            log.info("Creating Roles subdirectory {}.", roleDir);
            FileUtils.forceMkdir(roleDir);
        }
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        // Initialize the output maps.
        this.badRoles = new HashSet<String>();
        this.roleGenomeCounts = new CountMap<String>();
        this.genomeMap = new HashMap<String, CountMap<String>>();
        // Compute the main role map.
        log.info("Processing subsystems.");
        this.subsystemRoles = SubsystemRoleFactory.processArchive(this.subFile, this.roleInitFile);
        log.info("{} unique roles found in subsystems.", this.subsystemRoles.size());
        // Write the roles-in-subsystems file.
        File rolesOutFile = new File(this.outDir, "roles.in.subsystems");
        this.subsystemRoles.save(rolesOutFile);
        log.info("{} created.", rolesOutFile);
        // Read through the annotations of the authorized genomes.
        for (String genomeID : this.genomeIDs) {
            // Check for this genome.
            File annoFile = new File(this.annoDir, genomeID);
            if (! annoFile.exists())
                log.warn("Skipping genome {}-- not found in {}.", genomeID, this.annoDir);
            else try (TabbedLineReader annoReader = new TabbedLineReader(annoFile, 2)) {
                // Accumulate the genome's role counts in here.
                CountMap<String> gCounts = new CountMap<String>();
                // Loop through the file, counting roles.  Note that the feature ID does not matter.
                for (TabbedLineReader.Line line : annoReader)
                    Feature.usefulRoles(this.subsystemRoles, line.get(1)).stream().forEach(x -> gCounts.count(x.getId()));
                // Check for bad roles and count good roles.
                for (CountMap<String>.Count count : gCounts.counts()) {
                    if (count.getCount() > this.maxMulti)
                        this.badRoles.add(count.getKey());
                    else
                        this.roleGenomeCounts.count(count.getKey());
                }
                // Save this genome's counts.
                log.info("{} roles found in genome {}.", gCounts.size(), genomeID);
                this.genomeMap.put(genomeID, gCounts);
            }
        }
        log.info("{} roles were counted.  {} were bad.", this.roleGenomeCounts.size(), this.badRoles.size());
        // Produce the roles-to-use list.
        this.roleList = this.roleGenomeCounts.sortedCounts().stream()
                .filter(x -> (! this.badRoles.contains(x.getKey()) && x.getCount() >= this.minOccur))
                .map(x -> x.getKey()).collect(Collectors.toList());
        log.info("{} roles will be used from {} genomes.", this.roleList.size(), this.genomeMap.size());
        File useFile = new File(this.outDir, "roles.to.use");
        try (PrintWriter writer = new PrintWriter(useFile)) {
            for (String roleId : this.roleList) {
                // For each role we output the number of genomes containing the role, and then the
                // number of times it occurs at each multiplicity.
                int gCount = this.roleGenomeCounts.getCount(roleId);
                int[] occurrences = new int[this.maxMulti + 1];
                this.genomeMap.values().stream().forEach(x -> occurrences[x.getCount(roleId)]++);
                String oCounts = Arrays.stream(occurrences).mapToObj(k -> Integer.toString(k)).collect(Collectors.joining(","));
                writer.format("%s\t%d\t%s%n", roleId, gCount, oCounts);
            }
        }
        // Create the role-id header.
        String roles = this.roleList.stream().collect(Collectors.joining("\t"));
        // Now write the training set.
        try (PrintWriter outStream = new PrintWriter(new File(this.outDir, "training.tbl"))) {
            // First we write the header.
            outStream.format("genome\t%s%n", roles);
            // Shuffle the genomes so we get a random distribution.
            Shuffler<String> genomeList = new Shuffler<String>(this.genomeMap.keySet());
            genomeList.shuffle(genomeList.size());
            log.info("Writing training.tbl.");
            // Now we write the genomes.
            for (String genome : genomeList)
                outputGenome(outStream, genome);
        }
        // If we have a testing directory, write the testing set.
        if (this.testDir != null) {
            try (PrintWriter outStream = new PrintWriter(new File(this.outDir, "testing.tbl"))) {
                // First we write the header.
                outStream.format("genome\t%s%n", roles);
                // Loop through the genomes.
                GenomeDirectory genomes = new GenomeDirectory(this.testDir);
                for (Genome genome : genomes) {
                    log.info("Processing testing genome {}.", genome);
                    CountMap<String> gCounts = new CountMap<String>();
                    for (Feature feat : genome.getPegs())
                        feat.getUsefulRoles(this.subsystemRoles).stream().forEach(x -> gCounts.count(x.getId()));
                    this.writeCounts(outStream, genome.getId(), gCounts);
                }
            }
        }
    }

    /**
     * Write the genome's counts to the specified output stream.
     *
     * @param outStream		target output stream
     * @param genome		ID of the genome to write
     */
    public void outputGenome(PrintWriter outStream, String genome) {
        CountMap<String> genomeCounts = this.genomeMap.get(genome);
        if (genomeCounts == null)
            throw new IllegalArgumentException("Genome " + genome + " not found in raw.table.");
        writeCounts(outStream, genome, genomeCounts);
    }

    /**
     * Write the specified genome counts to a training/testing file.
     *
     * @param outStream		output writer
     * @param genomeId		ID of the relevant genome
     * @param genomeCounts	map of role occurrence counts
     */
    protected void writeCounts(PrintWriter outStream, String genomeId, CountMap<String> genomeCounts) {
        String counts = this.roleList.stream().map(r -> Integer.toString(genomeCounts.getCount(r)))
                .collect(Collectors.joining("\t"));
        outStream.format("%s\t%s%n", genomeId, counts);
    }

}
