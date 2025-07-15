/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;

/**
 * This command reads a RepGen list file and creates a complex of directories that can be used
 * for role training by the RoleTrainProcessor.  The list file is tab-delimited with headers,
 * and contains genome IDs in the "genome_id", genome names in the "genome_name" column, and
 * representative-genome IDs in the "rep_id" column.  We need to know the name of each representative
 * genome along with the list of genomes in each representative set.
 *
 * For each representative set that has sufficient size, we will build a subdirectory with the same name
 * as the representative-genome ID containing a "data.tbl" file that functions as an xmatrix for the
 * roles found in each represented genome.  Each subdirectory can be used to run the "RoleTrainProcessor"
 * and find missing roles.
 *
 * The positional parameters are the name of the role definition file, the name of the genome source containing
 * all the genomes, and the name of the output directory.  The RepGen list file should be presented as the standard
 * input.  In addition to the subdirectories for the representative sets, the output directory will contain an
 * "index.tbl" file that contains the ID, name, and set size of each representative genome for which there was
 * output and a copy of the role definition file called "roles.in.subsystems".  This last will be found as the
 * default file location by "RoleTrainProcessor".
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file name (if not STDIN)
 * -t	genome source type for master directory (default MASTER)
 * -m	minimum number of genomes for a representative set to qualify for output (default 100)
 * -M	maximum role occurrence count; higher occurrence counts are set to this value
 * 		(default 1, which leads to a presence/absence model)
 *
 * --clear	erase the output directory before processing
 *
 * @author Bruce Parrello
 *
 */
public class RoleBuildProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RoleBuildProcessor.class);
    /** role definition map */
    private RoleMap roleMap;
    /** input stream to use */
    private InputStream inputStream;
    /** map of genome IDs to names */
    private Map<String, String> genomeNames;
    /** map of representative genome IDs to representative set members */
    private Map<String, Set<String>> repSets;
    /** master genome source */
    private GenomeSource genomes;
    /** initial size for genome set hashes */
    private int hashSize;

    // COMMAND-LINE OPTIONS

    /** input file (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "rep.list.tbl", usage = "input file, if not STDIN")
    private File inFile;

    /** minimum number of genomes for a representative set to qualify for output */
    @Option(name = "--min", aliases = { "-m" }, metaVar = "500",
            usage = "minimum number of genomes for a representative set to qualify for output")
    private int minGenomes;

    /** maximum role occurrence count */
    @Option(name = "--maxRoles", aliases = { "-M" }, metaVar = "5",
            usage = "maximum permissible role-occurrence count (higher will be flattened to this value)")
    private int maxRoles;

    /** type of genome source for master genome directory */
    @Option(name = "--source", aliases = { "--type", "-t" }, usage = "genome source type (PATRIC, MASTER, DIR)")
    private GenomeSource.Type sourceType;

    /** if specified, the output directory will be erased before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** role definition file */
    @Argument(index = 0, metaVar = "roles.in.subsystems", usage = "role definition file", required = true)
    private File roleFile;

    /** genome source */
    @Argument(index = 1, metaVar = "genomeDir", usage = "genome source directory or file", required = true)
    private File genomeDir;

    /** output directory name */
    @Argument(index = 2, metaVar = "outDir", usage = "output directory name", required = true)
    private File outDir;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.maxRoles = 1;
        this.minGenomes = 100;
        this.sourceType = GenomeSource.Type.MASTER;
        this.clearFlag = false;
    }

    @Override
    protected void validateParms() throws IOException, ParseFailureException {
        // Insure the tuning parameters are reasonable.
        if (this.maxRoles < 1)
            throw new FileNotFoundException("Maximum number of role occurrences must be at least 1.");
        if (this.minGenomes < 10)
            throw new FileNotFoundException("Minimum set size must be at least 10.");
        // Compute the set size from the genome minimum.
        this.hashSize = this.minGenomes * 2;
        // Configure the input stream.
        if (this.inFile == null) {
            log.info("Genome specifications will be taken from the standard input.");
            this.inputStream = System.in;
        } else if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
        else {
            log.info("Genome specification will be taken from file {}.", this.inFile);
            this.inputStream = new FileInputStream(this.inFile);
        }
        // Process the role file.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " not found or unreadable.");
        log.info("Loading role definitions from {}.", this.roleFile);
        this.roleMap = RoleMap.load(this.roleFile);
        // Set up the genome source.
        if (! this.genomeDir.exists())
            throw new FileNotFoundException("Genome input source " + this.genomeDir + " does not exist.");
        this.genomes = this.sourceType.create(this.genomeDir);
        log.info("{} genomes found in {}.", this.genomes.size(), this.genomeDir);
        // Set up the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else
            log.info("Output will be to directory {}.", this.outDir);
        // Copy over the role definition file.
        File newRoleFile = new File(this.outDir, "roles.in.subsystems");
        FileUtils.copyFile(this.roleFile, newRoleFile);
        log.info("Role definitions copies to {}.", newRoleFile);
    }

    @Override
    protected void runCommand() throws Exception {
        try {
            // Initialize the hashes.  We want to process the representative sets in a recognizable sorted order,
            // so we use a tree map for them.
            this.repSets = new TreeMap<String, Set<String>>();
            this.genomeNames = new HashMap<String, String>(1000);
            // Loop through the input file, filling the hashes.
            try (TabbedLineReader genomeStream = new TabbedLineReader(this.inputStream)) {
                int genomeCount = 0;
                int genomeIdCol = genomeStream.findField("genome_id");
                int nameCol = genomeStream.findField("genome_name");
                int repCol = genomeStream.findField("rep_id");
                for (TabbedLineReader.Line line : genomeStream) {
                    String genomeId = line.get(genomeIdCol);
                    String repId = line.get(repCol);
                    // If the genome represents itself, save the name.
                    if (repId.contentEquals(genomeId))
                        this.genomeNames.put(repId, line.get(nameCol));
                    // Add the genome to the representative set.
                    Set<String> repSet = this.repSets.computeIfAbsent(repId, x -> new HashSet<String>(this.hashSize));
                    repSet.add(genomeId);
                    genomeCount++;
                    if (log.isInfoEnabled() && genomeCount % 100 == 0)
                        log.info("{} genomes read, {} representative sets found.", genomeCount, this.repSets.size());
                }
                log.info("{} total genomes in {} representative sets.", genomeCount, this.repSets.size());
            }
            // Create the index file.
            try (PrintWriter writer = new PrintWriter(new File(this.outDir, "index.tbl"))) {
                writer.println("rep_id\tname\tcount");
                // Now loop through the representative sets.
                for (Map.Entry<String, Set<String>> repSetEntry : this.repSets.entrySet()) {
                    var repSet = repSetEntry.getValue();
                    int repSize = repSet.size();
                    if (repSize >= this.minGenomes) {
                        String repId = repSetEntry.getKey();
                        String repName = this.genomeNames.computeIfAbsent(repId, x -> "Unknown genome " + x);
                        writer.format("%s\t%s\t%d%n", repId, repName, repSize);
                        this.buildRepSetDir(repId, repName, repSet);
                    }
                }
            }
        } finally {
            // Insure the input file is closed (if any).
            if (this.inFile != null)
                this.inputStream.close();
        }

    }

    /**
     * Create the classifier subdirectory for a single representative set.  Our basic strategy will be
     * to read through the genomes.  For each genome, we create a count-map for the roles.  For the
     * entire representative set, we create a set of role IDs found.  These are combined to produce the
     * output "data.tbl" file.
     *
     * @param repId		ID of the representative genome
     * @param repName	name of the representative genome
     * @param repSet	collection of genomes in the set
     *
     * @throws IOException
     */
    private void buildRepSetDir(String repId, String repName, Collection<String> repSet) throws IOException {
        log.info("Processing {}-genome representative set for {}: {}", repSet.size(), repId, repName);
        // Attempt to create the output directory.
        File outSubDir = new File(this.outDir, repId);
        if (outSubDir.isDirectory())
            log.info("Directory {} already exists-- skipping this set.", outSubDir);
        else {
            log.info("Creating output directory {}.", outSubDir);
            FileUtils.forceMkdir(outSubDir);
            // This hash will map each genome ID to a count-map keyed on role ID.
            var genomeRoleCounts = new HashMap<String, CountMap<String>>(repSet.size() * 4 / 3);
            // This will be the set of role IDs found.  We want these sorted.
            var roleSet = new TreeSet<String>();
            // Now loop through the genomes.
            for (String genomeId : repSet) {
                Genome genome = this.genomes.getGenome(genomeId);
                log.info("Scanning genome {}.", genome);
                var roleCounts = new CountMap<String>();
                genomeRoleCounts.put(genomeId, roleCounts);
                // Get the pegs for this genome and find interesting roles.
                genome.getPegs().stream().flatMap(x -> x.getUsefulRoles(this.roleMap).stream())
                        .forEach(x -> roleCounts.count(x.getId()));
                // Merge this genome's roles into the master set.
                roleSet.addAll(roleCounts.keys());
            }
            File xmatFile = new File(outSubDir, "data.tbl");
            log.info("Writing xmatrix to {}.", xmatFile);
            try (PrintWriter writer = new PrintWriter(xmatFile)) {
                writer.println("genome_id\t" + StringUtils.join(roleSet, '\t'));
                for (String genomeId : repSet) {
                    var roleCounts = genomeRoleCounts.get(genomeId);
                    String xLine = roleSet.stream()
                            .map(x -> Integer.toString(Math.min(roleCounts.getCount(x), this.maxRoles)))
                            .collect(Collectors.joining("\t", genomeId + "\t", ""));
                    writer.println(xLine);
                }
            }
        }
    }

}
