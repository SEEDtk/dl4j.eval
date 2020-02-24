/**
 *
 */
package org.theseed.dl4j.eval;

import java.util.List;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseProcessor;

/**
 * This command reads a raw.table and a roles.to.use and builds a training file.  Both these files are
 * output by one of the various build-role-tables scripts.
 *
 * In addition, the file test.set.tbl should contain the IDs of the genomes to be used as the testing
 * set.
 *
 * The training file contains the roles in roles.to.use across the top as column headers and has one row
 * per genome.  In each column of a genome's record is the number of times the role occurs in that genome.
 *
 * The positional parameter is the name of the directory containing the input files.  The output file
 * will be placed in that directory.
 *
 * The command-line options are as follows:
 *
 * -h	display usage information
 * -v	display more detailed progress messages on the log
 *
 *
 * @author Bruce Parrello
 *
 */
public class BuildProcessor extends BaseProcessor {

    // FIELDS
    /** input useful role list file */
    private File rolesToUse;
    /** input role occurrence table */
    private File rawTable;
    /** testing set genomes */
    private Set<String> testSet;
    /** list of roles to use */
    private List<String> roleList;
    /** role counts per genome */
    private Map<String, CountMap<String>> genomeMap;

    // COMMAND-LINE OPTIONS

    /** input / working directory */
    @Argument(index = 0, usage = "input directory, also used for output file", required = true)
    private File modelDir;

    @Override
    protected void setDefaults() { }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.modelDir.isDirectory())
            throw new FileNotFoundException(this.modelDir + " does not appear to be a directory.");
        // Get the input file names.
        this.rolesToUse = new File(this.modelDir, "roles.to.use");
        this.rawTable = new File(this.modelDir, "raw.table");
        // Verify they exist.
        if (! this.rolesToUse.canRead())
            throw new FileNotFoundException(this.rolesToUse + " is not found or unreadable.");
        if (! this.rawTable.canRead())
            throw new FileNotFoundException(this.rawTable + " is not found or unreadable.");
        // Read in the testing set.
        try (TabbedLineReader testStream = new TabbedLineReader(new File(this.modelDir, "test.set.tbl"))) {
            this.testSet = new HashSet<String>();
            for (TabbedLineReader.Line line : testStream) {
                this.testSet.add(line.get(0));
            }
            log.info("{} genomes set aside for the testing set.", this.testSet.size());
        }
        return true;
    }

    @Override
    public void run() {
        try {
            // Read in the roles to use.
            log.info("Reading {}.", this.rolesToUse);
            this.roleList = Evaluator.readRolesToUse(rolesToUse);
            // This will map each genome to its roles.
            this.genomeMap = new HashMap<String, CountMap<String>>();
            // Now read the raw table.
            log.info("Reading {}.", this.rawTable);
            try (TabbedLineReader rawStream = new TabbedLineReader(this.rawTable, 3)) {
                for (TabbedLineReader.Line line : rawStream) {
                    // Count this role as belonging to this genome.
                    String genome = Feature.genomeOf(line.get(2));
                    String role = line.get(1);
                    CountMap<String> genomeCounts = this.genomeMap.computeIfAbsent(genome, k -> new CountMap<String>());
                    genomeCounts.count(role);
                }
            }
            log.info("{} genomes found in input.", genomeMap.size());
            // Now write the training set.
            try (PrintWriter outStream = new PrintWriter(new File(this.modelDir, "training.tbl"))) {
                // First we write the header.
                String roles = this.roleList.stream().collect(Collectors.joining("\t"));
                outStream.format("genome\t%s%n", roles);
                log.info("Writing testing set.");
                // Now we write the testing set genomes.
                for (String genome : this.testSet) {
                    outputGenome(outStream, genome);
                }
                // Finally the rest of the genomes.
                log.info("Writing training set.");
                for (String genome : this.genomeMap.keySet()) {
                    if (! this.testSet.contains(genome))
                        outputGenome(outStream, genome);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
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
        String counts = this.roleList.stream().map(r -> Integer.toString(genomeCounts.getCount(r))).collect(Collectors.joining("\t"));
        outStream.format("%s\t%s%n", genome, counts);
    }

}
