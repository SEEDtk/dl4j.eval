/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.utils.ICommand;

/**
 * This processor performs the evaluations on PATRIC genomes.  It will take as input an evaluation directory and a tab-delimited
 * file of genome IDs.  Representative-genome technology will be used for completeness checking.
 *
 * The positional parameters are the name of the evaluation directory and the name of the input file.
 *
 * The input file should be tab-delimited.  The "-c" option specifies the name or index (1-based) of the column containing
 * the genome ID.
 *
 * When we are done, the output directory will contain a I<genomeID><code>.out</code> file for each genome.  This file will
 * contain the quality numbers plus the predicted and actual counts for each role.  In addition, a
 * summary of the results will be placed in "summary.tbl" in the same directory.
 *
 * The command-line options are as follows.
 *
 * -v	show more detailed progress messages
 * -o	the name of the output directory (defaults to the current directory)
 * -c	name or index (1-based) of the input column containing the genome ID, the default is column 1
 * -b	number of genomes to process in each batch
 * -i	name of the input file; if omitted, the standard input is assumed
 *
 * --terse		do not write the individual output files, only the summary
 *
 * @author Bruce Parrello
 *
 */
public class P3EvalProcessor extends Evaluator implements ICommand {

    // FIELDS
    /** stream for reading the input */
    private TabbedLineReader inStream;
    /** actual input column index */
    private int colIdx;
    /** activate PATRIC connection */
    private Connection p3;

    // COMMAND-LINE OPTIONS

    /** input file name */
    @Option(name = "-i", aliases = { "--input" }, metaVar = "genomed.tbl", usage = "input file name (if not STDIN)")
    private File inFile;

    /** key column index */
    @Option(name = "-c", aliases = { "--col" }, metaVar = "genome_id", usage = "name or index (1-based) of the input genome ID column")
    private String colId;

    /** number of genomes to process in each chunk */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "100", usage = "number of genomes to process in each batch")
    private int batchSize;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.colId = "1";
        this.inFile = null;
        this.batchSize = 200;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.isHelp()) {
                parser.printUsage(System.err);
            } else if (this.validateParms()) {
                // Open the input file.
                if (this.inFile != null) {
                    log.info("Genome IDs will be read from {}.", this.inFile);
                    this.inStream = new TabbedLineReader(this.inFile);
                } else {
                    log.info("Genome IDs will be read from standard input.");
                    this.inStream = new TabbedLineReader(System.in);
                }
                // Compute the input column index.
                this.colIdx = this.inStream.findField(this.colId);
               // Denote we're ready to run.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
        return retVal;
    }

    @Override
    public void run() {
        try {
            // Connect to PATRIC.
            this.p3 = new Connection();
            // Read in the role maps.
            this.initializeData();
            // Create the arrays.
            this.allocateArrays(this.batchSize);
            // Now we loop through the input genomes, processing them in batches.  The genomes are read into memory
            // one at a time, because of all the ancillary data needed.
            int iGenome = 0;
            for (TabbedLineReader.Line line : this.inStream) {
                if (iGenome >= this.batchSize) {
                    // No room for the next genome, so process this batch.
                    this.processBatch();
                    iGenome = 0;
                }
                // Add this genome to the data structures.  If it's not found, that's a warning.
                String genomeId = line.get(this.colIdx);
                P3Genome genome = P3Genome.Load(p3, genomeId, P3Genome.Details.STRUCTURE_ONLY);
                if (genome == null) {
                    log.debug("Could not find genome {} -- skipped.", genomeId);
                } else {
                    this.processGenome(iGenome, genome);
                    iGenome++;
                }
            }
            if (iGenome > 0) {
                this.setnGenomes(iGenome);
                this.processBatch();
            }
            // Clean up processing.
            this.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * Process all the genomes accumulated in the input arrays.
     *
     * @throws IOException
     */
    private void processBatch() throws IOException {
        this.evaluateConsistency();
        this.writeOutput();
    }

}
