/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Genome;
import org.theseed.utils.ICommand;

/**
 * This processor evaluates a single GTO.  This is less efficient than evaluating in a batch, but it
 * allows this program to be used in a genome-processing pipeline.  The GTO is read from the standard
 * input and written to the standard output, but both of these options can be overridden on the command
 * line.
 *
 * The output directory will contain the results of the evaluation, though they will also be stored in
 * the output GTO.  The main advantage here is that the output directory results can be parsed without
 * loading the entire genome.
 *
 * The positional parameter is the name of the evaluation directory.
 *
 * The command-line options are as follows.
 *
 * -v	show more detailed progress messages
 * -O	the name of the output directory (defaults to the current directory)
 * -s	maximum distance for a close protein in deep reporting (0 to 1)
 * -i	name of the input file (overrides STDIN)
 * -o	name of the output file (overrides STDOUT)
 *
 * --terse		do not write the individual output files, only the summary
 * --clear		clear the output directory before processing
 * --format		specify the output format-- HTML, DEEP, or TEXT
 * --ref		ID of a reference genome to use for all the evaluation reports
 *
 * @author Bruce Parrello
 *
 */
public class GtoEvalProcessor extends Evaluator implements ICommand {

    // COMMAND LINE

    /** input file name */
    @Option(name = "-i", aliases = { "--input" }, usage = "input file name (if not STDIN)")
    private File inFile;

    /** output file name */
    @Option(name = "-o", aliases = { "--output" }, usage = "output file name (if not STDOUT)")
    private File outFile;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.inFile = null;
        this.outFile = null;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.isHelp()) {
                parser.printUsage(System.err);
            } else if (this.validateParms()) {
                // Insure we can read the genome file, if there is one.
                if (this.inFile != null && ! this.inFile.canRead()) {
                    throw new FileNotFoundException("Input" + this.inFile + " does not exist or cannot be read.");
                }
                // Suppress the summary report.
                this.suppressSummary();
            }
            // Denote we're ready to run.
            retVal = true;
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
            // Read in the role maps.
            initializeData();
            // Read in the genome.
            log.info("Loading input genome.");
            Genome genome;
            if (this.inFile != null) {
                genome = new Genome(this.inFile);
            } else {
                genome = new Genome(System.in);
            }
            log.info("Analyzing genome {}.", genome.getId());
            // Allocate the arrays.
            this.allocateArrays(1);
            // Store the genome in the arrays.
            processGenome(0, genome);
            // Evaluate the consistency of the genomes.
            evaluateConsistency();
            // Write the results.
            writeOutput();
            // Now we need to write the GTO.  Get the version string.
            String version = this.getVersion();
            // Retrieve the evaluation report.
            GenomeStats gReport = this.getReport(0);
            log.info("Writing evaluated genome.");
            gReport.store(genome.getJson(), this.roleDefinitions, version);
            if (this.outFile != null) {
                genome.update(this.outFile);
            } else {
                genome.update(System.out);
            }
            // Finish processing.
            this.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
