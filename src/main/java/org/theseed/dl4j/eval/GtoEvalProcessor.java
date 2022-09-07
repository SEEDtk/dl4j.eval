/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.kohsuke.args4j.Option;
import org.theseed.dl4j.eval.reports.EvalReporter;
import org.theseed.dl4j.eval.reports.IRefReporter;
import org.theseed.dl4j.eval.reports.PatricRefGenomeComputer;
import org.theseed.dl4j.eval.reports.SingleRefGenomeComputer;
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
 * --ref		ID of a PATRIC genome to be used as the reference in a DEEP report
 * --home		home location of genome, to override the one in the GTO
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

    /** output directory */
    @Option(name = "-O", aliases = { "--outDir" }, metaVar = "outDir", usage = "output directory")
    private File outputDir;

    /** clear-output flag */
    @Option(name = "--clear", usage = "clear output directory before starting")
    private boolean clearOutputDir;

    /** file of reference genome GTO mappings */
    @Option(name = "--ref", usage = "ID of reference genome to use for a DEEP report")
    private String refGenomeID;

    /** overriding home location of genome */
    @Option(name = "--home", usage = "home location of genome (overrides GTO value)")
    private String homeName;


    @Override
    public void validateEvalParms() throws IOException {
        // Insure we can read the genome file, if there is one.
        if (this.inFile != null && ! this.inFile.canRead()) {
            throw new FileNotFoundException("Input " + this.inFile + " does not exist or cannot be read.");
        }
        // Set up the output directory.
        this.validateOutputDir(this.outputDir, this.clearOutputDir);
        // Set the reporting options for single-genome output.
        EvalReporter reporter = this.getReporter();
        reporter.setOption(EvalReporter.Option.NOSUMMARY);
        reporter.setOption(EvalReporter.Option.P3REPORT);
    }

    public void setDefaults() {
        this.inFile = null;
        this.outFile = null;
        this.outputDir = new File(System.getProperty("user.dir"));
        this.refGenomeID = null;
        this.homeName = null;
    }

    @Override
    public void runCommand() throws Exception {
        // Set up the reference-genome engine (if necessary).
        if (this.getReporter() instanceof IRefReporter) {
            IRefReporter refReporter = (IRefReporter) this.getReporter();
            if (refGenomeID != null)
                refReporter.setEngine(new SingleRefGenomeComputer(refGenomeID));
            else {
                // If we have a refGenomes.fa, we use the PATRIC-style reference genome computer.
                File modelDir = this.getModelDir();
                File refFasta = new File(modelDir, "refGenomes.fa");
                if (refFasta.canRead())
                    refReporter.setEngine(new PatricRefGenomeComputer(modelDir));
                else
                    log.info("No reference genomes available.");
            }
        }
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
        // Override the home, if required, and insure it is valid.
        if (this.homeName != null)
            genome.setHome(this.homeName);
        genome.checkHome();
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
        gReport.store(genome, this.getRoleDefinitions(), version, this.getOptions());
        if (this.outFile != null) {
            genome.save(this.outFile);
        } else {
            genome.save(System.out);
        }
        // Finish processing.
        this.close();
    }

}
