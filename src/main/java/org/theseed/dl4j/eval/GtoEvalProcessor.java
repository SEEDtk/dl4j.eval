/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.eval.reports.BinRefGenomeComputer;
import org.theseed.dl4j.eval.reports.EvalReporter;
import org.theseed.genome.Genome;
import org.theseed.utils.ICommand;

/**
 * This processor evaluates a single GTO.  This is less efficient than evaluating in a batch, but it
 * allows this program to be used in a genome-processing pipeline.  The GTO is read from the standard
 * input and written to the standard output, but both of these options can be overridden on the command
 * line.
 *
 * The output directory will contain the results of the evaluation, though they will also be stored in
 * the output GTO.  The main advantage here is that the output directory results can be viewed without
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
 * --improve	if specified, an attempt will be made to improve the genome by removing contigs
 * --bins		if specified, a bins.json file that can be used to compute coverage
 * --p3			force the HTML file name to GenomeReport.html
 *
 * @author Bruce Parrello
 *
 */
public class GtoEvalProcessor extends Evaluator implements ICommand {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GtoEvalProcessor.class);
    /** name of the subsystem projector file */
    private File improveFile;

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

    /** overriding home location of genome */
    @Option(name = "--home", usage = "home location of genome (overrides GTO value)")
    private String homeName;

    /** if specified, bad contigs will be removed if the genome is contaminated but mostly complete */
    @Option(name = "--improve", usage = "if specified, an attempt will be made to improve a bad genome")
    private boolean improveFlag;

    /** if specified, the location of a bins.json file that can be used to compute coverage and reference genomes */
    @Option(name = "--bins", usage = "if specified, the output (bins.json) file from a binning run used to create the genome")
    private File binFile;

    /** if specified, the HTML file name will be forced to GenomeReport.html */
    @Option(name = "--p3", usage = "if specified, the BV-BRC html output file name will be used")
    private boolean p3Flag;

    @Override
    public void setDefaults() {
        this.inFile = null;
        this.outFile = null;
        this.outputDir = new File(System.getProperty("user.dir"));
        this.homeName = null;
        this.improveFlag = false;
        this.binFile = null;
        this.p3Flag = false;
    }

    @Override
    public void validateEvalParms() throws IOException {
        // Insure we can read the genome file, if there is one.
        if (this.inFile != null && ! this.inFile.canRead()) {
            throw new FileNotFoundException("Input " + this.inFile + " does not exist or cannot be read.");
        }
        // Set up the output directory.
        this.validateOutputDir(this.outputDir, this.clearOutputDir);
        // If improvement is desired, verify we have a projector.
        this.improveFile = new File(this.getModelDir(), "projector.ser");
        if (this.improveFlag && ! this.improveFile.canRead())
            throw new FileNotFoundException("Improvement requested, but subsystem projector file {} is not found or unreadable.");
        // Validate the bins.json file.
        if (this.binFile != null && ! this.binFile.canRead())
            throw new FileNotFoundException("Binning analysis file " + this.binFile + " is not found or unreadable.");
        // Set the reporting options for single-genome output.
        EvalReporter reporter = this.getReporter();
        reporter.setOption(EvalReporter.Option.NOSUMMARY);
        if (this.p3Flag)
            reporter.setOption(EvalReporter.Option.P3REPORT);
    }

    @Override
    public void runCommand() throws Exception {
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
        // Create the genome analysis.
        var analyses = this.analyzeGenomes();
        var analysis = analyses[0];
        // Attempt to improve the genome.
        var gReport = this.getReport(0);
        boolean improved;
        if (! analysis.hasRefGenome()) {
            // Can't improve unless there is a reference genome.
            improved = false;
            if (this.improveFlag)
                log.info("Cannot improve genome:  no reference genome available.");
        } else if (this.improveFlag) {
            // Here the user wants us to try improvement.
            improved = this.improve(gReport, analysis, this.improveFile);
            if (improved) {
                log.info("Genome has improved.  Re-evaluating.");
                // We have changed the genome.  Re-evaluate it.
                processGenome(0, genome);
                evaluateConsistency();
                analyses = this.analyzeGenomes();
                // Refresh the genome report.
                gReport = this.getGReport(0);
            }
        }
        // If we are binning, store the coverage and the reference genome ID for the report.
        if (this.binFile != null) {
            double coverage = BinRefGenomeComputer.getCoverage(genome, this.binFile);
            BinRefGenomeComputer.storeBinData(genome, analysis.getRefGenomeId(), coverage);
        }
        // Write the results.
        writeOutput(analyses);
        // Now we need to write the GTO.  Get the version string.
        String version = this.getVersion();
        log.info("Writing evaluated genome.");
        // Store the evaluation report as the genome's quality information.
        gReport.store(genome, this.getRoleDefinitions(), version, this);
        if (this.outFile != null) {
            genome.save(this.outFile);
        } else {
            genome.save(System.out);
        }
        // Finish processing.
        this.close();
    }

}
