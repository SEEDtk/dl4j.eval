/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.dl4j.eval.reports.BinRefGenomeComputer;
import org.theseed.dl4j.eval.reports.IRefReporter;
import org.theseed.genome.BinFilter;
import org.theseed.genome.Genome;
import org.theseed.utils.ICommand;

/**
 * This class evaluates the binning of a metagenomic sample.  Each sample is stored in a directory, and the
 * individual genomes are in files named "binXX.gto", where "XX" is a number.  The evaluation results will
 * be stored in the "Eval" subdirectory of each sample's directory.
 *
 * The positional parameters are the name of the model directory and the name of the input sample directory.
 *
 * The command-line options are as follows.
 *
 * -v	show more detailed progress messages
 * -s	maximum distance for a close protein in deep reporting (0 to 1)
 * -R	in this case, the input directory is considered to be a directory of samples rather than a
 * 		single sample directory; all the samples will be processed
 *
 * --terse		do not write the individual output files, only the summary
 * --clear		clear each output directory before processing
 * --format		the format of the output reports
 * --ref		file of reference-genome mappings (the default is to use the reference genome of the bin)

 *
 * @author Bruce Parrello
 *
 */
public class BinProcessor extends Evaluator implements ICommand {

    // FIELDS
    /** list of input sample directories to process */
    private Collection<File> inputDirs;
    /** saved command-line parameters */
    private String[] options;

    // COMMAND-LINE OPTIONS

    /** directory-recursion flag */
    @Option(name = "-R", aliases = { "--recursive", "--sub" }, usage="process all subdirectories")
    private boolean recursive;

    /** input directory */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory for sample (or samples)", required = true)
    private File inDir;

    /** clear output directories */
    @Option(name = "--clear", usage = "clear the output subdirectory for each sample before processing")
    private boolean clearDir;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.recursive = false;
        this.clearDir = false;
        // Save the command-line options.
        this.options = args;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.isHelp()) {
                parser.printUsage(System.err);
            } else if (this.validateParms()) {
                if (! this.inDir.isDirectory()) {
                    throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
                } else if (! this.recursive) {
                    // Here we have a single input sample.  Insure it has bins.
                    File testBin = new File(this.inDir, "bin1.gto");
                    if (! testBin.exists()) {
                        throw new FileNotFoundException("Input directory " + this.inDir + " does not contain completed bins.");
                    }
                    // Set up to process the one input directory.
                    this.inputDirs = Collections.singleton(this.inDir);
                } else {
                    // Here we are recursing through the subdirectories of the input directory.
                    this.inputDirs = new ArrayList<File>();
                    for (File subFile : this.inDir.listFiles(File::isDirectory)) {
                        // Insure this directory is completed.
                        if (new File(subFile, "bins.report.txt").canRead() && new File(subFile, "bin1.gto").exists()) {
                            this.inputDirs.add(subFile);
                        } else {
                            log.debug("Skipping subdirectory {}:  no bins or not completed.", subFile);
                        }
                    }
                    log.info("{} subdirectories will be processed.", this.inputDirs.size());
                }
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
            // Connect the reference-genome computation engine.
            BinRefGenomeComputer binRefEngine = new BinRefGenomeComputer();
            if (this.getReporter() instanceof IRefReporter) {
                IRefReporter refReporter = (IRefReporter) this.getReporter();
                refReporter.setEngine(binRefEngine);
            }
            // Create a file-name filter for bin GTO files.
            FilenameFilter binFilter = new BinFilter();
            // Read in the role maps.
            initializeData();
            // Get the evaluator version.
            String version = this.getVersion();
            // Loop through the input directories.
            for (File sampleDir : this.inputDirs) {
                log.info("Processing sample in {}.", sampleDir);
                binRefEngine.setSampleDir(sampleDir);
                this.clearGenomes();
                // Create the output directory.
                File outDir = new File(sampleDir, "Eval");
                if (! outDir.exists()) {
                    log.info("Creating output directory {}.", outDir);
                    if (! outDir.mkdir())
                        throw new IOException("Error creating output directory for " + outDir);
                } else if (this.clearDir) {
                    log.info("Erasing output directory {}.", outDir);
                    FileUtils.cleanDirectory(outDir);
                }
                this.setOutDir(outDir);
                // Find all the bins in the input directory.
                log.info("Reading bins from {}.", sampleDir);
                File[] binFiles = sampleDir.listFiles(binFilter);
                log.info("{} bins found in {}.", binFiles.length, sampleDir);
                // Allocate the arrays and store the genomes.
                this.allocateArrays(binFiles.length);
                for (int i = 0; i < binFiles.length; i++) {
                    Genome bin = new Genome(binFiles[i]);
                    this.processGenome(i, bin);
                }
                // Evaluate the genomes and write the output.
                this.evaluateConsistency();
                this.writeOutput();
                // Update the GTOs with the quality information and create the index.tbl file.
                try (PrintWriter indexStream = new PrintWriter(new File(outDir, "index.tbl"))) {
                    indexStream.println("Bin ID\tBin Name\tGood");
                    for (int i = 0; i < binFiles.length; i++) {
                        GenomeStats gReport = this.getReport(i);
                        Genome gObject = gReport.getGenome();
                        String gId = gObject.getId();
                        log.debug("Updating GTO for {}.", gId);
                        gReport.store(gObject, this.getRoleDefinitions(), version, options);
                        gObject.update(binFiles[i]);
                        indexStream.format("%s\t%s\t%d%n",
                                gId, gObject.getName(), (gReport.isGood() ? 1 : 0));
                    }
                    indexStream.flush();
                }
                // Finish off the reporting.
                this.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
