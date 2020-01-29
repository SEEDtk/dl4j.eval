/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.utils.ICommand;


/**
 * This processor performs the evaluations on GTO files.  It will take as input an evaluation directory and an input directory of
 * GTO files.  The GTO files do not need to be fully-loaded.  Only the features and the genome ID are needed.  In addition,
 * the protein translation must be present for the seed protein.  Representative-genome technology will be used for
 * completeness checking, and only the length and ID is necessary for each contig.
 *
 * The positional parameters are the name of the evaluation directory and the name of the input directory.
 *
 * The input directory should contain a GTO file for each genome to evaluate.  Alternative, the second parameter can be the
 * name of an input GTO file, and then only that single file will be processed.
 *
 * When we are done, the output directory will contain a I<genomeID><code>.tsv</code> file for each genome.  This file will
 * contain the quality numbers plus the predicted and actual counts for each role.  In addition, a
 * summary of the results will be placed in "summary.tsv" in the same directory.
 *
 * The command-line options are as follows.
 *
 * -v	show more detailed progress messages
 * -o	the name of the output directory (defaults to the current directory)
 *
 * --terse		do not write the individual output files, only the summary
 * --update		store the quality information in the GTO and write it to the output directory
 * --format		specify the output format-- HTML, DEEP, or TEXT
 * --ref		ID of a reference genome to use for all the evaluation reports
 *
 * @author Bruce Parrello
 *
 */
public class EvalProcessor extends Evaluator implements ICommand {

    /** TRUE for single-gto mode */
    private boolean singleton;


    // COMMAND LINE

    /** update-GTO flag */
    @Option(name = "-u", aliases = { "--update" }, usage = "store results in the input GTO")
    private boolean update;

    /** input file/directory name */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory", required = true)
    private File inDir;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.update = false;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.isHelp()) {
                parser.printUsage(System.err);
            } else if (this.validateParms()) {
                // Check the input directory.
                if (! this.inDir.isDirectory()) {
                    if (this.inDir.canRead()) {
                        this.singleton = true;
                    } else {
                        throw new FileNotFoundException("Input " + this.inDir + " is neither a directory or a readable file.");
                    }
                }
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
            // Read in the role maps.
            initializeData();
            // Are we doing one genome or many?
            if (! this.singleton) {
                log.info("Genomes will be read from {}.", this.inDir);
                GenomeDirectory genomeDir = new GenomeDirectory(this.inDir);
                // We know the number of genomes, so we can allocate our arrays.
                this.allocateArrays(genomeDir.size());
                // Loop through the genomes.  Note we track the genome's index in genomeStats;
                int iGenome = 0;
                for (Genome genome : genomeDir) {
                    processGenome(iGenome, genome);
                    // Prepare for the next genome.
                    iGenome++;
                }
            } else {
                // Here we are processing only one genome.
                this.allocateArrays(1);
                // Read in the genome and process it.
                log.info("Reading genome from {}.", this.inDir);
                Genome genome = new Genome(this.inDir);
                processGenome(0, genome);
            }
            // Evaluate the consistency of the genomes.
            evaluateConsistency();
            // Write the results.
            writeOutput();
            // If we are updating GTOs, do it here.
            if (this.update) {
                // We need the version string.
                String version = this.getVersion();
                for (int g = 0; g < this.getGenomeCount(); g++) {
                    GenomeStats gReport = this.getReport(g);
                    Genome gObject = gReport.getGenome();
                    String gId = gObject.getId();
                    log.trace("Updating GTO for {}.", gId);
                    gReport.store(gObject.getJson(), this.roleDefinitions, version);
                    File outFile = new File(this.getOutDir(), gId + ".gto");
                    gObject.update(outFile);
                }
            }
            // Finish processing.
            this.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
