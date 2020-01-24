/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.preprocessor.DataNormalization;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.ICommand;

/**
 * This processor performs the actual consistency evaluations.  It will take as input a tabbed file of actual role counts and
 * use it to compute the predicted role counts.  The extent to which the predictions match the actual counts determines the
 * consistency.
 *
 * The positional parameters are the name of the model directory, the name of the output directory, and the name of the
 * input file.  The input file should have a genome ID in the first column and then one column per role in the order
 * expected by the models.  This order is given by the "roles.to.use" file in the model directory.  The model directory
 * should contain a <i>role</i><code>.ser</code> file for each role to be checked.  This may be a smaller set than the
 * total number of roles used as input.  The header of the input file should contain the role ID for each column.
 *
 * When we are done, the output directory will contain a I<genomeID><code>.out</code> file for each genome.  This file will
 * contain the coarse and fine consistency numbers plus the predicted and actual counts for each role.  In addition, a
 * summary of the results will be placed in "summary.tbl" in the same directory.
 *
 * The command-line options are as follows.
 *
 * -b	maximum number of genomes to process at one time; the default is 5000
 *
 * --terse	do not write the individual output files, only the summary
 *
 * @author Bruce Parrello
 *
 */
public class PredictProcessor implements ICommand {

    /** logging facility */
    private static Logger log = LoggerFactory.getLogger(PredictProcessor.class);

    // FIELDS
    /** array of headers from the input file, not including genome ID column */
    private String[] roles;
    /** array of genome IDs */
    private String[] genomes;
    /** number of genomes read in this batch */
    private int nGenomes;
    /** output matrix, first index is genome, second is role */
    private int[][] rolesPredicted;
    /** input matrix, first index is genome, second is role */
    private int[][] rolesActual;
    /** array of flags indicating which roles have output */
    private boolean[] rolesUsed;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    @Option(name="-b", aliases = { "--batchSize", "--batch" }, metaVar = "1000", usage = "input batch size")
    private int batchSize;

    /** summary-only flag */
    @Option(name = "--terse", usage = "do not output detail files")
    private boolean terse;

    /** model directory */
    @Argument(index = 0, metaVar = "modelDir", usage = "model directory", required = true)
    private File modelDir;

    /** output directory */
    @Argument(index = 1, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    /** input file name */
    @Argument(index = 2, metaVar = "inFile.tbl", usage = "input file name", required = true)
    private File inFile;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        this.batchSize = 5000;
        // Set the defaults.
        this.help = false;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else if (! this.modelDir.isDirectory()) {
                throw new FileNotFoundException("Model directory " + this.modelDir + " not found or invalid.");
            } else if (! this.inFile.canRead()) {
                throw new FileNotFoundException("Input file " + this.inFile + " not found or unreadable.");
            } else {
                if (! this.outDir.exists()) {
                    log.info("Creating directory {}.", this.outDir.getPath());
                    if (! this.outDir.mkdir()) {
                        throw new IOException("Could not create output directory.");
                    }
                } else if (! this.outDir.isDirectory()) {
                    throw new FileNotFoundException("Output directory " + this.outDir + " is invalid.");
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
            long start = System.currentTimeMillis();
            // Open the input file and get the role list.
            log.info("Opening input file {}.", this.inFile.getPath());
            try (TabbedLineReader inStream = new TabbedLineReader(this.inFile)) {
                // Extract the role headers.
                this.roles = Arrays.copyOfRange(inStream.getLabels(), 1, inStream.size());
                log.info("{} input roles found.", this.roles.length);
                // Allocate the input and output arrays.
                this.rolesActual = new int[this.batchSize][this.roles.length];
                this.rolesPredicted = new int[this.batchSize][this.roles.length];
                // Allocate the array of roles-used flags.   They default to false.
                this.rolesUsed = new boolean[this.roles.length];
                // Finally, create an array for the genome IDs.
                this.genomes = new String[this.batchSize];
                // We are now ready to begin.  Loop through the input batches.
                int batchCount = 0;
                int genomeCount = 0;
                while (inStream.hasNext()) {
                    batchCount++;
                    log.info("Reading input batch {}.", batchCount);
                    for (this.nGenomes = 0; this.nGenomes < this.batchSize && inStream.hasNext(); this.nGenomes++) {
                        // Get the next input line.
                        TabbedLineReader.Line line = inStream.next();
                        // Save the genome ID.
                        this.genomes[this.nGenomes] = line.get(0);
                        // Read the counts.
                        for (int i = 0; i < this.roles.length; i++)
                            rolesActual[this.nGenomes][i] = line.getInt(i+1);
                        genomeCount++;
                    }
                    log.info("{} lines read into batch.", this.nGenomes);
                    // Create the holding area for the model input.
                    INDArray features = Nd4j.zeros(this.nGenomes, this.roles.length - 1);
                    // Now we loop through the roles, computing predictions for each role.
                    for (int iRole = 0; iRole < this.roles.length; iRole++) {
                        String role = this.roles[iRole];
                        File modelFile = new File(this.modelDir, role + ".ser");
                        if (modelFile.exists()) {
                            this.rolesUsed[iRole] = true;
                            log.info("Processing role {}.", role);
                            // Create the input matrix.  It contains all the columns but the one for our target role.
                            for (int i = 0; i < this.nGenomes; i++) {
                                for (int j = 0; j < iRole; j++) {
                                    features.put(i, j, this.rolesActual[i][j]);
                                }
                                for (int j = iRole; j < this.roles.length - 1; j++) {
                                    features.put(i, j, this.rolesActual[i][j+1]);
                                }
                            }
                            // Read the model and get the normalizer.
                            MultiLayerNetwork model = ModelSerializer.restoreMultiLayerNetwork(modelFile, false);
                            DataNormalization normalizer = ModelSerializer.restoreNormalizerFromFile(modelFile);
                            // Normalize the inputs.
                            normalizer.transform(features);
                            // Compute the predictions for this role.
                            INDArray output = model.output(features);
                            // Convert the predictions from one-hots to numbers.
                            for (int i = 0; i < this.nGenomes; i++) {
                                int jBest = 0;
                                double vBest = output.getDouble(i, 0);
                                for (int j = 1; j < output.size(1); j++) {
                                    double v = output.getDouble(i, j);
                                    if (v > vBest) {
                                        vBest = v;
                                        jBest = j;
                                    }
                                }
                                this.rolesPredicted[i][iRole] = jBest;
                            }
                        }
                    }
                    // Write the summary file and the output files.
                    try (PrintWriter outStream = new PrintWriter(new File(this.outDir, "summary.tbl"))) {
                        log.info("Writing output files for {} genomes in batch {}.", this.nGenomes, batchCount);
                        outStream.println("Genome\tCoarse\tFine\tOver\tUnder");
                        for (int g = 0; g < this.nGenomes; g++) {
                            String genome = this.genomes[g];
                            File outFile = new File(this.outDir, genome + ".out");
                            PrintWriter genomeStream = null;
                            try {
                                if (! this.terse) {
                                    log.info("Writing {}.", outFile.getPath());
                                    genomeStream = new PrintWriter(outFile);
                                    genomeStream.println("Role\tactual\tpredicted");
                                }
                                int coarse = 0;
                                int fine = 0;
                                int total = 0;
                                int over = 0;
                                int under = 0;
                                for (int i = 0; i < this.roles.length; i++) {
                                    if (this.rolesUsed[i]) {
                                        String role = this.roles[i];
                                        int actual = this.rolesActual[g][i];
                                        int predicted = this.rolesPredicted[g][i];
                                        if (! this.terse)
                                            genomeStream.format("%s\t%d\t%d%n", role, actual, predicted);
                                        if (actual == predicted) {
                                            // Here the prediction is perfect.
                                            coarse++;
                                            fine++;
                                        } else if (actual < predicted) {
                                            // The role occurs too few times, but if it occurs at least once, we are coarse-consistent.
                                            under++;
                                            if (actual > 0) coarse++;
                                        } else {
                                            // The role occurs too many times, but if it is supposed to be there, we are coarse-consistent.
                                            over++;
                                            if (predicted > 0) coarse++;
                                        }
                                        total++;
                                    }
                                }
                                double coarsePct = (double) coarse * 100.0 / total;
                                double finePct = (double) fine * 100.0 / total;
                                if (! this.terse) {
                                    genomeStream.println();
                                    genomeStream.format("*\tCoarse Consistency\t%8.2f%n", coarsePct);
                                    genomeStream.format("*\tFine Consistency\t%8.2f%n", finePct);
                                    genomeStream.format("*\tOverrepresented\t%8d%n", over);
                                    genomeStream.format("*\tUnderrepresented\t%8d%n", under);
                                }
                                outStream.format("%s\t%8.2f\t%8.2f\t%8d\t%8d%n", genome, coarsePct, finePct, over, under);
                            } finally {
                                if (genomeStream != null)
                                    genomeStream.close();
                            }
                        }
                    }
                }
                String rate = String.format("%6.2f", (double) (System.currentTimeMillis() - start) / (genomeCount * 1000));
                log.info("{} batches processed. {} genomes evaluated. {} seconds/genome.", batchCount, genomeCount, rate);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
