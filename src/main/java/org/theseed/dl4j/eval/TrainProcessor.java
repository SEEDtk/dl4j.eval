/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.counters.RegressionStatistics;
import org.theseed.dl4j.train.CrossValidateProcessor;
import org.theseed.io.ParmDescriptor;
import org.theseed.io.ParmFile;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This object trains classifiers for genome annotation consistency evaluation.  The input file has role IDs along the top and a genome ID
 * in the first column.  A full pass is made through the file to get a list of the roles and the counts found for each.  A model is then
 * trained for every role.  These models can be used to create a consistency profile for new genomes.
 *
 * The models are trained by org.dlj4j.run.RandomForestTrainingProcessor.  The positional parameters are the name of the parameter file
 * and the name of the model directory.  The training set must be in a file named "training.tbl" in the model directory.
 *
 * The following command-line options are supported.
 *
 * -m	minimum acceptable accuracy for a model to be output; the default is 0.93
 * -q	maximum acceptable IQR during cross-validation
 * -k	number of folds to use for cross-validation
 *
 * --resume		name of the output file from an interrupted run; the run will be resumed
 * --parms		parameter file name in model directory
 *
 * @author Bruce Parrello
 *
 */
public class TrainProcessor extends BaseProcessor {

    /** logging facility */
    private static Logger log = LoggerFactory.getLogger(TrainProcessor.class);

    // FIELDS

    /** map of role IDs to maximum labels */
    private CountMap<String> roleMap;
    /** training file */
    private File trainingFile;
    /** array of training file labels */
    private String[] labels;
    /** label file */
    private File labelFile;
    /** role file subdirectory */
    private File rolesDir;
    /** output stream */
    private PrintStream output;
    /** set of roles already processed */
    private Set<String> processedRoles;
    /** parameter file */
    private File parmFile;

    // COMMAND LINE

    /** minimum accuracy */
    @Option(name = "--min", aliases = { "-m" }, metaVar = "0.95", usage = "minimum acceptable accuracy")
    private double minAcc;

    /** maximum cross-validation IQR */
    @Option(name = "--iqr", aliases = { "-q" }, metaVar = "0.07", usage = "maximum acceptable cross-validation IQR")
    private double maxIqr;

    /** number of folds for cross-validation */
    @Option(name = "--foldK", aliases = { "-k" }, metaVar = "10", usage = "number of folds for cross-validation")
    private int foldK;

    /** old output file if we are resuming */
    @Option(name = "--resume", metaVar = "output.log", usage = "if we are resuming, the output file from the interrupted run")
    private File resumeFile;

    /** parameter file name */
    @Option(name = "--parms", metaVar = "alt.parms.prm", usage = "parameter file name")
    private String parmFileName;

    /** model directory */
    @Argument(index = 0, metaVar = "modelDir", usage = "model directory", required = true)
    private File modelDir;

    @Override
    public void setDefaults() {
        this.help = false;
        this.minAcc = 0.93;
        this.maxIqr = 0.05;
        this.foldK = 5;
        this.parmFileName = "parms.prm";
        this.resumeFile = null;
    }

    @Override
    public boolean validateParms() throws IOException, ParseFailureException {
        if (! this.modelDir.isDirectory()) {
            throw new FileNotFoundException("Model directory " + this.modelDir + " not found or invalid.");
        } else {
            // Compute the parm file name.
            this.parmFile = new File(this.modelDir, this.parmFileName);
            if (! this.parmFile.canRead())
                throw new FileNotFoundException("Parameter file " + this.parmFile + " not found or unreadable.");
            // Compute the training file name.
            this.trainingFile = new File(this.modelDir, "training.tbl");
            if (! this.trainingFile.exists()) {
                throw new FileNotFoundException("Model directory " + this.modelDir + " does not contain a training file.");
            }
            // Compute the label file name.
            this.labelFile = new File(this.modelDir, "labels.txt");
            // Insure we have the Roles subdirectory.
            this.rolesDir = new File(this.modelDir, "Roles");
            if (! rolesDir.isDirectory()) {
                log.info("Creating Roles directory in {}.", this.modelDir);
                boolean ok = rolesDir.mkdir();
                if (! ok)
                    throw new IOException("Error creating Roles directory.");
            } else if (this.resumeFile == null) {
                // Clean out any previous results.  We are starting fresh.
                log.info("Erasing roles directory {}.", this.rolesDir);
                FileUtils.cleanDirectory(this.rolesDir);
            }
            // Verify the number of folds.
            if (this.foldK < 2)
                throw new ParseFailureException("Invalid fold count.  Must be 2 or more.");
            // Check for a resume situation.
            if (this.resumeFile == null) {
                // Normal processing.  Put the log to the standard output.
                System.out.println("#\trole_id\taccuracy\tIQR\tgood");
                this.output = System.out;
                // Denote no roles were pre-processed.
                this.processedRoles = Collections.emptySet();
                // Clear the trial log.
                File trials = new File(this.modelDir, "trials.log");
                if (trials.exists())
                    FileUtils.forceDelete(trials);
            } else {
                // Resume processing.  Save the roles we've already seen.
                this.processedRoles = TabbedLineReader.readSet(this.resumeFile, "role_id");
                // Open the resume file for append-style output with autoflush.
                FileOutputStream outStream = new FileOutputStream(this.resumeFile, true);
                this.output = new PrintStream(outStream, true);
            }
        }
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        // Here we must analyze the training file.  For each role, we need the distribution of labels.
        this.roleMap = new CountMap<String>();
        try (TabbedLineReader trainingStream = new TabbedLineReader(this.trainingFile)) {
            log.info("Analyzing training data in {}.", this.trainingFile.toString());
            // Get the label array.
            this.labels = trainingStream.getLabels();
            // Now read the training data.
            for (TabbedLineReader.Line line : trainingStream) {
                for (int i = 1; i < trainingStream.size(); i++) {
                    String role = this.labels[i];
                    int current = this.roleMap.getCount(role);
                    int newValue = line.getInt(i);
                    if (newValue > current) this.roleMap.setCount(role, newValue);
                }
            }
        }
        log.info("{} roles found for processing.", this.roleMap.size());
        // Create the cross-validation processor.
        CrossValidateProcessor processor = new CrossValidateProcessor();
        // Create the parameters for it.
        String[] crossParms = new String[] { "--type", "DECISION", "--folds", Integer.toString(this.foldK),
                "--parms", this.parmFile.toString(), this.modelDir.toString() };
        // Load the parameter file.
        ParmFile parms = new ParmFile(this.parmFile);
        // Update the ID column.
        parms.add(new ParmDescriptor("id", this.labels[0], "ID column for accuracy reports"));
        parms.add(new ParmDescriptor("meta", this.labels[0], "meta-data columns in input"));
        // Now loop through the roles, creating the models.
        for (int i = 1; i < this.labels.length; i++) {
            // Get the role ID.
            String role = this.labels[i];
            if (this.processedRoles.contains(role)) {
                log.info("{} already processed.", role);
            } else {
                // Get the maximum label for this role.
                int maxLabel = this.roleMap.getCount(role);
                // Create the label file.
                try (PrintWriter labelStream = new PrintWriter(this.labelFile)) {
                    for (int j = 0; j <= maxLabel; j++) {
                        labelStream.println(j);
                    }
                }
                // Now we need to update the parameter file for this role.
                parms.add(new ParmDescriptor("col", role, "role being predicted"));
                // Add the model file name.
                File modelFile = new File(this.rolesDir, role + ".ser");
                parms.add(new ParmDescriptor("name", modelFile.toString(), "output model file"));
                // Add the comment.
                String comment = String.format("Role %d: %s, maximum count %d", i, role, maxLabel);
                parms.add(new ParmDescriptor("comment", comment, "Comment for trial log"));
                // Save the parm file.
                parms.save(this.parmFile);
                double rating = 0.0;
                double iqr = 1.0;
                // Create the final argument buffer.
                boolean ok = processor.parseCommand(crossParms);
                if (ok) {
                    processor.run();
                    RegressionStatistics stats = processor.getResults();
                    rating = 1.0 - stats.trimean();
                    iqr = stats.iqr();
                } else {
                    throw new RuntimeException("Error processing role.");
                }
                String goodFlag = "1";
                if (rating < this.minAcc || iqr > this.maxIqr) {
                    // Here we can't keep the model.
                    log.info("Model {} rejected due to insufficient accuracy (trimean = {}, iqr = {}.", modelFile, rating, iqr);
                    FileUtils.forceDelete(modelFile);
                    goodFlag = "";
                }
                this.output.format("%d\t%s\t%8.5f\t%8.5f\t%s%n", i, role, rating, iqr, goodFlag);
            }
        }
    }
}
