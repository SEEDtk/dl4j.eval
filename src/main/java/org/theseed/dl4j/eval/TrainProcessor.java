/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.text.TextStringBuilder;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.dl4j.train.ClassTrainingProcessor;
import org.theseed.dl4j.train.RunStats;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.ICommand;
import org.theseed.utils.Parms;

/**
 * This object trains classifiers for genome annotation consistency evaluation.  The input file has role IDs along the top and a genome ID
 * in the first column.  A full pass is made through the file to get a list of the roles and the counts found for each.  A model is then
 * trained for every role.  These models can be used to create a consistency profile for new genomes.
 *
 * The models are trained by org.dlj4j.run.ClassTrainingProcessor.  The positional parameters are the name of the parameter file, the name
 * of the model directory, and the size of the testing set.  The parameter file should not contain the --testSize, --col, --name, --input,
 * --meta, or --comment parameters, as these will be added by this method.  The training set must be in a file named "training.tbl" in the
 * model directory.
 *
 * The following command-line options are supported.
 *
 * -m	minimum acceptable accuracy for a model to be output; the default is 0.95
 *
 * @author Bruce Parrello
 *
 */
public class TrainProcessor implements ICommand {

    /** logging facility */
    private static Logger log = LoggerFactory.getLogger(ClassTrainingProcessor.class);

    // FIELDS

    /** constant argument list */
    private List<String> basicArgs;
    /** map of role IDs to maximum labels */
    CountMap<String> roleMap;
    /** training file */
    File trainingFile;
    /** array of training file labels */
    String[] labels;
    /** array of accuracy scores */
    double[] ratings;
    /** label file */
    File labelFile;
    /** role file subdirectory */
    File rolesDir;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** minimum accuracy */
    @Option(name = "--min", aliases = { "-m" }, metaVar = "0.93", usage = "minimum acceptable accuracy")
    private double minAcc;

    /** model directory */
    @Argument(index = 0, metaVar = "modelDir", usage = "model directory", required = true)
    private File modelDir;

    /** parameter file */
    @Argument(index = 1, metaVar = "parms.prm", usage = "parameter file", required = true)
    private File parmFile;

    /** testing set size */
    @Argument(index = 2, metaVar = "120", usage = "testing set size", required = true)
    private int testSize;


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
        this.help = false;
        this.minAcc = 0.95;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Read in the parms file.
                this.basicArgs = Parms.fromFile(this.parmFile);
                if (! this.modelDir.isDirectory()) {
                    throw new FileNotFoundException("Model directory " + this.modelDir + " not found or invalid.");
                } else {
                    // Compute the training file name.
                    this.trainingFile = new File(this.modelDir, "training.tbl");
                    if (! this.trainingFile.exists()) {
                        throw new FileNotFoundException("Model directory " + this.modelDir + " does not contain a training file.");
                    }
                    // Compute the label file name.
                    this.labelFile = new File(this.modelDir, "labels.txt");
                    // Add the testing set size to the argument list.
                    this.basicArgs.add("--testSize");
                    this.basicArgs.add(Integer.toString(this.testSize));
                    // Insure we have the Roles subdirectory.
                    this.rolesDir = new File(this.modelDir, "Roles");
                    if (! rolesDir.isDirectory()) {
                        log.info("Creating Roles directory in {}.", this.modelDir);
                        boolean ok = rolesDir.mkdir();
                        if (! ok)
                            throw new IOException("Error creating Roles directory.");
                    }
                    // We made it this far, we can run the application.
                    retVal = true;
                }
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

    public void run() {
        try {
            // Here we must analyze the training file.  For each role, we need the distribution of labels.
            this.roleMap = new CountMap<String>();
            try (TabbedLineReader trainingStream = new TabbedLineReader(this.trainingFile)) {
                log.info("Analyzing training data in {}.", this.trainingFile.toString());
                // Get the label array.
                this.labels = trainingStream.getLabels();
                // Save the metadata label from the first column.
                this.basicArgs.add("--meta");
                this.basicArgs.add(this.labels[0]);
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
            // Allocate the accuracy array.
            this.ratings = new double[this.labels.length];
            // Create the training processor.
            ClassTrainingProcessor processor = new ClassTrainingProcessor();
            // Suppress saving the model.  We will do that later if we like it.
            processor.setSearchMode();
            // Now loop through the roles, creating the models.
            for (int i = 1; i < this.labels.length; i++) {
                // Start with a copy of the argument list.
                ArrayList<String> theseParms = new ArrayList<String>(this.basicArgs);
                // Add the role as the label column.
                String role = this.labels[i];
                theseParms.add("--col");
                theseParms.add(role);
                // Get the maximum label for this role.
                int maxLabel = this.roleMap.getCount(role);
                // Add the comment.
                theseParms.add("--comment");
                theseParms.add(String.format("Role %d: %s, maximum count %d", i, role, maxLabel));
                // Finally, push on the model directory.
                theseParms.add(this.modelDir.getPath());
                // Now we need to create the label file.
                try (PrintWriter labelStream = new PrintWriter(this.labelFile)) {
                    for (int j = 0; j <= maxLabel; j++) {
                        labelStream.println(j);
                    }
                }
                // We are finally ready.  Create the model.
                String[] argBuffer = new String[theseParms.size()];
                theseParms.toArray(argBuffer);
                boolean ok = processor.parseCommand(argBuffer);
                if (ok) {
                    processor.run();
                    ratings[i] = processor.getRating();
                    if (ratings[i] >= this.minAcc) {
                        // Here we can keep this model.
                        processor.saveModelForced(new File(this.rolesDir, role + ".ser"));
                    }
                } else {
                    throw new RuntimeException("Error processing role.");
                }
            }
            // Now log all the scores.
            String boundary = StringUtils.repeat('=', 34);
            // We will build the report in here.
            TextStringBuilder buffer = new TextStringBuilder(33 * (this.labels.length + 4));
            buffer.appendNewLine();
            buffer.appendln(boundary);
            // Write out the heading line.
            buffer.appendln(String.format("%-21s %11s", "Role", "Score"));
            // Write out a space.
            buffer.appendln("");
            // Write out the data lines.
            for (int i = 1; i < this.labels.length; i++) {
                char flag = (this.ratings[i] < this.minAcc ? '*' : ' ');
                buffer.appendln(String.format("%-21s %11.4f%c", this.labels[i], this.ratings[i], flag));
            }
            // Write out a trailer line.
            buffer.appendln(boundary);
            // Write it to the output log.
            String report = buffer.toString();
            log.info(report);
            // Write it to the trial file.
            try {
                RunStats.writeTrialReport(processor.getTrialFile(), "Summary of Model Ratings", report);
            } catch (IOException e) {
                System.err.println("Error writing trial log:" + e.getMessage());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
