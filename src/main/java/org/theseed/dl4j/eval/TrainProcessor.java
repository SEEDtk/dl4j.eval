/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.time.Duration;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.nd4j.linalg.dataset.DataSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.dl4j.decision.RandomForest;
import org.theseed.dl4j.decision.RandomForest.Method;
import org.theseed.dl4j.decision.RandomForest.Parms;
import org.theseed.io.ParmDescriptor;
import org.theseed.io.ParmFile;
import org.theseed.io.TabbedLineReader;

/**
 * This object trains classifiers for genome annotation consistency evaluation.  The input file has role IDs along the top and a genome ID
 * in the first column.  A full pass is made through the file to get a list of the roles and the counts found for each.  A model is then
 * trained for every role.  These models can be used to create a consistency profile for new genomes.
 *
 * The models are trained by org.dlj4j.run.RandomForestTrainingProcessor.  The positional parametere
 * is the name of the model directory.  The training set must be in a file named "training.tbl" in the model directory.
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
    protected static Logger log = LoggerFactory.getLogger(TrainProcessor.class);

    // FIELDS

    /** training file */
    private File trainingFile;
    /** testing file */
    private File testingFile;
    /** roles-to-use file */
    private File roleListFile;
    /** role file subdirectory */
    private File rolesDir;
    /** output stream */
    private PrintStream output;
    /** set of roles already processed */
    private Set<String> processedRoles;
    /** parameter file */
    private File parmFile;
    /** list of metadata columns */
    private static final List<String> META_COLS = Arrays.asList("genome");
    /** number of roles saved */
    private int saveCount;

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
    @Argument(index = 0, metaVar = "evalDir", usage = "evaluation directory", required = true)
    private File evalDir;

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
        if (! this.evalDir.isDirectory()) {
            throw new FileNotFoundException("Model directory " + this.evalDir + " not found or invalid.");
        } else {
            // Compute the parm file name.
            this.parmFile = new File(this.evalDir, this.parmFileName);
            if (! this.parmFile.canRead())
                throw new FileNotFoundException("Parameter file " + this.parmFile + " not found or unreadable.");
            // Compute the training file name.
            this.trainingFile = new File(this.evalDir, "training.tbl");
            if (! this.trainingFile.exists())
                throw new FileNotFoundException("Model directory " + this.evalDir + " does not contain a training file.");
            // Compute the testing file name.
            this.testingFile = new File(this.evalDir, "testing.tbl");
            if (! this.testingFile.exists())
                throw new FileNotFoundException("Model directory " + this.evalDir + " does not contain a testing file.");
            // Finally, compute the roles-to-use file name/
            this.roleListFile = new File(this.evalDir, "roles.to.use");
            if (! this.roleListFile.canRead())
                throw new FileNotFoundException("Role list file " + this.roleListFile + " not found or unreadable.");
            // Insure we have the Roles subdirectory.
            this.rolesDir = new File(this.evalDir, "Roles");
            if (! rolesDir.isDirectory()) {
                log.info("Creating Roles directory in {}.", this.evalDir);
                boolean ok = rolesDir.mkdir();
                if (! ok)
                    throw new IOException("Error creating Roles directory.");
            } else if (this.resumeFile == null) {
                // Clean out any previous results.  We are starting fresh.
                log.info("Erasing roles directory {}.", this.rolesDir);
                FileUtils.cleanDirectory(this.rolesDir);
            }
            // Verify the validation parameters.
            if (this.minAcc > 1.0 || this.minAcc < 0)
                throw new ParseFailureException("Minimum accuracy must be between 0 and 1.");
            if (this.maxIqr > 1.0 || this.maxIqr < 0)
                throw new ParseFailureException("Maximum IQR must be between 0 and 1.");
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
                this.saveCount = 0;
            } else {
                // Resume processing.  Save the roles we've already seen.
                this.processedRoles = new HashSet<String>();
                try (TabbedLineReader reader = new TabbedLineReader(this.resumeFile)) {
                    int idCol = reader.findField("role_id");
                    int goodCol = reader.findField("good");
                    for (TabbedLineReader.Line line : reader) {
                        this.processedRoles.add(line.get(idCol));
                        if (line.getFlag(goodCol)) this.saveCount++;
                    }
                }
                // Open the resume file for append-style output with autoflush.
                FileOutputStream outStream = new FileOutputStream(this.resumeFile, true);
                this.output = new PrintStream(outStream, true);
            }
        }
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        // These are used for timing computations.
        long activeTime = 0;
        int processCount = 0;
        // Read in the two datasets.
        DataSet trainingSet = EvalUtilities.readDataSet(this.trainingFile, META_COLS);
        log.info("{} records in training set.", trainingSet.getFeatures().rows());
        DataSet testingSet = EvalUtilities.readDataSet(this.testingFile, META_COLS);
        log.info("{} records in testing set.", testingSet.getFeatures().rows());
        // Now we need to get the hyper-parameters from the parm file.
        RandomForest.Parms hParms = new RandomForest.Parms(trainingSet);
        this.parseParameters(hParms);
        // Loop through the roles.
        try (TabbedLineReader roleReader = new TabbedLineReader(this.roleListFile, 3)) {
            int roleI = 1;
            for (TabbedLineReader.Line line : roleReader) {
                String roleId = line.get(0);
                String labelString = line.get(2);
                if (this.processedRoles.contains(roleId))
                    log.info("Skipping role {}: {} (already processed).", roleI, roleId);
                else {
                    long start = System.currentTimeMillis();
                    log.info("Processing role {}: {} with occurrence counts {}.", roleI, roleId, labelString);
                    // Create the training and testing sets for this role.
                    DataSet roleTrainingSet = EvalUtilities.prepareData(trainingSet, roleI-1, labelString);
                    DataSet roleTestingSet = EvalUtilities.prepareData(testingSet, roleI-1, labelString);
                    log.info("Generating classifier.");
                    RandomForest classifier = new RandomForest(roleTrainingSet, hParms);
                    double accuracy = classifier.getAccuracy(roleTestingSet);
                    log.info("Accuracy for {} classifier is {}.", roleId, accuracy);
                    if (accuracy < this.minAcc)
                        output.format("%d\t%s\t%8.4f\t\t%n", roleI, roleId, accuracy);
                    else if (this.maxIqr == 1.0) {
                        // We have good accuracy and we are in fast mode.  Save the classifier.
                        File classFile = new File(this.rolesDir, roleId + ".ser");
                        this.saveCount++;
                        log.info("Saving classifier to {}. {} saved out of {} so far.",
                                classFile, this.saveCount, roleI);
                        classifier.save(classFile);
                        output.format("%d\t%s\t%8.4f\t\t1%n", roleI, roleId, accuracy);
                    } else {
                        // We have good accuracy but we must verify the stability of the data.
                        int kFold = this.foldK;
                        double iqr = EvalUtilities.crossValidate(hParms, roleTrainingSet, kFold);
                        if (iqr > this.maxIqr) {
                            log.info("Role {} rejected due to IQR = {}.", roleId, iqr);
                            output.format("%d\t%s\t%8.4f\t%8.4f\t%n", roleI, roleId, accuracy, iqr);
                        } else {
                            File classFile = new File(this.rolesDir, roleId + ".ser");
                            this.saveCount++;
                            log.info("Saving classifier to {}. IQR = {}. {} saved out of {} so far.",
                                    classFile, iqr, this.saveCount, roleI);
                            classifier.save(classFile);
                            output.format("%d\t%s\t%8.4f\t%8.4f\t1%n", roleI, roleId, accuracy, iqr);
                        }
                    }
                    activeTime += (System.currentTimeMillis() - start);
                    processCount++;
                    log.info("Average time per role = {}.", Duration.ofMillis(activeTime).dividedBy(processCount).toString());
                }
                // Increment the role counter.
                roleI++;
            }
        }
    }

    /**
     * Parse the parameter file and use it to override the specified hyper-parameters.
     *
     * @param hParms	hyper-parameters to override
     *
     * @throws IOException
     */
    private void parseParameters(Parms hParms) throws IOException {
        ParmFile parms = new ParmFile(this.parmFile);
        ParmDescriptor parm = parms.get("maxDepth");
        if (parm != null && ! parm.isCommented())
            hParms.setMaxDepth(Integer.valueOf(parm.getValue()));
        parm = parms.get("minSplit");
        if (parm != null && ! parm.isCommented())
            hParms.setLeafLimit(Integer.valueOf(parm.getValue()));
        parm = parms.get("maxFeatures");
        if (parm != null && ! parm.isCommented())
            hParms.setNumFeatures(Integer.valueOf(parm.getValue()));
        parm = parms.get("nEstimators");
        if (parm != null && ! parm.isCommented())
            hParms.setNumTrees(Integer.valueOf(parm.getValue()));
        parm = parms.get("sampleSize");
        if (parm != null && ! parm.isCommented())
            hParms.setNumExamples(Integer.valueOf(parm.getValue()));
        parm = parms.get("seed");
        if (parm != null && ! parm.isCommented())
            RandomForest.setSeed(Integer.valueOf(parm.getValue()));
        parm = parms.get("method");
        if (parm != null && ! parm.isCommented())
            hParms.setMethod(Method.valueOf(parm.getValue()));
    }

}
