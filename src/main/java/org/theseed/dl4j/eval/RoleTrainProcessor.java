/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.TabbedDataSetReader;
import org.theseed.dl4j.decision.RandomForest;
import org.theseed.dl4j.train.ClassPredictError;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command builds role classifiers for the roles in a set of genomes.  The input directory
 * should contain a "data.tbl" file that has a record for each genome and a column of counts for
 * each role.   A subdirectory called "Roles" will be created to hold the classifiers. In addition,
 * a "roles.to.run" file will be created that identifies the roles for which we have sufficiently
 * accurate predictors and the roles for which there is variable input along with the predictive
 * capability of each role.  This information will be used to determine how to run the predictors.
 *
 * The missing roles in each genome will be output to a report in the file "missingRoles.tbl".
 *
 * The positional parameter is the name of the input directory, which must contain "data.tbl" and
 * which will receive the "roles.to.run" file, the "missingRoles.tbl" file and the "Roles" subdirectory.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -m	minimum acceptable accuracy for a model to be output (default is 0.93
 * -q	maximum acceptable IQR during cross-validation (default 0.05)
 * -k	number of folds to use for cross-validation (default 5)
 * -t	testing set size, as a fraction of the total number of genomes (default 0.1)
 * -r	role definition file (default is "roles.in.subsystems" in the parent of the input directory)
 *
 * --max		maxmimum acceptable role count (default 1)
 * --clear		if specified, the Roles subdirectory will be erased before processing and the "missingRole.tbl"
 * 				and "roles.to.run" files will be overwritten; otherwise, we will resume from the previous run
 *
 * @author Bruce Parrello
 *
 */
public class RoleTrainProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RoleTrainProcessor.class);
    /** data input file name */
    private File dataFile;
    /** roles subdirectory */
    private File rolesDir;
    /** role-occurrence map */
    private Map<String, int[]> rolesToUse;
    /** master dataset of role counts */
    private DataSet masterTable;
    /** master list of role names, in column order */
    private List<String> roleNames;
    /** set of roles that contain no information (only one count), encoded as column index numbers */
    private BitSet uselessRoles;
    /** paraneters for the classifiers */
    RandomForest.Parms hParms;
    /** training set for current role */
    private DataSet trainingSet;
    /** testing set for current role */
    private DataSet testingSet;
    /** number of rows to put in the testing set */
    private int nTestRows;
    /** map of testing set rows to genome IDs */
    private String[] testGenomes;
    /** map of training set rows to genome IDs */
    private String[] trainGenomes;
    /** number of possible role-count values */
    private int nLabels;
    /** number of input columns */
    private int nInputs;
    /** role definition map */
    private RoleMap roleMap;
    /** counter of missing roles */
    private int missingRoles;
    /** counter of predictable roles */
    private int predictable;
    /** counter of unstable roles */
    private int unstable;
    /** counter of useful roles */
    private int useful;

    /** set of previously-run roles */
    private Set<String> rolesRun;
    /** print writer for missing-roles report */
    private PrintWriter missingRolesStream;
    /** print writer for roles-to-run report */
    private PrintWriter rolesToRunStream;
    /** file name for missing-roles report */
    private File missingRolesFile;
    /** file name for roles-to-run report */
    private File rolesToRunFile;
    /** header for roles-to-run report */
    private static final String ROLE_REPORT_HEADER = "role_id\tuseful\taccuracy\tiqr\tcommon\tpredictable\trole_name";
    /** data line format for roles-to-run report */
    private static final String ROLE_REPORT_FORMAT = "%s\t%s\t%6.4f\t%6.4f\t%3d\t%s\t%s%n";
    /** header line for missing-roles report */
    private static final String MISSING_ROLE_HEADER = "genome\tmissing_role\tdescription";

    // COMMAND-LINE OPTIONS

    /** minimum accuracy */
    @Option(name = "--min", aliases = { "-m" }, metaVar = "0.95", usage = "minimum acceptable accuracy")
    private double minAcc;

    /** maximum cross-validation IQR */
    @Option(name = "--iqr", aliases = { "-q" }, metaVar = "0.07", usage = "maximum acceptable cross-validation IQR")
    private double maxIqr;

    /** number of folds for cross-validation */
    @Option(name = "--foldK", aliases = { "-k" }, metaVar = "10", usage = "number of folds for cross-validation")
    private int foldK;

    /** maximum acceptable role count */
    @Option(name= "--max", metaVar = "5", usage = "maximum acceptable role count in the input")
    private int maxRoles;

    /** TRUE if the roles directory should be cleared */
    @Option(name = "--clear", usage = "if specified, the Roles subdirectory will be cleared before processing")
    private boolean clearFlag;

    /** fraction of input set to hold out for testing */
    @Option(name = "--testSize", aliases = { "-t" }, usage = "fraction of input set to use for testing")
    private double testSize;

    /** role definition file */
    @Option(name = "--roleMap", aliases = { "-r" }, usage = "role definition file")
    private File roleFile;

    /** input directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory for this genome/role set")
    private File inDir;

    @Override
    protected void setDefaults() {
        this.minAcc = 0.93;
        this.maxIqr = 0.05;
        this.foldK = 5;
        this.maxRoles = 1;
        this.testSize = 0.1;
        this.roleFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Validate the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        this.dataFile = new File(this.inDir, "data.tbl");
        if (! this.dataFile.canRead())
            throw new FileNotFoundException("Input directory " + this.inDir + " does not have a readable data.tbl file.");
        // Check the role file.
        if (this.roleFile == null) {
            File parentDir = this.inDir.getAbsoluteFile().getParentFile();
            this.roleFile = new File(parentDir, "roles.in.subsystems");
        }
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " is not found or unreadable.");
        log.info("Loading role definitions from {}.", this.roleFile);
        this.roleMap = RoleMap.load(this.roleFile);
        // Verify the validation parameters.
        if (this.minAcc > 1.0 || this.minAcc < 0)
            throw new ParseFailureException("Minimum accuracy must be between 0 and 1.");
        if (this.maxIqr > 1.0 || this.maxIqr < 0)
            throw new ParseFailureException("Maximum IQR must be between 0 and 1.");
        // Verify the number of folds.
        if (this.foldK < 2)
            throw new ParseFailureException("Invalid fold count.  Must be 2 or more.");
        // Verify the max-count.
        if (this.maxRoles < 1)
            throw new ParseFailureException("Maximum role occurrence count must be positive.");
        this.nLabels = this.maxRoles + 1;
        // Verify the testing-set fraction.
        if (this.testSize <= 0.0 || this.testSize >= 1)
            throw new ParseFailureException("Testing-set fraction must be between 0 and 1.");
        // Define the output files.
        this.rolesToRunFile = new File(this.inDir, "roles.to.run");
        this.missingRolesFile = new File(this.inDir, "missingRoles.tbl");
        // Denote no roles have been run yet.
        this.missingRoles = 0;
        this.predictable = 0;
        this.unstable = 0;
        this.useful = 0;
        this.rolesRun = new HashSet<String>(3000);
        // Set up the roles subdirectory.  If the directory is new, or there is no roles-to-run file,
        // we start the output files from scratch.  Otherwise, we have to compute our previous progress.
        this.rolesDir = new File(this.inDir, "Roles");
        if (! this.rolesDir.exists()) {
            log.info("Creating Roles subdirectory.");
            FileUtils.forceMkdir(this.rolesDir);
            this.initOutput();
        } else if (this.clearFlag) {
            log.info("Erasing Roles subdirectory.");
            FileUtils.cleanDirectory(this.rolesDir);
            this.initOutput();
        } else if (! this.rolesToRunFile.exists())
            this.initOutput();
        else
            this.restartOutput();
        return true;
    }

    /**
     * Initialize the two output streams.  Here we are resuming a run.  The roles-to-run file is read to
     * determine what roles have already been run and the missing-roles file is rebuilt to contain only the
     * data for those roles.  We then append to roles-to-run.
     *
     * @throws IOException
     */
    private void restartOutput() throws IOException {
        log.info("Restarting job.");
        // Get all the missing-roles records sorted by role.
        var missingRoleMap = new TreeMap<String, List<String>>();
        if (this.missingRolesFile.exists()) {
            try (TabbedLineReader missedRoleStream = new TabbedLineReader(this.missingRolesFile)) {
                int roleCol = missedRoleStream.findField("missing_role");
                for (TabbedLineReader.Line line : missedRoleStream) {
                    String roleId = line.get(roleCol);
                    List<String> missingLines = missingRoleMap.computeIfAbsent(roleId, x -> new ArrayList<String>());
                    missingLines.add(line.getAll());
                }
            }
        }
        // Now, read roles-to-run to get all of the roles fully processed.  We know the file exists
        // because we wouldn't be here if it didn't.
        try (TabbedLineReader roleStream = new TabbedLineReader(rolesToRunFile)) {
            int roleCol = roleStream.findField("role_id");
            int usefulCol = roleStream.findField("useful");
            int predictCol = roleStream.findField("predictable");
            int iqrCol = roleStream.findField("iqr");
            // Open the missing-roles file for output.  We will rebuild it here.
            this.missingRolesStream = new PrintWriter(this.missingRolesFile);
            this.missingRolesStream.println(MISSING_ROLE_HEADER);
            // Now we loop through the file, filling rolesRun and the missing-roles file,
            // and updating the counters.
            for (TabbedLineReader.Line line : roleStream) {
                String role = line.get(roleCol);
                this.rolesRun.add(role);
                if (line.getDouble(iqrCol) > this.maxIqr) this.unstable++;
                if (line.getFlag(usefulCol)) this.useful++;
                if (line.getFlag(predictCol)) this.predictable++;
                // Output the missing roles.
                List<String> missingLines = missingRoleMap.get(role);
                if (missingLines != null) {
                    for (String mLine : missingLines) {
                        this.missingRolesStream.println(mLine);
                        this.missingRoles++;
                    }
                }
            }
        }
        // Insure the missing-roles file is up-to-date.
        this.missingRolesStream.flush();
        // Open the role-to-run report for appending.
        this.rolesToRunStream = new PrintWriter(new FileOutputStream(this.rolesToRunFile, true));
        log.info("Job restarted with {} already run.  {} missing roles, {} predictable, {} useful.",
                this.rolesRun.size(), this.missingRoles, this.predictable, this.useful);
    }

    /**
     * Initialize the two output streams.  Both are open for writing.
     *
     * @throws IOException
     */
    private void initOutput() throws IOException {
        // Open the two output writers for normal writing and write the headers.
        this.rolesToRunStream = new PrintWriter(this.rolesToRunFile);
        this.rolesToRunStream.println(ROLE_REPORT_HEADER);
        this.rolesToRunStream.flush();
        // Since now the role file exists, insure the header is present.
        this.rolesToRunStream.flush();
        this.missingRolesStream = new PrintWriter(this.missingRolesFile);
        this.missingRolesStream.println(MISSING_ROLE_HEADER);
        log.info("Initialized for new job.");
    }

    @Override
    protected void runCommand() throws Exception {
        try {
            // Read in the master dataset.
            this.setupDataTable();
            log.info("Processing individual roles.");
            final int nRoles = this.roleNames.size();
            // Set up the timing facility.
            int runCount = 0;
            long start = System.currentTimeMillis();
            // Loop through the roles.
            for (int i = 0; i < nRoles; i++) {
                String roleId = this.roleNames.get(i);
                if (! this.rolesRun.contains(roleId)) {
                    String roleName = this.roleMap.get(roleId);
                    if (roleName == null) roleName = "<unknown>";
                    this.rolesRun.add(roleId);
                    // Compute the role's most common value.
                    int[] counts = this.rolesToUse.get(roleId);
                    int bestI = IntStream.range(0, this.maxRoles + 1).reduce(0, (i1,i2) -> (counts[i1] < counts[i2] ? i2 : i1));
                    if (this.uselessRoles.get(i)) {
                        log.info("{} has count {} for all genomes.", roleId, bestI);
                        this.rolesToRunStream.format(ROLE_REPORT_FORMAT, roleId, "", 1.0, 0.0, bestI, "Y", roleName);
                        this.predictable++;
                    } else {
                        this.useful++;
                        // Split master dataset into training and testing, and remove useless roles.
                        this.buildDataSets(i);
                        log.info("Generating classifier for {}.", roleId);
                        RandomForest classifier = new RandomForest(this.trainingSet, this.hParms);
                        double accuracy = classifier.getAccuracy(this.testingSet);
                        log.info("Accuracy found was {}.", accuracy);
                        if (accuracy < this.minAcc) {
                            // Here the role is unpredictable.
                            this.rolesToRunStream.format(ROLE_REPORT_FORMAT, roleId, "Y", accuracy, 1.0, bestI, "", roleName);
                        } else {
                            runCount++;
                            // Use cross-validation to insure the data is stable.
                            log.info("Checking stability of {}.", roleId);
                            String predictFlag;
                            // Concatenate the training set to the testing set.
                            double iqr = EvalUtilities.crossValidate(this.hParms, this.trainingSet.copy(), this.foldK);
                            if (iqr > this.maxIqr) {
                                log.info("{} is unstable.  IQR = {}.", roleId, iqr);
                                this.unstable++;
                                predictFlag = "";
                            } else {
                                File classFile = new File(this.rolesDir, roleId + ".ser");
                               this.predictable++;
                                predictFlag = "Y";
                                // Save the classifier.
                                classifier.save(classFile);
                                log.info("{} with IQR {} and accuracy {} saved to {}.", roleId, iqr, accuracy, classFile);
                                // Now predict everything and pull out the ones where the predicted role count is lower
                                // than the actual.
                                this.checkForMissing(roleId, roleName, classifier, this.trainingSet,
                                        this.trainGenomes);
                                this.checkForMissing(roleId, roleName, classifier, this.testingSet,
                                        this.testGenomes);
                                this.missingRolesStream.flush();
                            }
                            this.rolesToRunStream.format(ROLE_REPORT_FORMAT, roleId, "Y", accuracy, iqr, bestI, predictFlag, roleName);
                            this.rolesToRunStream.flush();
                        }
                    }
                    if (log.isInfoEnabled() && runCount > 0) {
                        double rate = (System.currentTimeMillis() - start) / (runCount * 1000.0);
                        log.info("{} roles processed.  {} predictable, {} useful, {} unstable, {} missing role occurrences, {} seconds/role.",
                                this.rolesRun.size(), this.predictable, this.useful, this.unstable, this.missingRoles,
                                rate);
                    }
                }
            }
        } finally {
            // Insure the files close.
            this.missingRolesStream.close();
            this.rolesToRunStream.close();
        }
    }

    /**
     * Run the classifier and determine if any roles are missing.
     *
     * @param roleId		ID of role being predicted
     * @param roleName		name of the role for the reports
     * @param classifier	classifier for making predictions
     * @param mySet			dataset to run through predictor
     * @param myGenomeIds	array of genome IDs for the dataset rows
     */
    private void checkForMissing(String roleId, String roleName,
            RandomForest classifier, DataSet mySet, String[] myGenomeIds) {
        INDArray predictions = classifier.predict(mySet.getFeatures());
        INDArray actuals = mySet.getLabels();
        for (int g = 0; g < myGenomeIds.length; g++) {
            int actualCount = ClassPredictError.computeBest(actuals, g);
            int predictCount = ClassPredictError.computeBest(predictions, g);
            if (actualCount < predictCount) {
                log.info("Actual is {}, predicted is {} for {} missing in {}.", actualCount, predictCount,
                        roleId, myGenomeIds[g]);
                this.missingRolesStream.println(myGenomeIds[g] + "\t" + roleId + "\t" + roleName);
                this.missingRoles++;
            }
        }
    }

    /**
     * Build testing and training datasets for the indicated role.  The role is described by its index number.
     * This process is complicated.  First, we need to distribute the different role counts into both sets.
     * For now, we do this by insuring the testing set has one of each value and then randomizing the rest.
     * Second, we need to remove the columns for the useless roles and convert the indicated column to labels.
     *
     * @param iRole		index of the role column to use as the output label
     */
    protected void buildDataSets(int iRole) {
        // Our first task will be to determine the rows that go in the testing set.  We need to scramble
        // the rows (which we will do indirectly), and then put the first one we find of each output
        // value in the testing set.  First, the array with the scrambling.
        final int nRows = this.masterTable.numExamples();
        int[] shuffler = IntStream.range(0, nRows).toArray();
        ArrayUtils.shuffle(shuffler);
        // This boolean array tells us which values have a presence in the testing set.
        boolean[] found = new boolean[this.nLabels];
        // Finally, get the array of rows.
        INDArray rows = this.masterTable.getFeatures();
        // Now shuffle one of each of the desired values to the top.
        for (int i = 0; i < nRows; i++) {
            int r = shuffler[i];
            // Get the role count for this row.
            int count = rows.getInt(r, iRole);
            if (! found[count]) {
                // Here we have found the first instance of a value.  Swap it to the testing-set area.
                found[count] = true;
                if (count != i) {
                    shuffler[i] = shuffler[count];
                    shuffler[count] = r;
                }
            }
        }
        // Now we create the two datasets.  We also need directories for the genome IDs.
        this.testGenomes = new String[this.nTestRows];
        this.trainGenomes = new String[nRows - this.nTestRows];
        this.testingSet = this.buildDataSet(shuffler, iRole, this.testGenomes, 0, this.nTestRows);
        this.trainingSet = this.buildDataSet(shuffler, iRole, this.trainGenomes, this.nTestRows, nRows);
    }



    /**
     * This method creates a single dataset for classifying the specified role from the specified
     * range of the shuffle array that points into the master data table.
     *
     * @param shuffler		shuffle array containing the indices of the rows to use
     * @param iRole			column index of the output role to use for labels
     * @param gArray		array for storing the genome IDs of each output row
     * @param i1			position of the first index from the shuffe array to copy
     * @param i2			position past the last index from the suffle array to copy
     *
     * @return a dataset formed from the specified master table rows
     */
    protected DataSet buildDataSet(int[] shuffler, int iRole, String[] gArray, int i1, int i2) {
        // Get the original metadata array.
        List<String> masterGenomeIds = this.masterTable.getExampleMetaData(String.class);
        // Get the features of the master data table.
        INDArray features = this.masterTable.getFeatures();
        // Get the number of columns in the features array.
        int nCols = this.masterTable.numInputs();
        // Build the output arrays.
        final int outRows = i2 - i1;
        INDArray labelsOut = Nd4j.zeros(outRows, this.nLabels);
        INDArray featuresOut = Nd4j.zeros(outRows, this.nInputs);
        int outR = 0;
        // Loop through the shuffle array, extracting rows.
        for (int r = i1; r < i2; r++) {
            int inR = shuffler[r];
            gArray[r - i1] = masterGenomeIds.get(inR);
            // This will be the next output column.
            int outC = 0;
            for (int c = 0; c < nCols; c++) {
                if (c == iRole) {
                    // Here we have the output column.
                    int count = features.getInt(inR, c);
                    labelsOut.put(outR, count, 1.0);
                } else if (! this.uselessRoles.get(c)) {
                    // Here we have an output column.
                    double count = features.getDouble(inR, c);
                    featuresOut.put(outR, outC, count);
                    outC++;
                }
            }
            outR++;
        }
        DataSet retVal = new DataSet(featuresOut, labelsOut);
        return retVal;
    }

    /**
     * Read in the data files to get the master dataset and the list of role names.
     *
     * @throws IOException
     */
    protected void setupDataTable() throws IOException {
        // Survey the role counts.
        log.info("Computing role counts in {}.", this.dataFile);
        this.rolesToUse = EvalUtilities.buildRolesToUse(this.dataFile, this.maxRoles);
        // Initialize the reader with the first column as the only metadata column.
        try (TabbedDataSetReader dataReader = new TabbedDataSetReader(this.dataFile, List.of("1"))) {
            this.masterTable = dataReader.readAll();
            RandomForest.flattenDataSet(this.masterTable);
            this.roleNames = dataReader.getFeatureNames();
            int nRoles = this.roleNames.size();
            log.info("{} roles and {} genomes found in {}.", nRoles, this.masterTable.numExamples(),
                    this.dataFile);
            // Determine the index numbers for the columns that contain no information.
            // We count the number of non-zero occurrence counts.  If there is more than one, the
            // column is useful, otherwise it is useless.
            this.uselessRoles = new BitSet(this.roleNames.size());
            int uselessCount = 0;
            for (int i = 0; i < nRoles; i++) {
                int[] counts = this.rolesToUse.get(this.roleNames.get(i));
                long used = Arrays.stream(counts).filter(x -> x > 0).count();
                if (used <= 1) {
                    this.uselessRoles.set(i);
                    uselessCount++;
                }
            }
            log.info("{} roles are useless.", uselessCount);
            this.nInputs = nRoles - uselessCount;
            // Set up the hyper-parameters.
            this.hParms = new RandomForest.Parms(this.masterTable.numExamples(), this.nInputs);
            // Compute the testing-set size.
            this.nTestRows = (int) (this.masterTable.numExamples() * this.testSize);
            if (this.nTestRows <= this.maxRoles) this.nTestRows = this.nLabels;
        }
    }

    /**
     * @return the masterTable
     */
    protected DataSet getMasterTable() {
        return this.masterTable;
    }

    /**
     * @return the trainingSet
     */
    protected DataSet getTrainingSet() {
        return this.trainingSet;
    }

    /**
     * @return the testingSet
     */
    protected DataSet getTestingSet() {
        return this.testingSet;
    }

    /**
     * Specify the input directory.
     *
     * @param dir		input directory to use
     */
    protected void setInDir(File dir) {
        this.inDir = dir;
    }

    /**
     * @return the testGenomes array
     */
    protected String[] getTestGenomes() {
        return this.testGenomes;
    }

    /**
     * @return the trainGenomes array
     */
    protected String[] getTrainGenomes() {
        return this.trainGenomes;
    }

    /**
     * @return the roleNames
     */
    protected List<String> getRoleNames() {
        return this.roleNames;
    }

    /**
     * @return the uselessRoles
     */
    protected BitSet getUselessRoles() {
        return this.uselessRoles;
    }



}
