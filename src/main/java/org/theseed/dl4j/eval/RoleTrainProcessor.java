/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.nd4j.linalg.dataset.DataSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.TabbedDataSetReader;
import org.theseed.dl4j.decision.RandomForest;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command builds role classifiers for the roles in a set of genomes.  The input directory
 * should contain a "data.tbl" file that has a record for each genome and a column of counts for
 * each role.   A subdirectory called "Roles" will be created to hold the classifiers. In addition,
 * a "roles.to.run" file will be created that identifies the roles for which we have sufficiently
 * accurate predictors and the roles for which there is variable input along with the predictive
 * capability of each role.  This information will be used to determine how to run the predictors.
 *
 * The missing roles in each genome will be output to a report on the standard output.
 *
 * The positional parameter is the name of the input directory, which must contain "data.tbl" and
 * which will receive the "roles.to.run" file and the "Roles" subdirectory.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for missing roles (if not STDOUT)
 * -m	minimum acceptable accuracy for a model to be output (default is 0.93
 * -q	maximum acceptable IQR during cross-validation (default 0.05)
 * -k	number of folds to use for cross-validation (default 5)
 *
 * --max	maxmimum acceptable role count (default 1)
 * --clear	if specified, the Roles subdirectory will be erased before processing
 *
 * @author Bruce Parrello
 *
 */
public class RoleTrainProcessor extends BaseReportProcessor {

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

    /** input directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory for this genome/role set")
    private File inDir;

    @Override
    protected void setReporterDefaults() {
        this.minAcc = 0.93;
        this.maxIqr = 0.05;
        this.foldK = 5;
        this.maxRoles = 1;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Validate the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        this.dataFile = new File(this.inDir, "data.tbl");
        if (! this.dataFile.canRead())
            throw new FileNotFoundException("Input directory " + this.inDir + " does not have a readable data.tbl file.");
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
        // Set up the roles subdirectory.
        this.rolesDir = new File(this.inDir, "Roles");
        if (! this.rolesDir.exists()) {
            log.info("Creating Roles subdirectory.");
            FileUtils.forceMkdir(this.rolesDir);
        } else if (this.clearFlag) {
            log.info("Erasing Roles subdirectory.");
            FileUtils.cleanDirectory(this.rolesDir);
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Survey the role counts.
        log.info("Computing role counts in {}.", this.dataFile);
        this.rolesToUse = EvalUtilities.buildRolesToUse(this.dataFile, this.maxRoles);
        // Read in the master dataset.
        this.setupDataTable();
        // Write the output header for the false-negative report.
        writer.println("genome\tmissing_role");
        // Set up the roles-to-run file.
        try (PrintWriter roleWriter = new PrintWriter(new File(this.inDir, "roles.to.run"))) {
            roleWriter.println("num\trole_id\tuseful\taccuracy\tiqr");
            // TODO process each role, producing report
        }
    }

    /**
     * Read in the data file to get the master dataset and the list of role names.
     *
     * @throws IOException
     */
    private void setupDataTable() throws IOException {
        // Initialize the reader with the first column as the only metadata column.
        try (TabbedDataSetReader dataReader = new TabbedDataSetReader(this.dataFile, List.of("1"))) {
            this.masterTable = dataReader.readAll();
            RandomForest.flattenDataSet(this.masterTable);
            this.roleNames = dataReader.getFeatureNames();
            log.info("{} roles  found in {}.", this.roleNames.size(), this.dataFile);
            // Determine the index numbers for the columns that contain no information.
            // We count the number of non-zero occurrence counts.  If there is more than one, the
            // column is useful, otherwise it is useless.
            this.uselessRoles = new BitSet(this.roleNames.size());
            final int n = this.masterTable.numInputs();
            int uselessCount = 0;
            for (int i = 0; i < n; i++) {
                int[] counts = this.rolesToUse.get(this.roleNames.get(i));
                long used = Arrays.stream(counts).filter(x -> x > 0).count();
                if (used <= 1) {
                    this.uselessRoles.set(i);
                    uselessCount++;
                }
            }
            log.info("{} roles are useless.", uselessCount);
        }
    }

}
