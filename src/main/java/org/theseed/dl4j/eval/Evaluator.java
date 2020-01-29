package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import org.apache.commons.io.FileUtils;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.preprocessor.DataNormalization;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.LoggerContext;

/**
 *
 * This is a subclass that evaluates genomes.  The various subclasses get the actual genome data into memory so it can be processed here.
 * The key file structure that powers this object is the evaluation directory.  The evaluation directory must contain the following files.
 *
 * 	roles.in.subsystems		a role definition file, tab-delimited without headers, containing role IDs in the first column
 * 							and role names in the third
 *  roles.to.use			a list of the roles to use in consistency checking, in the order they appear in the consistency
 *  						models
 *  XXXXXX.ser				the occurrence-prediction model for the role with ID "XXXXXX"
 *  comp.tbl				a completeness definition file (as produced by the kmer.reps RoleProcessor class); each completeness
 *  						group consists of a tab-delimited header line containing (0) the genome ID, (1) the minimum similarity score,
 *  						(2) the group name, and (3) the seed protein sequence, followed by one or more data lines beginning
 *  						with three spaces and containing the ID of a universal role for the group, followed by a trailer
 *  						line containing "//"; the final group is "root" and matches everything.
 *
 * @author Bruce Parrello
 *
 */
public class Evaluator {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(EvalProcessor.class);
    /** array of roles to use */
    protected ArrayList<String> roles;
    /** array of genome statistics objects */
    private GenomeStats[] reports;
    /** number of genomes read in this batch */
    private int nGenomes;
    /** output matrix, first index is genome, second is role */
    private int[][] rolesActual;
    /** array of flags indicating which roles have output */
    private boolean[] rolesUsed;
    /** role definition map */
    protected RoleMap roleDefinitions;
    /** completeness database */
    private List<UniversalRoles> compList;
    /** start time in milliseconds */
    private long start;
    /** evaluation version string */
    private String version;
    /** total number of genomes evaluated */
    private int gCount;
    /** reporting object */
    private EvalReporter reporter;
    /** directory of consistency role models */
    private File roleDir;

    // COMMAND LINE OPTIONS

    /** help option */
    @Option(name = "-h", aliases = { "--help" }, help = true)
    private boolean help;

    /** summary-only flag */
    @Option(name = "--terse", usage = "do not output detail files")
    private boolean terse;

    /** clear-output flag */
    @Option(name = "--clear", usage = "clear output directory before starting")
    private boolean clearOutDir;

    /** debug-message flag */
    @Option(name = "-v", aliases = { "--verbose", "--debug" }, usage = "show more detailed progress messages")
    private boolean debug;

    /** reporting format */
    @Option(name = "--format", usage = "format for output reports")
    private EvalReporter.Type format;

    /** override reference genome */
    @Option(name = "--ref", usage = "reference genome for deep evaluation reports")
    private String refGenomeId;

    /** sensitivity of protein comparison on reports */
    @Option(name = "-s", aliases = { "--sensitivity" }, metaVar = "0.9", usage = "maximum distance for a close protein")
    private double sensitivity;

    /** output directory */
    @Option(name = "-O", aliases = { "--outDir" }, metaVar = "outDir", usage = "output directory")
    private File outDir;

    /** model directory */
    @Argument(index = 0, metaVar = "modelDir", usage = "model directory", required = true)
    private File modelDir;

    /**
     * Construct an evaluator.
     */
    public Evaluator() {
        // Initialize the flags and the output directory value.
        this.help = false;
        this.terse = false;
        this.debug = false;
        this.clearOutDir = false;
        this.outDir = new File(System.getProperty("user.dir"));
        this.format = EvalReporter.Type.TEXT;
        this.sensitivity = 0.8;
    }

    /**
     * @return the detail level needed in genomes read from PATRIC (can be overridden by subclasses)
     */
    public P3Genome.Details getDetailLevel() {
        return this.reporter.getDetailLevel();
    }

    /**
     * Validate the common command-line parameters.
     *
     * @return TRUE if we can run, else FALSE
     *
     * @throws IOException
     */
    protected boolean validateParms() throws IOException {
        boolean retVal = false;
        if (! this.modelDir.isDirectory()) {
            throw new FileNotFoundException("Model directory " + this.modelDir + " not found or invalid.");
        } else {
            log.info("Evaluation database is in directory {}.", this.modelDir);
            // Get the consistency role directory.
            this.roleDir = new File(this.modelDir, "Roles");
            if (! this.roleDir.isDirectory())
                throw new FileNotFoundException("Roles subdirectory not found in " + this.modelDir + ".");
            // Check the output directory.
            if (! this.outDir.exists()) {
                log.info("Creating directory {}.", this.outDir);
                if (! this.outDir.mkdir()) {
                    throw new IOException("Could not create output directory.");
                }
            } else if (! this.outDir.isDirectory()) {
                throw new FileNotFoundException("Output directory " + this.outDir + " is invalid.");
            } else if (this.clearOutDir) {
                log.info("Erasing output directory.");
                FileUtils.cleanDirectory(this.outDir);
            }
            log.info("Output will be in directory {}.", this.outDir);
            // Create the reporting object.
            this.reporter = EvalReporter.Create(this.outDir, this.format);
            // Check for terse mode.
            if (this.terse)
                this.reporter.setOption(EvalReporter.Option.NODETAILS);
            // Validate the sensitivity.
            if (this.sensitivity < 0 || this.sensitivity >= 1.0)
                throw new IllegalArgumentException("Invalid sensitivity: must be between 0 and 1, exclusive.");
            retVal = true;
        }
        return retVal;
    }

    /**'
     * Turn off the summary report.
     */
    public void suppressSummary() {
        this.reporter.setOption(EvalReporter.Option.NOSUMMARY);
    }

    /**
     * Save the genome count and allocate the arrays.
     *
     * @param nGenomes	number of genomes to process
     */
    protected void allocateArrays(int nGenomes) {
        this.nGenomes = nGenomes;
        log.info("{} genomes per batch.", this.nGenomes);
        this.rolesActual = new int[this.nGenomes][this.roles.size()];
        this.reports = new GenomeStats[this.nGenomes];
    }

    /**
     * Read in the support structures.
     *
     * @throws IOException
     */
    protected void initializeData() throws IOException {
        if (this.debug) {
            // To get more progress messages, we set the log level in logback.
            LoggerContext loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory();
            ch.qos.logback.classic.Logger logger = loggerContext.getLogger("org.theseed");
            logger.setLevel(Level.toLevel("TRACE"));
        }
        // Start the clock.
        this.start = System.currentTimeMillis();
        this.gCount = 0;
        // Compute the version.
        this.loadVersion();
        log.info("Evaluation database version is {}.", this.version);
        // Read in the consistency roles.
        this.roles = new ArrayList<String>(2800);
        File rolesToUseFile = new File(this.modelDir, "roles.to.use");
        log.info("Reading consistency roles from {}.", rolesToUseFile);
        try (TabbedLineReader rolesToUse = new TabbedLineReader(rolesToUseFile, 3)) {
            for (TabbedLineReader.Line line : rolesToUse) {
                this.roles.add(line.get(0));
            }
        }
        log.info("{} roles will be used for consistency check.", this.roles.size());
        // Read in the role definition file.
        File roleDefineFile = new File(this.modelDir, "roles.in.subsystems");
        log.info("Reading role definitions from {}.", roleDefineFile);
        this.roleDefinitions = RoleMap.load(roleDefineFile);
        // Read the universal role group definitions.
        File compFile = new File(this.modelDir, "comp.tbl");
        log.info("Reading universal roles from {}.", compFile);
        this.compList = UniversalRoles.Load(compFile);
        // Create the roles-used array for the consistency checker.
        this.rolesUsed = new boolean[this.roles.size()];
    }

    /**
     * Parse out and process all the data for a single genome.
     *
     * @param iGenome	index of this genome
     * @param genome	the Genome object to process
     */
    protected void processGenome(int iGenome, Genome genome) {
        this.gCount++;
        log.trace("Processing #{} {}: {}.", gCount, genome.getId(), genome.getName());
        // Create the reporting object.
        GenomeStats gReport = new GenomeStats(genome);
        this.reports[iGenome] = gReport;
        // This will contain the role counts.
        CountMap<String> roleCounts = new CountMap<String>();
        // Loop through the genome's features.
        for (Feature feat : genome.getPegs()) {
            // Count it as a peg.
            gReport.countPeg(feat);
            // Loop through the useful roles.
            for (Role role : feat.getUsefulRoles(roleDefinitions)) {
                if (role.getId().contentEquals("PhenTrnaSyntAlph") && feat.getProteinLength() > 0)
                    gReport.countSeed(feat.getProteinTranslation());
                roleCounts.count(role.getId());
            }
        }
        // Compute the contig metrics.
        gReport.computeMetrics(genome);
        // Find the universal role group for this genome.
        UniversalRoles uniGroup = UniversalRoles.findGroup(gReport.getSeed(), this.compList);
        gReport.setGroup(uniGroup.getName());
        // Compute completeness and contamination.
        for (String uniRole : uniGroup.getRoles()) {
            gReport.completeRole(uniRole, roleCounts.getCount(uniRole));
        }
        // Store the rolesActual counts.
        for (int i = 0; i < this.roles.size(); i++) {
            this.rolesActual[iGenome][i] = roleCounts.getCount(this.roles.get(i));
        }
    }

    /**
     * Evaluate the consistency of the genomes whose role information is currently in memory.
     *
     * @throws IOException
     * @throws FileNotFoundException
     */
    protected void evaluateConsistency() throws IOException, FileNotFoundException {
        int fWidth = this.roles.size() - 1;
        INDArray features = Nd4j.zeros(this.nGenomes, fWidth);
        // Now we loop through the roles, computing predictions for each role.
        log.info("Analyzing roles for {} genomes.", this.nGenomes);
        for (int iRole = 0; iRole < this.roles.size(); iRole++) {
            String role = this.roles.get(iRole);
            File modelFile = new File(this.roleDir, role + ".ser");
            if (modelFile.exists()) {
                this.rolesUsed[iRole] = true;
                log.trace("Processing role #{} {}.", iRole, role);
                // Create the input matrix.  It contains all the columns but the one for our target role.
                for (int i = 0; i < this.nGenomes; i++) {
                    for (int j = 0; j < iRole; j++) {
                        features.put(i, j, this.rolesActual[i][j]);
                    }
                    for (int j = iRole; j < fWidth; j++) {
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
                    this.reports[i].consistentRole(role, jBest, this.rolesActual[i][iRole]);
                }
            }
        }
    }

    /**
     * Write the output from the evaluations.
     *
     * @throws IOException
     */
    protected void writeOutput() throws IOException {
        log.info("Writing output for batch.");
        // If this is our first time, initialize the reports.
        if (this.gCount == this.nGenomes) {
            this.reporter.open(this.version, this.roleDefinitions, this.modelDir);
            // Here we tune the deep report.
            if (this.reporter instanceof EvalDeepReporter) {
                EvalDeepReporter deepReporter = ((EvalDeepReporter) this.reporter);
                if (this.refGenomeId != null)
                    deepReporter.setRefGenomeOverride(this.refGenomeId);
                deepReporter.setSensitivity(this.sensitivity);
            }
        }
        for (int g = 0; g < this.nGenomes; g++) {
            GenomeStats gReport = this.reports[g];
            this.reporter.writeGenome(gReport);
        }
    }

    // Finish processing and clean up.
    public void close() {
        // Close the report object.
        this.reporter.close();
        // Output the time spent.
        String rate = String.format("%6.2f", (double) (System.currentTimeMillis() - start) / (this.gCount * 1000));
        log.info("{} genomes evaluated at {} seconds/genome.", this.gCount, rate);
    }


    /**
     * @return the number of genomes in the current batch
     */
    protected int getGenomeCount() {
        return nGenomes;
    }

    /**
     * @return the version string for the evaluation database
     */
    protected String getVersion() {
        return this.version;
    }

    /**
     * Compute the version of this evaluation model
     */
    private void loadVersion() {
        File vFile = new File(this.modelDir, "VERSION");
        try (Scanner vScanner = new Scanner(vFile)) {
            this.version = vScanner.nextLine();
        } catch (IOException e) {
            // Cannot read version file, use the directory name.
            this.version = String.format("%s (%s)", this.modelDir.getName(), e.getMessage());
        }
    }

    /**
     * @param nGenomes 	specify the number of genomes in the current batch
     */
    protected void setnGenomes(int nGenomes) {
        this.nGenomes = nGenomes;
    }

    /**
     * @return the specified genome quality report
     *
     * @param iGenome	index of the desired genome
     */
    protected GenomeStats getReport(int iGenome) {
        return this.reports[iGenome];
    }

    /**
     * @return TRUE if the user wants help instead of running the program
     */
    protected boolean isHelp() {
        return help;
    }

    /**
     * @return the output directory
     */
    protected File getOutDir() {
        return outDir;
    }


}
