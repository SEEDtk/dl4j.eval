package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.eval.reports.EvalDeepReporter;
import org.theseed.dl4j.eval.reports.EvalReporter;
import org.theseed.dl4j.eval.reports.FileRefGenomeComputer;
import org.theseed.dl4j.eval.reports.IRefReporter;
import org.theseed.dl4j.eval.reports.PatricRefGenomeComputer;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Genome;
import org.theseed.utils.ParseFailureException;

/**
 *
 * This is a subclass that evaluates genomes.  The various subclasses get the actual genome data into memory so it can be processed here.
 * The key file structure that powers this object is the evaluation directory.  The evaluation directory must contain the following files.
 *
 * 	roles.in.subsystems		a role definition file, tab-delimited without headers, containing role IDs in the first column
 * 							and role names in the third
 *  roles.to.use			a list of the roles to use in consistency checking, in the order they appear in the consistency
 *  						models
 *  XXXXXX.ser				the occurrence-prediction model for the role with ID "XXXXXX" (in subdirectory "Roles")
 *  comp.tbl				a completeness definition file (as produced by the kmer.reps RoleProcessor class); each completeness
 *  						group consists of a tab-delimited header line containing (0) the genome ID, (1) the minimum similarity score,
 *  						(2) the group name, and (3) the seed protein sequence, followed by one or more data lines beginning
 *  						with three spaces and containing the ID of a universal role for the group, followed by a trailer
 *  						line containing "//"; the final group is "root" and matches everything.
 *
 * If comp.tbl is missing, no completeness/contamination check will be performed.
 *
 * @author Bruce Parrello
 *
 */
public abstract class Evaluator extends BaseEvaluator implements IConsistencyChecker {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(Evaluator.class);
    /** reporting object */
    private EvalReporter reporter;
    /** current output directory */
    private File outDir;
    /** TRUE if the reporter has been initialized, else FALSE */
    private boolean reportsOpened;

    // COMMAND LINE OPTIONS

    /** summary-only flag */
    @Option(name = "--terse", usage = "do not output detail files")
    private boolean terse;

    /** reporting format */
    @Option(name = "--format", usage = "format for output reports")
    private EvalReporter.Type format;

    /** sensitivity of protein comparison on reports */
    @Option(name = "-s", aliases = { "--sensitivity" }, metaVar = "0.9", usage = "maximum distance for a close protein")
    private double sensitivity;

    /**
     * Construct an evaluator.
     */
    public Evaluator() {
        super();
        // Initialize the output values.
        this.terse = false;
        this.reportsOpened = false;
        this.outDir = new File(System.getProperty("user.dir"));
        this.format = EvalReporter.Type.TEXT;
        this.sensitivity = 0.8;
    }

    /**
     * @return the detail level needed in genomes read from PATRIC (can be overridden by subclasses)
     */
    public P3Genome.Details getDetailLevel() {
        return this.getReporter().getDetailLevel();
    }

    /**
     * Validate the output-related parameters.
     *
     * @throws ParseFailureException
     */
    @Override
    protected void validateOutputParms() throws ParseFailureException {
        // Create the reporting object.
        this.reporter = EvalReporter.create(this.format);
        // Check for terse mode.
        if (this.terse)
            this.getReporter().setOption(EvalReporter.Option.NODETAILS);
        // Here we tune the deep report.
        if (this.getReporter() instanceof EvalDeepReporter) {
            // Validate the sensitivity.
            if (this.sensitivity < 0 || this.sensitivity >= 1.0)
                throw new ParseFailureException("Invalid sensitivity: must be between 0 and 1, exclusive.");
            // Set the special parameters.
            EvalDeepReporter deepReporter = ((EvalDeepReporter) this.getReporter());
            deepReporter.setSensitivity(this.sensitivity);
        }
    }

    /**
     * Set up the reference-genome computation engine.
     *
     * @param refGenomeFile		reference-genome file (or NULL if PATRIC is to be used)
     *
     * @throws IOException
     */
    protected void setupRefGenomeEngine(File refGenomeFile) throws IOException {
        if (this.getReporter() instanceof IRefReporter) {
            IRefReporter refReporter = (IRefReporter) this.getReporter();
            if (refGenomeFile != null) {
                refReporter.setEngine(new FileRefGenomeComputer(refGenomeFile));
            } else {
                refReporter.setEngine(new PatricRefGenomeComputer(this.getModelDir()));
            }
        }
    }

    /**
     * Validate the output directory and optionally clear it.  This also sets the directory,
     * so it must be called after "validateParms".
     *
     * @param outputDir			proposed output directory
     * @param clearOutputDir	TRUE if the directory should be cleared
     *
     * @throws IOException
     */
    protected void validateOutputDir(File outputDir, boolean clearOutputDir) throws IOException {
        // Check the output directory.
        if (! outputDir.exists()) {
            log.info("Creating directory {}.", outputDir);
            if (! outputDir.mkdirs()) {
                throw new IOException("Could not create output directory.");
            }
        } else if (! outputDir.isDirectory()) {
            throw new FileNotFoundException("Output directory " + outputDir + " is invalid.");
        } else if (clearOutputDir) {
            log.info("Erasing output directory.");
            FileUtils.cleanDirectory(outputDir);
        }
        this.setOutDir(outputDir);
        log.info("Output will be in directory {}.", this.outDir);
    }

    /**
     * Inform the reporting object of whether or not we have completeness data.
     *
     * @param exists	TRUE if we have completeness data
     */
    protected void setHaveCompleteness(boolean exists) {
        this.getReporter().setHaveCompleteness(exists);
    }

    /**
     * @return the list of roles in the roles-to-use file
     *
     * @param rolesToUseFile	file containing role IDs in the first column
     * @throws IOException
     */
    public static List<String> readRolesToUse(File rolesToUseFile) throws IOException {
        List<String> retVal = new ArrayList<String>(2800);
        try (TabbedLineReader rolesToUse = new TabbedLineReader(rolesToUseFile, 3)) {
            for (TabbedLineReader.Line line : rolesToUse) {
                retVal.add(line.get(0));
            }
        }
        return retVal;
    }

    /**
     * Reset for another evaluation.
     */
    protected void clearGenomes() {
        this.reportsOpened = false;
    }

    /**
     * Write the output from the evaluations.
     *
     * @throws IOException
     */
    protected void writeOutput() throws IOException {
        log.info("Writing output.");
        // If this is our first time, initialize the reports.
        if (! this.reportsOpened) {
            this.getReporter().open(this.getVersion(), this.getRoleDefinitions(), this.getModelDir());
            this.reportsOpened = true;
        }
        // Set up the reference genomes.
        this.getReporter().setupGenomes(this.getGReports());
        // Write the output.
        for (int g = 0; g < this.getGenomeCount(); g++) {
            GenomeStats gReport = this.getGReport(g);
            this.getReporter().writeGenome(gReport);
        }
    }

    // Finish processing and clean up.
    public void close() {
        // Close the report object.
        this.getReporter().close();
        // Output the time spent.
        log.info("{} genomes evaluated at {} seconds/genome.", this.getGenomesProcessed(), this.getSpeed());
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

    /**
     * Specify the output directory.
     *
     * @param outDir	new output directory
     */
    protected void setOutDir(File outDir) {
        this.outDir = outDir;
        this.getReporter().setOutDir(outDir);
    }

    /**
     * @return the reporter
     */
    public EvalReporter getReporter() {
        return reporter;
    }

}
