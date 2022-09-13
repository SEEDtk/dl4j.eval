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
import org.theseed.dl4j.eval.reports.EvalReporter;
import org.theseed.dl4j.eval.reports.FileRefGenomeComputer;
import org.theseed.dl4j.eval.reports.IRefReporter;
import org.theseed.dl4j.eval.reports.NullRefGenomeComputer;
import org.theseed.dl4j.eval.reports.PatricRefGenomeComputer;
import org.theseed.dl4j.eval.reports.RefGenomeComputer;
import org.theseed.dl4j.eval.reports.SingleRefGenomeComputer;
import org.theseed.dl4j.eval.stats.GenomeAnalysis;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
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
    /** reference-genome engine */
    private RefGenomeComputer refEngine;

    // COMMAND LINE OPTIONS

    /** summary-only flag */
    @Option(name = "--terse", usage = "do not output detail files")
    private boolean terse;

    /** file of reference genome GTO mappings */
    @Option(name = "--ref", usage = "file of taxon ID to reference-genome GTO mappings")
    private File refGenomeFile;

    /** reference genome ID to use */
    @Option(name = "--refId", usage = "ID of reference genome to use if all genomes should use the same one")
    private String refGenomeID;

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
        this.refGenomeFile = null;
        this.refGenomeID = null;
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
     * @throws IOException
     */
    @Override
    protected void validateOutputParms() throws ParseFailureException, IOException {
        // Create the reporting object.
        this.reporter = EvalReporter.create(this.format);
        // Check for terse mode.
        if (this.terse)
            this.getReporter().setOption(EvalReporter.Option.NODETAILS);
        // Validate the sensitivity.
        if (this.sensitivity < 0 || this.sensitivity >= 1.0)
            throw new ParseFailureException("Invalid sensitivity: must be between 0 and 1, exclusive.");
        // Set up the reference engine.
        if (this.reporter instanceof IRefReporter) {
            if (this.refGenomeID != null) {
                // Here, the user has specified and explicit reference genome for everything.
                this.refEngine = new SingleRefGenomeComputer(this.refGenomeID);
                log.info("Reference genome forced to {}.", this.refGenomeID);
            } else if (refGenomeFile != null) {
                // Here, we have a file telling us where to find the reference genomes for each taxonomic ID.
                this.refEngine = new FileRefGenomeComputer(this.refGenomeFile);
                log.info("Reference genomes will be computed using {}.", refGenomeFile);
            } else {
                // Here we are pulling reference genomes from PATRIC using the model's reference-genome file.
                this.refEngine = new PatricRefGenomeComputer(this.getModelDir());
                log.info("Reference genomes will be taken from PATRIC.");
            }
        } else {
            this.refEngine = new NullRefGenomeComputer();
            log.info("No reference genomes will be used.");
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
     * Get analysis of the genomes.
     */
    protected GenomeAnalysis[] analyzeGenomes() {
        // Get the genome reports.
        GenomeStats[] gReports = this.getGReports();
        // Allocate the output array.
        GenomeAnalysis[] retVal = new GenomeAnalysis[gReports.length];
        // Fill in the reference genomes.
        this.refEngine.setupReferences(gReports);
        // Loop through the genome reports, creating the analyses.
        for (int i = 0; i < gReports.length; i++) {
            GenomeStats gReport = gReports[i];
            Genome genome = gReport.getGenome();
            Genome refGenome = this.refEngine.ref(genome);
            log.info("Analyzing genome {}.", genome);
            retVal[i] = new GenomeAnalysis(gReport, refGenome, this.getRoleDefinitions(), this.sensitivity);
        }
        return retVal;
    }


    /**
     * Write the output from the evaluations.
     *
     * @throws IOException
     */
    protected void writeOutput(GenomeAnalysis[] analyses) throws IOException {
        log.info("Writing output.");
        // If this is our first time, initialize the reports.
        if (! this.reportsOpened) {
            this.getReporter().open(this.getVersion(), this.getRoleDefinitions(), this.getModelDir());
            this.reportsOpened = true;
        }
       // Write the output.
        for (int g = 0; g < this.getGenomeCount(); g++) {
            GenomeStats gReport = this.getGReport(g);
            this.getReporter().writeGenome(gReport, analyses[g]);
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

    /**
     * @return the protein sensitivity for deep analysis
     */
    public double getSensitivity() {
        return sensitivity;
    }

}
