package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import org.kohsuke.args4j.Argument;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.dl4j.decision.RandomForest;
import org.theseed.dl4j.eval.reports.RefGenomeComputer;
import org.theseed.dl4j.eval.stats.GenomeAnalysis;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.AnalysisEvent;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.subsystems.SubsystemProjector;
import org.theseed.utils.BaseProcessor;
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
public abstract class BaseEvaluator extends BaseProcessor implements IConsistencyChecker {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BaseEvaluator.class);
    /** array of roles to use */
    protected List<String> roles;
    /** array of genome statistics objects */
    private GenomeStats[] reports;
    /** number of genomes read in this batch */
    private int nGenomes;
    /** output matrix, first index is genome, second is role */
    private int[][] rolesActual;
    /** role definition map */
    private RoleMap roleDefinitions;
    /** completeness database */
    private List<UniversalRoles> compList;
    /** start time in milliseconds */
    private long start;
    /** evaluation version string */
    private String version;
    /** total number of genomes evaluated */
    private int gCount;
    /** directory of consistency role models */
    private File roleDir;
    /** array of consistency role models */
    private RandomForest[] models;
    /** sensitivity level for protein comparisons (maximum allowable distance) */
    private double sensitivity;
    /** reference-genome engine */
    private RefGenomeComputer refEngine;

    // COMMAND LINE OPTIONS

    /** model directory */
    @Argument(index = 0, metaVar = "modelDir", usage = "model directory", required = true)
    private File modelDir;

    /**
     * Construct an evaluator.
     */
    public BaseEvaluator() {
        // Initialize the flags.
        this.compList = null;
        this.sensitivity = 0.8;
    }

    /**
     * @return the detail level needed in genomes read from PATRIC (can be overridden by subclasses)
     */
    public abstract P3Genome.Details getDetailLevel();

    /**
     * Validate the common command-line parameters.
     *
     * @return TRUE if we can run, else FALSE
     *
     * @throws IOException
     */
    @Override
    protected final boolean validateParms() throws IOException, ParseFailureException {
        validateModelParms();
        validateOutputParms();
        // Validate the subclass parms.
        this.validateEvalParms();
        // Denote we can run.
        return true;
    }

    /**
     * Validate the output-related parameters.
     *
     * @throws ParseFailureException
     */
    protected abstract void validateOutputParms() throws ParseFailureException, IOException;

    /**
     * Validate the evaluation parameters not related to reporting.
     *
     * @throws FileNotFoundException
     */
    private void validateModelParms() throws FileNotFoundException {
        if (! this.modelDir.isDirectory())
            throw new FileNotFoundException("Model directory " + this.modelDir + " not found or invalid.");
        log.info("Evaluation database is in directory {}.", this.modelDir);
        // Get the consistency role directory.
        this.roleDir = new File(this.modelDir, "Roles");
        if (! this.roleDir.isDirectory())
            throw new FileNotFoundException("Roles subdirectory not found in " + this.modelDir + ".");
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
        // Start the clock.
        this.start = System.currentTimeMillis();
        this.gCount = 0;
        // Compute the version.
        this.loadVersion();
        log.info("Evaluation database version is {}.", this.version);
        // Read in the consistency roles.
        File rolesToUseFile = new File(this.modelDir, "roles.to.use");
        log.info("Reading consistency roles from {}.", rolesToUseFile);
        this.roles = readRolesToUse(rolesToUseFile);
        log.info("{} roles will be used for consistency check.", this.roles.size());
        // Create the role model array.
        this.models = new RandomForest[this.roles.size()];
        int rCount = 0;
        for (int i = 0; i < this.models.length; i++) {
            String role = this.roles.get(i);
            File modelFile = new File(this.roleDir, role + ".ser");
            if (! modelFile.exists())
                this.models[i] = null;
            else {
                this.models[i] = RandomForest.load(modelFile);
                rCount++;
            }
        }
        log.info("{} consistency models read from {}.", rCount, this.roleDir);
        // Read in the role definition file.
        File roleDefineFile = new File(this.modelDir, "roles.in.subsystems");
        log.info("Reading role definitions from {}.", roleDefineFile);
        this.roleDefinitions = RoleMap.load(roleDefineFile);
        // Read the universal role group definitions.
        File compFile = new File(this.modelDir, "comp.tbl");
        if (! compFile.exists()) {
            log.warn("No completeness checks will be performed.");
            this.setHaveCompleteness(false);
        } else {
            log.info("Reading universal roles from {}.", compFile);
            this.compList = UniversalRoles.Load(compFile);
            this.setHaveCompleteness(true);
        }
    }

    /**
     * @return the processing speed
     */
    public Duration getSpeed() {
        Duration retVal = Duration.ZERO;
        if (gCount > 0)
            retVal = Duration.ofMillis(System.currentTimeMillis() - start).dividedBy(gCount);
        return retVal;
    }

    /**
     * Inform the reporting software of whether or not we have completeness data.
     *
     * @param exists	TRUE if we have completeness data
     */
    protected abstract void setHaveCompleteness(boolean exists);

    /**
     * Parse out and process all the data for a single genome.
     *
     * @param iGenome	index of this genome
     * @param genome	the Genome object to process
     */
    protected void processGenome(int iGenome, Genome genome) {
        this.gCount++;
        log.debug("Processing #{} {}: {}.", gCount, genome.getId(), genome.getName());
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
        if (this.compList != null) {
            UniversalRoles uniGroup = UniversalRoles.findGroup(gReport.getSeed(), this.compList);
            gReport.setGroup(uniGroup.getName());
            // Compute completeness and contamination.
            for (String uniRole : uniGroup.getRoles()) {
                gReport.completeRole(uniRole, roleCounts.getCount(uniRole));
            }
        }
        // Store the SSU rRNA flag.
        gReport.setHasSsuRRna(! genome.getSsuRRna().isEmpty());
        // Store the rolesActual counts.
        for (int i = 0; i < this.roles.size(); i++) {
            this.rolesActual[iGenome][i] = roleCounts.getCount(this.roles.get(i));
        }
    }

    /**
     * Evaluate the consistency of the genomes whose role information is currently in memory.
     *
     * @throws IOException
     */
    protected void evaluateConsistency() throws IOException {
        int fWidth = roles.size() - 1;
        INDArray features = Nd4j.zeros(nGenomes, fWidth);
        // Now we loop through the roles, computing predictions for each role.
        log.info("Analyzing roles for {} genomes.", nGenomes);
        for (int iRole = 0; iRole < roles.size(); iRole++) {
            String role = roles.get(iRole);
            // Fetch the model.
            RandomForest model = this.models[iRole];
            if (model != null) {
                log.debug("Processing role #{} {}.", iRole, role);
                // Create the input matrix.  It contains all the columns but the one for our target role.
                for (int i = 0; i < nGenomes; i++) {
                    for (int j = 0; j < iRole; j++) {
                        features.put(i, j, rolesActual[i][j]);
                    }
                    for (int j = iRole; j < fWidth; j++) {
                        features.put(i, j, rolesActual[i][j+1]);
                    }
                }
                // Compute the predictions for this role.
                INDArray output = model.predict(features);
                // Convert the predictions from one-hots to numbers.
                for (int i = 0; i < nGenomes; i++) {
                    int jBest = 0;
                    double vBest = output.getDouble(i, 0);
                    for (int j = 1; j < output.size(1); j++) {
                        double v = output.getDouble(i, j);
                        if (v > vBest) {
                            vBest = v;
                            jBest = j;
                        }
                    }
                    this.storeActual(iRole, role, i, jBest);
                }
            }
        }
    }

    /**
     * Evaluate the consistency of the genomes whose role information is currently in memory.
     *
     * @param checker		consistency-checker object requesting the analysis
     * @param roles			list of role IDs
     * @param rolesActual	matrix of role counts in the genomes, first index is genome, second is role
     * @param nGenomes		number of genomes in the matrix
     *
     * @return an array of flags indicating which roles were used in the evaluation
     *
     * @throws IOException
     */
    public static boolean[] runConsistency(IConsistencyChecker checker, List<String> roles, int[][] rolesActual,
            int nGenomes) throws IOException {
        boolean[] retVal = new boolean[roles.size()];
        File roleDir = checker.getRoleDir();
        int fWidth = roles.size() - 1;
        INDArray features = Nd4j.zeros(nGenomes, fWidth);
        // Now we loop through the roles, computing predictions for each role.
        log.info("Analyzing roles for {} genomes.", nGenomes);
        for (int iRole = 0; iRole < roles.size(); iRole++) {
            String role = roles.get(iRole);
            File modelFile = new File(roleDir, role + ".ser");
            if (modelFile.exists()) {
                retVal[iRole] = true;
                log.debug("Processing role #{} {}.", iRole, role);
                // Create the input matrix.  It contains all the columns but the one for our target role.
                for (int i = 0; i < nGenomes; i++) {
                    for (int j = 0; j < iRole; j++) {
                        features.put(i, j, rolesActual[i][j]);
                    }
                    for (int j = iRole; j < fWidth; j++) {
                        features.put(i, j, rolesActual[i][j+1]);
                    }
                }
                // Read the model.
                RandomForest model = RandomForest.load(modelFile);
                // Compute the predictions for this role.
                INDArray output = model.predict(features);
                // Convert the predictions from one-hots to numbers.
                for (int i = 0; i < nGenomes; i++) {
                    int jBest = 0;
                    double vBest = output.getDouble(i, 0);
                    for (int j = 1; j < output.size(1); j++) {
                        double v = output.getDouble(i, j);
                        if (v > vBest) {
                            vBest = v;
                            jBest = j;
                        }
                    }
                    checker.storeActual(iRole, role, i, jBest);
                }
            }
        }
        return retVal;
    }

    /**
     * Store the actual role count for a role.
     *
     * @param iRole		index of the role being counted
     * @param role		ID of the role being counted
     * @param iGenome	index of the relevant genome
     * @param count		number of role occurrences
     */
    @Override
    public void storeActual(int iRole, String role, int iGenome, int count) {
        this.reports[iGenome].consistentRole(role, count, this.rolesActual[iGenome][iRole]);
    }

    /**
     * @return the array of genome report objects
     */
    public GenomeStats[] getGReports() {
        return this.reports;
    }

    /**
     * @return the genome report object for the genome in the specified row.
     */
    public GenomeStats getGReport(int idx) {
        return this.reports[idx];
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
            this.version = String.format("%s (%s)", this.modelDir.getName(), e.toString());
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
     * @return the role definition map
     */
    protected RoleMap getRoleDefinitions() {
        return roleDefinitions;
    }

    /**
     * @return the number of genomes processed
     */
    public int getGenomesProcessed() {
        return this.gCount;
    }

    /**
     * Validate the subclass parameters.
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    public abstract void validateEvalParms() throws IOException, ParseFailureException;

    /**
     * @return the evaluation directory
     */
    protected File getModelDir() {
        return this.modelDir;
    }

    /**
     * @return the role-classifier model directory
     */
    public File getRoleDir() {
        return this.roleDir;
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
     * Get analysis of the genomes.
     */
    protected GenomeAnalysis[] analyzeGenomes() {
        // Get the genome reports.
        GenomeStats[] gReports = this.getGReports();
        // Allocate the output array.
        GenomeAnalysis[] retVal = new GenomeAnalysis[gReports.length];
        // Fill in the reference genomes.
        this.getRefEngine().setupReferences(gReports);
        // Loop through the genome reports, creating the analyses.
        for (int i = 0; i < gReports.length; i++) {
            GenomeStats gReport = gReports[i];
            Genome genome = gReport.getGenome();
            Genome refGenome = this.getRefEngine().ref(genome);
            log.info("Analyzing genome {}.", genome);
            retVal[i] = new GenomeAnalysis(gReport, refGenome, this.getRoleDefinitions(), this.sensitivity);
        }
        return retVal;
    }

    /**
     * @return the protein sensitivity for deep analysis
     */
    public double getSensitivity() {
        return sensitivity;
    }

    /**
     * Specify the sensitivity for protein comparisons.
     *
     * @param newValue	proposed new sensitivity
     */
    public void setSensitivity(double newValue) {
        this.sensitivity = newValue;
    }

    /**
     * Attempt to improve a genome by deleting suspicious contigs.
     *
     * @param genome	genome to improve
     * @param gReport	evaluation report for the genome
     * @param analysis	genome quality analysis
     * @param subFile	subsystem projector file
     *
     * @return TRUE if the genome was modified, else FALSE
     *
     * @throws IOException
     */
    protected boolean improve(Genome genome, GenomeStats gReport, GenomeAnalysis analysis, File subFile) throws IOException {
        boolean retVal = false;
        // Only bother to do this if the genome is contaminated.
        if (gReport.isClean())
            log.info("Genome is not contaminated.  No improvement attempted.");
        else {
            log.info("Attempting to improve genome {}.", genome);
            boolean subsOK = true;
            int deletes = 0;
            int dnaDelete = 0;
            // Loop through the bad contigs.
            for (Contig contig : analysis.getBadContigs()) {
                // Delete this contig.
                log.info("Deleting contig {}.", contig.getId());
                boolean ok = genome.deleteContig(contig);
                if (! ok) subsOK = false;
                // Denote the genome has changed.
                deletes++;
                dnaDelete += contig.length();
            }
            log.info("{} contigs with total length {} deleted.", deletes, dnaDelete);
            if (! subsOK) {
                // Here a subsystem has gone invalid.  Re-project the subsystems.
                log.info("Subsystems must be regenerated after improvement.  Loading projector from {}.", subFile);
                var projector = SubsystemProjector.load(subFile);
                log.info("Projecting subsystems onto {}.", genome);
                projector.project(genome);
                log.info("{} subsystems in genome.", genome.getSubsystems().size());
            }
            retVal = (deletes > 0);
            if (retVal) {
                // Since we updated the genome, we need to add an analysis event.
                AnalysisEvent event = new AnalysisEvent("dl4j.eval improve", this);
                genome.addEvent(event);
            }
        }
        return retVal;
    }

    /**
     * @return the refEngine
     */
    public RefGenomeComputer getRefEngine() {
        return refEngine;
    }

    /**
     * @param refEngine the refEngine to set
     */
    public void setRefEngine(RefGenomeComputer refEngine) {
        this.refEngine = refEngine;
    }

}
