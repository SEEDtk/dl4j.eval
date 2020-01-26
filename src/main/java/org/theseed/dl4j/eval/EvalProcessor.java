/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

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
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.ICommand;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.LoggerContext;


/**
 * This processor performs the evaluations.  It will take as input an evaluation directory and an input directory of
 * GTO files.  The GTO files do not need to be fully-loaded.  Only the features and the genome ID are needed.  In addition,
 * the protein translation must be present for the seed protein.  Representative-genome technology will be used for
 * completeness checking, and only the length and ID is necessary for each contig.
 *
 * The positional parameters are the name of the evaluation directory, the name of the input directory, and the name of the
 * output directory.  The evaluation directory must contain the following files.
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
 * The input directory should contain a GTO file for each genome to evaluate.  Alternative, the second parameter can be the
 * name of an input GTO file, and then only that single file will be processed.
 *
 * When we are done, the output directory will contain a I<genomeID><code>.out</code> file for each genome.  This file will
 * contain the quality numbers plus the predicted and actual counts for each role.  In addition, a
 * summary of the results will be placed in "summary.tbl" in the same directory.
 *
 * The command-line options are as follows.
 *
 * -v	show more detailed progress messages
 *
 * --terse		do not write the individual output files, only the summary
 * --consume	delete the genome files as they are processed (only applies when the input is a directory)
 * --update		store the quality information in the GTO (mutually exclusive with consume)
 *
 * @author Bruce Parrello
 *
 */
public class EvalProcessor implements ICommand {

    /** logging facility */
    private static Logger log = LoggerFactory.getLogger(EvalProcessor.class);

    // FIELDS
    /** array of roles to use */
    private ArrayList<String> roles;
    /** array of genome statistics objects */
    private GenomeStats[] reports;
    /** in update mode, these are the saved Genomes so we can update */
    private Genome[] gtos;
    /** number of genomes read in this batch */
    private int nGenomes;
    /** output matrix, first index is genome, second is role */
    private int[][] rolesActual;
    /** array of flags indicating which roles have output */
    private boolean[] rolesUsed;
    /** role definition map */
    private RoleMap roleDefinitions;
    /** TRUE for single-gto mode */
    private boolean singleton;
    /** completeness database */
    private List<UniversalRoles> compList;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    @Option(name="-b", aliases = { "--batchSize", "--batch" }, metaVar = "1000", usage = "input batch size")
    private int batchSize;

    /** summary-only flag */
    @Option(name = "--terse", usage = "do not output detail files")
    private boolean terse;

    /** consume-input flag */
    @Option(name = "--consume", usage = "delete input files after they are processed")
    private boolean consume;

    /** debug-message flag */
    @Option(name = "-v", aliases = { "--verbose", "--debug" }, usage = "show more detailed progress messages")
    private boolean debug;

    /** update-GTO flag */
    @Option(name = "-u", aliases = { "--update" }, usage = "store results in the input GTO")
    private boolean update;

    /** model directory */
    @Argument(index = 0, metaVar = "modelDir", usage = "model directory", required = true)
    private File modelDir;

    /** input file name */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory", required = true)
    private File inDir;

    /** output directory */
    @Argument(index = 2, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

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
            } else {
                if (! this.inDir.isDirectory()) {
                    if (this.inDir.canRead()) {
                        this.singleton = true;
                        if (this.consume) {
                            log.error("WARNING: input file will not be deleted.");
                        }
                    } else {
                        throw new FileNotFoundException("Input " + this.inDir + " is neither a directory or a readable file.");
                    }
                }
                if (this.update && this.consume)
                    throw new RuntimeException("Cannot specify both --consume and --update.");
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
            if (this.debug) {
                // To get more progress messages, we set the log level in logback.
                LoggerContext loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory();
                ch.qos.logback.classic.Logger logger = loggerContext.getLogger("org.theseed");
                logger.setLevel(Level.toLevel("TRACE"));
            }
            long start = System.currentTimeMillis();
            // Read in the role maps.
            initializeData();
            // Are we doing one genome or many?
            if (! this.singleton) {
                log.info("Genomes will be read from {}.", this.inDir);
                GenomeDirectory genomeDir = new GenomeDirectory(this.inDir);
                // We know the number of genomes, so we can allocate our arrays.
                this.nGenomes = genomeDir.size();
                log.info("{} genomes found in directory.", this.nGenomes);
                this.rolesActual = new int[this.nGenomes][this.roles.size()];
                this.reports = new GenomeStats[this.nGenomes];
                this.gtos = new Genome[this.nGenomes];
                // Loop through the genomes.  Note we track the genome's index in genomeStats;
                int iGenome = 0;
                for (Genome genome : genomeDir) {
                    processGenome(iGenome, genome);
                    if (this.consume) {
                        File genomeFile = genome.getFile();
                        log.debug("Deleting {}.", genomeFile.getPath());
                        Files.delete(genomeFile.toPath());
                    }
                    // Prepare for the next genome.
                    iGenome++;
                }
            } else {
                // Here we are processing only one genome.
                this.nGenomes = 1;
                this.rolesActual = new int[1][this.roles.size()];
                this.reports = new GenomeStats[1];
                this.gtos = new Genome[1];
                // Read in the genome and process it.
                log.info("Reading genome from {}.", this.inDir.getPath());
                Genome genome = new Genome(this.inDir);
                processGenome(0, genome);
            }
            // Evaluate the consistency of the genomes.
            evaluateConsistency();
            // Write the results.
            writeOutput();
            String rate = String.format("%6.2f", (double) (System.currentTimeMillis() - start) / (this.nGenomes * 1000));
            log.info("{} genomes evaluated at {} seconds/genome.", this.nGenomes, rate);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Read in the support structures.
     *
     * @throws IOException
     */
    protected void initializeData() throws IOException {
        this.roles = new ArrayList<String>(2800);
        File rolesToUseFile = new File(this.modelDir, "roles.to.use");
        log.info("Reading consistency roles from {}.", rolesToUseFile.getPath());
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
        log.trace("Processing #{} {}: {}.", iGenome, genome.getId(), genome.getName());
        // Create the reporting object.
        GenomeStats gReport = new GenomeStats(genome.getId(), genome.getDomain());
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
        // Save the GTO information if we're updating.
        if (this.update) {
            this.gtos[iGenome] = genome;
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
        log.info("Analyzing roles.");
        for (int iRole = 0; iRole < this.roles.size(); iRole++) {
            String role = this.roles.get(iRole);
            File modelFile = new File(this.modelDir, role + ".ser");
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
     * @throws IOException
     */
    protected void writeOutput() throws IOException {
        // Write the summary file and the output files.
        try (PrintWriter outStream = new PrintWriter(new File(this.outDir, "summary.tbl"))) {
            log.info("Writing output for {} genomes.", this.nGenomes);
            outStream.println("Genome\tCoarse\tFine\tCompleteness\tContamination\tHypothetical\tContigs\tGood_Seed\tScore");
            for (int g = 0; g < this.nGenomes; g++) {
                String genome = this.reports[g].getId();
                GenomeStats gReport = this.reports[g];
                double coarsePct = gReport.getCoarsePercent();
                double finePct = gReport.getFinePercent();
                double completePct = gReport.getCompletePercent();
                double contamPct = gReport.getContaminationPercent();
                double hypoPct = gReport.getHypotheticalPercent();
                int contigs = gReport.getContigCount();
                double score = gReport.getScore();
                String goodSeed = (gReport.isGoodSeed() ? "Y" : "");
                File outFile = new File(this.outDir, genome + ".out");
                if (! this.terse) {
                    log.trace("Writing {}.", outFile.getPath());
                    try (PrintWriter genomeStream = new PrintWriter(outFile)) {
                        genomeStream.println("Role\tactual\tpredicted\tuniversal");
                        Collection<String> pprs = gReport.getProblematicRoles();
                        for (String role : pprs) {
                            GenomeStats.ProblematicRole ppr = gReport.getReport(role);
                            String uFlag = (ppr.isUniversal() ? "Y" : "");
                            genomeStream.format("%s\t%d\t%d\t%s%n", role, ppr.getActual(), ppr.getPredicted(), uFlag);
                        }
                        genomeStream.println();
                        genomeStream.format("*\tCompleteness Group\t%s\t%n", gReport.getGroup());
                        genomeStream.format("*\tCompleteness\t%2.2f\t%n", completePct);
                        genomeStream.format("*\tContamination\t%2.2f\t%n", contamPct);
                        genomeStream.format("*\tCoarse Consistency\t%2.2f\t%n", coarsePct);
                        genomeStream.format("*\tFine Consistency\t%2.2f\t%n", finePct);
                        genomeStream.format("*\tHypothetical Rate\t%2.2f\t%n", hypoPct);
                        genomeStream.format("*\tPLFAM Coverage\t%2.2f\t%n", gReport.getPlfamPercent());
                        genomeStream.format("*\tCDS Coverage\t%2.2f\t%n", gReport.getCdsPercent());
                        genomeStream.format("*\tDNA Length\t%d\t%n", gReport.getDnaSize());
                        genomeStream.format("*\tContig Count\t%d\t%n", gReport.getContigCount());
                        genomeStream.format("*\tHypothetical Roles\t%d\t%n", gReport.getHypoCount());
                        genomeStream.format("*\tContig Length L50\t%d\t%n", gReport.getL50());
                        genomeStream.format("*\tContig Length N50\t%d\t%n", gReport.getN50());
                        genomeStream.format("*\tCDS Count\t%d\t%n", gReport.getPegCount());
                        genomeStream.format("*\tPLFAM Protein Count\t%d\t%n", gReport.getPlfamCount());
                        genomeStream.format("*\tGood PheS found\t%s\t%n", (gReport.isGoodSeed() ? "Yes" : "No"));
                    }
                }
                outStream.format("%s\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8d\t%-8s\t%8.2f%n",
                        genome, coarsePct, finePct, completePct, contamPct, hypoPct, contigs, goodSeed, score);
                if (this.update) {
                    log.trace("Updating GTO for {}.", genome);
                    Genome gObject = this.gtos[g];
                    gReport.store(gObject.getJson(), this.roleDefinitions);
                    gObject.update(gObject.getFile());
                }
            }
        }
    }

}
