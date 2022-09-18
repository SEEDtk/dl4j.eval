/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.eval.reports.BinRefGenomeComputer;
import org.theseed.dl4j.eval.reports.EvalDeepReporter;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.utils.ParseFailureException;
import org.threeten.bp.Duration;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This command evaluates binning results.  The binning directory must contain a bins.json describing the bins, and a series
 * of Bin.X.XXXXXX.genome files containing the annotated bins.  The bins.json file will be used to associate reference genomes
 * with the bin genomes, which will then be loaded from the XXXXXX.X.json files in the binning directory.  The bins will then
 * be deep-evaluated, and optionally improved.  The output directory will contain the evaluated/improved bins as GTOs plus
 * the HTML evaluation reports.  Since we know all the reference genomes are in BV-BRC, many simplifying assumptions can be
 * made.
 *
 * The positional parameters are the name of the model directory, the name of the input binning directory and the
 * name of the target output directory.
 *
 * The command-line options are
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -s	maximum distance for a protein to be considered close (default 0.8)
 *
 * --clear		erase the output directory before processing
 *
 * @author Bruce Parrello
 *
 */
public class BinEvalProcessor extends BaseEvaluator {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinEvalProcessor.class);
    /** reporting facility */
    private EvalDeepReporter reporter;
    /** list of input genome files */
    private File[] binFiles;
    /** subsystem projector file */
    private File projectorFile;

    /** sensitivity of protein comparison on reports */
    @Option(name = "-s", aliases = { "--sensitivity" }, metaVar = "0.9", usage = "maximum distance for a close protein")
    private double sensitivityLevel;

    /** erase the output directory before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** input binning directory */
    @Argument(index = 1, metaVar = "binDir", usage = "input binning directory containing GTOs and binning data")
    private File binDir;

    /** output directory for results */
    @Argument(index = 2, metaVar = "outDir", usage = "output directory for updated GTOs and reports")
    private File outDir;

    /**
     * File name filter for the annotated bin (.genome) files.
     */
    private static class GenomeFileFilter implements FilenameFilter {

        /** file name pattern for bin files */
        private static final Pattern BIN_FILE_NAME = Pattern.compile("Bin\\.\\d+\\.\\d+\\.genome");

        @Override
        public boolean accept(File dir, String name) {
            return BIN_FILE_NAME.matcher(name).matches();
        }

    }

    @Override
    protected void setDefaults() {
        this.sensitivityLevel = 0.8;
        this.clearFlag = false;
    }

    @Override
    public Details getDetailLevel() {
        return Details.PROTEINS;
    }

    @Override
    protected void validateOutputParms() throws ParseFailureException, IOException {
        // Insure the output directory is valid.
        this.validateOutputDir(this.outDir, this.clearFlag);
        // Initialize the reporter.
        this.reporter = new EvalDeepReporter();
        this.reporter.setOutDir(this.outDir);
        log.info("Output will be to {}.", this.outDir);
    }

    @Override
    protected void setHaveCompleteness(boolean exists) {
        this.reporter.setHaveCompleteness(exists);
    }

    @Override
    public void validateEvalParms() throws IOException, ParseFailureException {
        // Insure the input directory is valid.
        if (! this.binDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.binDir + " is not found or invalid.");
        // Insure we have a bins.json file.
        File binsJsonFile = new File(this.binDir, "bins.json");
        if (! binsJsonFile.canRead())
            throw new FileNotFoundException("Input directory " + this.binDir + " does not appear to be a binning directory: no bins.json file.");
        // Verify the subsystem projector.
        this.projectorFile = new File(this.getModelDir(), "projector.ser");
        if (! this.projectorFile.canRead())
            throw new FileNotFoundException("Subsystem projector in " + this.getModelDir() + " not found or unreadable.");
        // Create the reference engine and save the sample directory.
        var refEngine = new BinRefGenomeComputer();
        refEngine.setSampleDir(this.binDir);
        this.setRefEngine(refEngine);
        // Get the bin files themselves, sorted by bin number.
        this.binFiles = this.binDir.listFiles(new GenomeFileFilter());
        log.info("{} bin files found in {}.", this.binFiles.length, this.binDir);
    }

    @Override
    protected void runCommand() throws Exception {
        // Here we build a list of the GTOs for the bins.  The bin number is stored in the GTO as a quality field.
        // This number corresponds to the bin's position in the bins.json file.
        log.info("Loading bin genomes from {}.", this.binDir);
        List<Genome> binGenomes = new ArrayList<Genome>(this.binFiles.length);
        for (File binFile : this.binFiles) {
            // Compute the bin number.  The second piece is the bin number as a string.  This is enforced by the
            // file name pattern in the filter.
            String[] binNumString = StringUtils.split(binFile.getName(), '.');
            int binNum = Integer.valueOf(binNumString[1]);
            // Read the genome.
            Genome binGenome = new Genome(binFile);
            log.info("Bin {} is {}.", binNum, binGenome);
            // Set the bin number and denote the genome is BV-BRC.
            JsonObject quality = binGenome.getQuality();
            quality.put("bin_num", binNum);
            binGenome.setHome("BV-BRC");
            // Save the bin as a genome to evaluate.
            binGenomes.add(binGenome);
        }
        // Initialize the data structurs.
        this.initializeData();
        // Set up the evaluation batch.
        final int nGenomes = binGenomes.size();
        if (nGenomes == 0)
            log.info("No bins to evaluate.");
        else {
            this.allocateArrays(nGenomes);
            long start = System.currentTimeMillis();
            for(int i = 0; i < nGenomes; i++)
                this.processGenome(i, binGenomes.get(i));
            this.setnGenomes(nGenomes);
            // Evaluate the consistency of the genomes.
            evaluateConsistency();
            // Analyze the genomes.
            var analyses = this.analyzeGenomes();
            long goodCount = Arrays.stream(analyses).filter(x -> x.isGood()).count();
            log.info("{} good bins found in initial analysis.", goodCount);
            // Check for possible improvements.
            int improveCount = 0;
            for (int i = 0; i < nGenomes; i++) {
                var gReport = this.getGReport(i);
                var analysis = analyses[i];
                if (analysis.hasRefGenome()) {
                    Genome genome = gReport.getGenome();
                    log.info("Attempting improvement of {}.", genome);
                    boolean improved = this.improve(gReport, analysis, this.projectorFile);
                    if (improved) {
                        this.processGenome(i, genome);
                        improveCount++;
                    }
                }
            }
            log.info("{} bins improved.", improveCount);
            // If we made improvements, rerun the evaluation.
            if (improveCount > 0) {
                this.evaluateConsistency();
                analyses = this.analyzeGenomes();
                goodCount = Arrays.stream(analyses).filter(x -> x.isGood()).count();
                log.info("{} good bins found after improvement.", goodCount);
            }
            // Write the genomes to the output directory and output the reports.
            log.info("Initializing reports.");
            this.reporter.open(this.getVersion(), this.getRoleDefinitions(), this.getModelDir());
           // Write the reports and the genomes.
            for (int i = 0; i < nGenomes; i++) {
                GenomeStats gReport = this.getGReport(i);
                Genome genome = gReport.getGenome();
                File outFile = new File(this.outDir, genome.getId() + ".gto");
                log.info("Writing genome {} to {}.", genome, outFile);
                genome.save(outFile);
                this.reporter.writeGenome(gReport, analyses[i]);
            }
            this.reporter.close();
            Duration time = Duration.ofMillis((System.currentTimeMillis() - start) / nGenomes);
            log.info("{} bins processed.  {} per genome.", nGenomes, time);
        }
    }

}
