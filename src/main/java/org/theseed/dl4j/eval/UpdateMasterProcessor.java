/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;
import java.util.HashSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.GenomeEval;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeMultiDirectory;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Genome;
import org.theseed.utils.ParseFailureException;

/**
 * This command refreshes the PATRIC master evaluation table after the PATRIC master GTO directory has been altered.
 * This should only be used if the evaluator has not changed. The old evaluation file is read from the standard input
 * and echoed to the standard output.  Any genomes no longer in the master directory are skipped.  Then, the remaining
 * genomes are evaluated and the results added to the end.
 *
 * The positional parameters are the name of the evaluation directory and the name of the PATRIC master GTO directory.
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	old master evaluation table (if not STDIN)
 * -o	new master evaluation table (if not STDOUT)
 * -b	batch size (default 200)
 *
 * @author Bruce Parrello
 *
 */
public class UpdateMasterProcessor extends BaseEvaluator {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(UpdateMasterProcessor.class);
    /** input stream for old evaluation master table */
    private TabbedLineReader reader;
    /** output stream for new evaluation master table */
    private PrintWriter writer;
    /** genome master directory source */
    private GenomeMultiDirectory genomes;
    /** TRUE if we have a completeness evaluator, else FALSE */
    private boolean haveCompleteness;

    // COMMAND-LINE OPTIONS

    /** input file (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "inFile.tbl", usage = "input file (if not STDIN)")
    private File inFile;

    /** output file (if not STDOUT) */
    @Option(name = "--output", aliases = { "-o" }, usage = "output file for report (if not STDOUT)")
    private File outFile;

    /** number of genomes per batch */
    @Option(name = "--batch", aliases = { "-b", "--batchSize" }, metaVar = "100", usage = "number of genomes to process at once")
    private int batchSize;

    @Argument(index = 1, metaVar = "P3Master", usage = "name of the PATRIC master genome directory")
    private File genomeDir;
    private int iGenome;


    @Override
    public P3Genome.Details getDetailLevel() {
        return P3Genome.Details.STRUCTURE_ONLY;
    }

    @Override
    protected void setDefaults() {
        this.batchSize = 200;
        this.inFile = null;
        this.reader = null;
        this.outFile = null;
        this.writer = null;
    }

    @Override
    protected void setHaveCompleteness(boolean exists) {
        this.haveCompleteness = exists;
    }

    @Override
    protected void validateOutputParms() throws ParseFailureException, IOException {
        // Set up the genome source.
        if (! this.genomeDir.isDirectory())
            throw new FileNotFoundException("Genome master directory " + this.genomeDir + " is not found or invalid.");
        else {
            log.info("Genomes are stored in master directory {}.", this.genomeDir);
            this.genomes = new GenomeMultiDirectory(this.genomeDir);
            log.info("{} genomes found.", this.genomes.size());
        }
    }

    @Override
    public void validateEvalParms() throws IOException, ParseFailureException {
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be greater than 0.");
        // Set up the input file.
        if (this.inFile == null) {
            log.info("Old evaluation data will be taken from the standard input.");
            this.reader = new TabbedLineReader(System.in);
        } else if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
        else {
            log.info("Old evaluation data will be read from {}.", this.inFile);
            this.reader = new TabbedLineReader(this.inFile);
        }
        // Set up the output file.  Note we turn on auto-flush.
        if (this.outFile == null) {
            log.info("New evaluation data will be written to the standard output.");
            this.writer = new PrintWriter(System.out, true);
        } else {
            log.info("New evaluation data will be written to {}.", this.outFile);
            var fileWriter = new FileWriter(this.outFile);
            this.writer = new PrintWriter(fileWriter, true);
        }
    }

    @Override
    protected void runCommand() throws Exception {
        try {
            // Read in the role maps.  This also sets the all-important haveCompleteness flag.
            this.initializeData();
            // Now we copy and filter the input file.  First, we need to output a header line.
            this.writer.println(GenomeEval.getHeader(this.haveCompleteness));
            // We will accumulate the genome IDs from the input file here.
            var genomesProcessed = new HashSet<String>(2000);
            // Now, loop through the input.
            int linesIn = 0;
            int linesDeleted = 0;
            log.info("Processing genomes from the previous run.");
            for (var line : this.reader) {
                linesIn++;
                // Get the genome ID.
                String genomeId = line.get(GenomeEval.GENOME_COL);
                if (! this.genomes.contains(genomeId)) {
                    linesDeleted++;
                    log.debug("Genome {} ({}) is no longer present.", genomeId, line.get(GenomeEval.NAME_COL));
                } else {
                    // Echo the line and save the genome ID.
                    this.writer.println(line.getAll());
                    genomesProcessed.add(genomeId);
                }
            }
            log.info("{} of {} genomes retained from input.  {} deleted.", genomesProcessed.size(), linesIn, linesDeleted);
            // Allocate our arrays.
            this.allocateArrays(this.batchSize);
            this.iGenome = 0;
            // Loop through the master directory.
            long start = System.currentTimeMillis();
            int linesAdded = 0;
            int genomesIn = 0;
            int batchCount = 0;
            for (String genomeId : this.genomes.getIDs()) {
                genomesIn++;
                if (! genomesProcessed.contains(genomeId)) {
                    // Here we need to evaluate a new genome.  Insure there is room in the current batch.
                    if (this.iGenome >= this.batchSize) {
                        this.processBatch();
                        batchCount++;
                        if (log.isInfoEnabled()) {
                            long millisPerGenome = ((System.currentTimeMillis() - start)) / linesAdded;
                            Duration d = Duration.ofMillis(millisPerGenome);
                            log.info("{} genomes scanned of {}. {} per genome.", genomesIn, this.genomes.size(), d.toString());
                        }
                    } else {
                        linesAdded++;
                        Genome genome = this.genomes.get(genomeId);
                        this.processGenome(this.iGenome, genome);
                        this.iGenome++;
                    }
                }
            }
            if (iGenome > 0) {
                this.processBatch();
                batchCount++;
            }
            log.info("{} new genomes evaluated of {} in directory.  {} batches required.", linesAdded, this.genomes.size(), batchCount);
        } finally {
            // Insure the output and input files are closed.
            this.writer.close();
            this.reader.close();
        }

    }

    /**
     * Evaluate the current batch of genomes.
     *
     * @throws IOException
     */
    private void processBatch() throws IOException {
        log.info("Processing genome batch with {} genomes.", this.iGenome);
        this.setnGenomes(this.iGenome);
        // Evaluate the consistency of the genomes.
        this.evaluateConsistency();
        // Write the results.
        for (int i = 0; i < this.getGenomeCount(); i++) {
            GenomeStats gReport = this.getGReport(i);
            String outputLine = gReport.formatStandardOutputLine(this.haveCompleteness);
            this.writer.println(outputLine);
        }
        // Empty the batch.
        this.iGenome = 0;
    }

}
