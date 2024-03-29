/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.time.Duration;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.stats.GenomeEval;

/**
 * This command performs a mass evaluation on a genome source (PATRIC ID file, master genome directory,
 * GTO directory) and produces a sequential flat-file output.  It is designed for situations where a
 * very large number of genomes needs to be processed.  It provides minimal information, but is restartable
 * and efficient.
 *
 * The positional parameters are the name of the evaluation directory and the name of the input directory.
 * The output is normally to the standard output, but is to the resume file if "--resume" is specified.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed log messages
 * -b	batch size (default 200)
 *
 * --resume		if specified, the name of a previous run's output file; it will be assumed the run was interrupted
 * 				and it will be continued with the next genome not yet processed
 * --source		type of input-- master genome directory, GTO directory, or file of PATRIC genome IDs
 *
 * @author Bruce Parrello
 *
 */
public class MassEvalProcessor extends BaseEvaluator {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MassEvalProcessor.class);
    /** input genome source */
    private GenomeSource master;
    /** output stream */
    private PrintStream outStream;
    /** TRUE if we have completeness data */
    private boolean haveCompleteness;
    /** set of genome IDs to skip */
    private Set<String> skipSet;
    /** number of completeness-related columns */
    private static final int COMPLETENESS_COLUMNS = 3;


    // COMMAND-LINE OPTIONS

    /** number of genomes per batch */
    @Option(name = "--batch", aliases = { "-b", "--batchSize" }, metaVar = "100", usage = "number of genomes to process at once")
    private int batchSize;

    /** old output file if we are resuming */
    @Option(name = "--resume", metaVar = "output.log", usage = "if we are resuming, the output file from the interrupted run")
    private File resumeFile;

    /** type of genome source */
    @Option(name = "--source", aliases = { "-t", "--type" }, usage = "type of genome input (master genome directory, GTO directory, patric ID file)")
    private GenomeSource.Type inType;

    /** input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "input genome file or directory", required = true)
    private File inDir;

    @Override
    protected void setDefaults() {
        this.inType = GenomeSource.Type.MASTER;
        this.haveCompleteness = false;
        this.batchSize = 200;
    }

    /**
     * Check for a resume situation.
     *
     * @return the number of items already processed
     *
     * @throws IOException
     */
    protected void setup() throws IOException {
        // Check for a resume situation.
        if (this.resumeFile == null) {
            // Normal processing.  Put the report to the standard output.
            // We also need a header.
            int cols = GenomeEval.DEFAULT_HEADERS.length;
            if (! this.haveCompleteness) cols -= COMPLETENESS_COLUMNS;
            String header = IntStream.range(0,  cols).mapToObj(i -> GenomeEval.DEFAULT_HEADERS[i]).collect(Collectors.joining("\t"));
            System.out.println(header);
            this.outStream = System.out;
            this.skipSet = Collections.emptySet();
        } else {
            // Resume processing.  Save the genomes we've already seen.
            try (TabbedLineReader reader = new TabbedLineReader(this.resumeFile)) {
                int idColIdx = reader.findField("Genome");
                this.skipSet = new HashSet<String>(this.batchSize);
                for (TabbedLineReader.Line line : reader) {
                    this.skipSet.add(line.get(idColIdx));
                }
            }
            // Open the resume file for append-style output with autoflush.
            FileOutputStream outStream = new FileOutputStream(this.resumeFile, true);
            this.outStream = new PrintStream(outStream, true);
        }
    }

    @Override
    protected void runCommand() throws Exception {
        // Read in the role maps.
        this.initializeData();
        // Set up the input.
        this.setup();
        // Allocate our arrays.
        this.allocateArrays(this.batchSize);
        // Get the list of genome IDs to process.
        Set<String> genomeIDs = this.master.getIDs();
        genomeIDs.removeAll(this.skipSet);
        // Loop through the genomes.  Note we track the genome's index in genomeStats.
        long start = System.currentTimeMillis();
        int done = 0;
        int total = genomeIDs.size();
        int iGenome = 0;
        for (String genomeID : genomeIDs) {
            Genome genome = this.master.getGenome(genomeID);
            // Insure there is room for this genome.
            if (iGenome >= this.batchSize) {
                processBatch();
                done += iGenome;
                if (log.isInfoEnabled()) {
                    double millisPerGenome = ((double) (System.currentTimeMillis() - start)) / done;
                    Duration d = Duration.ofMillis((long) ((total - done) * millisPerGenome));
                    log.info("{} of {} genomes processed. {} remaining.", done, total, d.toString());
                }
                iGenome = 0;
            }
            // Store the genome.
            processGenome(iGenome, genome);
            // Prepare for the next one.
            iGenome++;
        }
        // Process the residual batch.
        this.setnGenomes(iGenome);
        processBatch();
    }

    /**
     * Process the current batch of genomes.
     *
     * @throws IOException
     */
    public void processBatch() throws IOException {
        log.info("Processing genome batch with {} genomes.", this.getGenomeCount());
        // Evaluate the consistency of the genomes.
        evaluateConsistency();
        // Write the results.
        for (int i = 0; i < this.getGenomeCount(); i++) {
            GenomeStats gReport = this.getGReport(i);
            String outputLine = gReport.formatStandardOutputLine(this.haveCompleteness);
            this.outStream.println(outputLine);
        }
        log.info("{} genomes processed, {} per genome.", this.getGenomesProcessed(), this.getSpeed());
    }

    @Override
    public Details getDetailLevel() {
        return P3Genome.Details.STRUCTURE_ONLY;
    }

    @Override
    protected void validateOutputParms() throws ParseFailureException, IOException {
    }

    @Override
    protected void setHaveCompleteness(boolean exists) {
        this.haveCompleteness = exists;
    }

    @Override
    public void validateEvalParms() throws IOException, ParseFailureException {
        // Set up the input stream.
        log.info("Genomes will be read from {} source {}.", this.inType, this.inDir);
        this.master = this.inType.create(this.inDir);
        this.master.setDetailLevel(this.getDetailLevel());
    }

}
