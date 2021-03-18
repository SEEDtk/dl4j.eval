/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.utils.ParseFailureException;


/**
 * This processor performs the evaluations on GTO files.  It will take as input an evaluation directory and an input directory of
 * GTO files.  The GTO files do not need to be fully-loaded.  Only the structural information is needed.  In addition,
 * the protein translation must be present for the seed protein and the SSU rRNA must be filled in.  Representative-genome
 * technology will be used for completeness checking, and only the length and ID is necessary for each contig.
 *
 * The positional parameters are the name of the evaluation directory and the name of the input directory.
 *
 * The input can be a file containing PATRIC genome IDs, a genome master directory, or a directory containing plain-text
 * GTOs.  The --source parameter determines which type of input is expected.
 *
 * When we are done, the output directory will contain a report file for each genome (named using the genome
 * ID) and a summary report (named "summary").  The file type depends on the type of reporting requested.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 * -O	the name of the output directory (defaults to the current directory)
 * -s	maximum distance for a close protein in deep reporting (0 to 1)
 * -b   maximum number of genomes to process at one time
 *
 * --terse		do not write the individual output files, only the summary
 * --update		store the quality information in the GTO and write it to the output directory
 * --format		specify the output format-- HTML, DEEP, or TEXT
 * --ref		name of a reference genome file to use for the deep-format reports
 * --clear		clear the output directory before processing
 * --source		type of input (PATRIC, DIR, MASTER)
 *
 * @author Bruce Parrello
 *
 */
public class EvalProcessor extends Evaluator  {

    // FIELDS
    /** genome source iterator */
    private GenomeSource genomeDir;

    // COMMAND LINE

    /** update-GTO flag */
    @Option(name = "-u", aliases = { "--update" }, usage = "store results in the input GTO")
    private boolean update;

    /** input file/directory name */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory", required = true)
    private File inDir;

    /** output directory */
    @Option(name = "-O", aliases = { "--outDir" }, metaVar = "outDir", usage = "output directory")
    private File outputDir;

    /** file of reference genome GTO mappings */
    @Option(name = "--ref", usage = "file of taxon ID to reference-genome GTO mappings")
    private File refGenomeFile;

    /** clear-output flag */
    @Option(name = "--clear", usage = "clear output directory before starting")
    private boolean clearOutputDir;

    /** number of genomes to process in each chunk */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "100", usage = "number of genomes to process in each batch")
    private int batchSize;

    /** type of input */
    @Option(name = "--source", usage = "type of genome input")
    private GenomeSource.Type inType;

    @Override
    public void validateEvalParms() throws IOException, ParseFailureException {
        // Check the input directory.
        if (! this.inDir.isDirectory()) {
            throw new FileNotFoundException("Input " + this.inDir + " is neither a directory or a readable file.");
        }
        // Set up the input stream.
        log.info("Genomes will be read from {} source {}.", this.inType, this.inDir);
        this.genomeDir = this.inType.create(this.inDir);
        // Set up the output directory.
        this.validateOutputDir(this.outputDir, this.clearOutputDir);
    }

    @Override
    public void setDefaults() {
        this.update = false;
        this.refGenomeFile = null;
        this.outputDir = new File(System.getProperty("user.dir"));
        this.batchSize = 200;
        this.inType = GenomeSource.Type.DIR;
    }

    @Override
    public void runCommand() throws Exception {
        // Set up the reference-genome engine (if necessary).
        this.setupRefGenomeEngine(this.refGenomeFile);
        // Read in the role maps.
        initializeData();
        // Allocate our arrays.
        this.allocateArrays(this.batchSize);
        // Loop through the genomes.  Note we track the genome's index in genomeStats;
        int iGenome = 0;
        for (Genome genome : genomeDir) {
            // Insure there is room for this genome.
            if (iGenome >= this.batchSize) {
                processBatch();
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
        // Finish processing.
        this.close();
    }

    /**
     * Process the current batch of genomes.
     *
     * @throws IOException
     * @throws FileNotFoundException
     * @throws NoSuchAlgorithmException
     * @throws UnsupportedEncodingException
     */
    public void processBatch()
            throws IOException, FileNotFoundException, NoSuchAlgorithmException, UnsupportedEncodingException {
        log.info("Processing genome batch with {} genomes.", this.getGenomeCount());
        // Evaluate the consistency of the genomes.
        evaluateConsistency();
        // Write the results.
        writeOutput();
        // If we are updating GTOs, do it here.
        if (this.update) {
            // We need the version string.
            String version = this.getVersion();
            for (int g = 0; g < this.getGenomeCount(); g++) {
                GenomeStats gReport = this.getReport(g);
                Genome gObject = gReport.getGenome();
                String gId = gObject.getId();
                log.debug("Updating GTO for {}.", gId);
                gReport.store(gObject, this.getRoleDefinitions(), version, this.getOptions());
                File outFile = new File(this.getOutDir(), gId + ".gto");
                gObject.save(outFile);
            }
        }
    }

}
