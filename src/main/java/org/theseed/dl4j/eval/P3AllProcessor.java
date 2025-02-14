/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.GenomeMultiDirectory;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;

import com.github.cliftonlabs.json_simple.JsonObject;


/**
 * This command creates a directory of evaluation GTOs for all the public, prokaryotic genomes in PATRIC.  This creates an
 * enormous number of files, so they are stored in compressed form and split between multiple directories.  A master index is
 * kept at the top directory level.
 *
 * This is a resumable operation.  The default is to add genomes that are not already present (unless --clear is specified).
 *
 * The positional parameter is the output directory name.  The command-line options are as follows.
 *
 * -h	show command-line usage
 * -v	show more detailed log messages
 *
 * --clear	erase the output directory before starting
 *
 *
 * @author Bruce Parrello
 *
 */
public class P3AllProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(P3AllProcessor.class);
    /** output directory controller */
    private GenomeMultiDirectory gOutDir;
    /** connection to PATRIC */
    private P3Connection p3;

    // COMMAND-LINE OPTIONS

    /** TRUE to erase any old genomes */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before starting")
    private boolean clearFlag;

    /** output directory name */
    @Argument(index = 0, metaVar = "outDir", usage = "output directory name")
    private File outDir;

    /**
     * This class sorts genome records by ID, putting records already in the database first.
     */
    private class GenomeSorter implements Comparator<JsonObject> {

        @Override
        public int compare(JsonObject o1, JsonObject o2) {
            int retVal = 0;
            String k1 = KeyBuffer.getString(o1, "genome_id");
            String k2 = KeyBuffer.getString(o2, "genome_id");
            boolean b1 = P3AllProcessor.this.gOutDir.contains(k1);
            boolean b2 = P3AllProcessor.this.gOutDir.contains(k2);
            if (b1 == b2)
                retVal = k1.compareTo(k2);
            else if (b1)
                retVal = -1;
            else
                retVal = 1;
            return retVal;
        }

    }

    @Override
    protected void setDefaults() {
        this.clearFlag = false;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Prepare the output directory.
        if (! this.outDir.isDirectory()) {
            this.gOutDir = GenomeMultiDirectory.create(this.outDir, false);
        } else if (this.clearFlag) {
            this.gOutDir = GenomeMultiDirectory.create(this.outDir, true);
        } else {
            this.gOutDir = new GenomeMultiDirectory(this.outDir);
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Connect to PATRIC.
        this.p3 = new P3Connection();
        // Get a list of the public, prokaryotic genomes.
        SortedSet<JsonObject> genomes = new TreeSet<JsonObject>(this.new GenomeSorter());
        this.p3.addAllProkaryotes(genomes);
        log.info("{} genomes found in PATRIC.  {} already in output directory.", genomes.size(), this.gOutDir.size());
        // Loop through the genomes, removing the ones already processed.
        int processed = 0;
        int downloaded = 0;
        long start = System.currentTimeMillis();
        for (JsonObject gRequest : genomes) {
            String genomeId = KeyBuffer.getString(gRequest, "genome_id");
            String name = KeyBuffer.getString(gRequest, "genome_name");
            processed++;
            if (this.gOutDir.contains(genomeId))
                log.info("Skipping {}. {} ({}):  already in output directory.", processed, genomeId, name);
            else {
                log.info("Processing {}. {} ({}).", processed, genomeId, name);
                P3Genome genome = P3Genome.load(p3, genomeId, P3Genome.Details.STRUCTURE_ONLY);
                if (genome == null)
                    log.warn("WARNING: genome {} not found in PATRIC.", genomeId);
                else {
                    this.gOutDir.add(genome);
                    downloaded++;
                }
                if (downloaded > 0) {
                    Duration duration = Duration.ofMillis(System.currentTimeMillis() - start);
                    log.info("{} of {} genomes processed, {} per genome downloaded.", processed, genomes.size(),
                             duration.dividedBy(downloaded));
                }
            }
        }
    }

}
