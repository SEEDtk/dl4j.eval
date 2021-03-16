/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.util.Arrays;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.GenomeMultiDirectory;
import org.theseed.p3api.Connection;
import org.theseed.p3api.Connection.Table;
import org.theseed.p3api.Criterion;
import org.theseed.p3api.P3Genome;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

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
    private Connection p3;
    /** list of domains to use */
    private static final List<String> DOMAINS = Arrays.asList("Bacteria", "Archaea");

    // COMMAND-LINE OPTIONS

    /** TRUE to erase any old genomes */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before starting")
    private boolean clearFlag;

    /** output directory name */
    @Argument(index = 0, metaVar = "outDir", usage = "output directory name")
    private File outDir;

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
        this.p3 = new Connection();
        // Get a list of the public, prokaryotic genomes.
        List<JsonObject> genomes = this.p3.getRecords(Table.GENOME, "kingdom", DOMAINS, "genome_id,genome_name", Criterion.EQ("public", "1"));
        log.info("{} genomes found in PATRIC.  {} already in output directory.", genomes.size(), this.gOutDir.size());
        // Loop through the genomes.
        int processed = 0;
        int downloaded = 0;
        long start = System.currentTimeMillis();
        for (JsonObject gRequest : genomes) {
            String genomeId = Connection.getString(gRequest, "genome_id");
            String name = Connection.getString(gRequest, "genome_name");
            processed++;
            if (this.gOutDir.contains(genomeId))
                log.info("Skipping {}. {} ({}):  already in output directory.", processed, genomeId, name);
            else {
                log.info("Processing {}. {} ({}).", processed, genomeId, name);
                P3Genome genome = P3Genome.Load(p3, genomeId, P3Genome.Details.STRUCTURE_ONLY);
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
