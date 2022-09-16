/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.utils.ParseFailureException;

/**
 * This command evaluates binning results.  The binning directory must contain a bins.json describing the bins, and a series
 * of binX.XXXXXX.genome files containing the annotated bins.  The bins.json file will be used to associate reference genomes
 * with the bin genomes, which will then be loaded from the XXXXXX.X.json files in the binning directory.  The bins will then
 * be deep-evaluated, and optionally improved.  The output directory will contain the evaluated/improved bins as GTOs plus
 * the HTML evaluation reports.  Since we know all the reference genomes are in BV-BRC, many simplifying assumptions can be
 * made.
 *
 * The positional parameters are the name of the input binning directory and the name of the target output directory.
 *
 * The command-line options are
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -s	maximum distance for a protein to be considered close
 *
 *
 * @author Bruce Parrello
 *
 */
public class BinEvalProcessor extends BaseEvaluator {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinEvalProcessor.class);

    @Override
    public Details getDetailLevel() {
        // TODO code for getDetailLevel
        return null;
    }

    @Override
    protected void validateOutputParms() throws ParseFailureException, IOException {
        // TODO code for validateOutputParms

    }

    @Override
    protected void setHaveCompleteness(boolean exists) {
        // TODO code for setHaveCompleteness

    }

    @Override
    public void validateEvalParms() throws IOException, ParseFailureException {
        // TODO code for validateEvalParms

    }

    @Override
    protected void setDefaults() {
        // TODO code for setDefaults

    }

    @Override
    protected void runCommand() throws Exception {
        // TODO code for runCommand

    }
    // FIELDS
    // TODO data members for BinEvalProcessor

    // TODO constructors and methods for BinEvalProcessor
}
