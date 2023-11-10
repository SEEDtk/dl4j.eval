/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3SubsystemProjector;

/**
 * This is a simple command that saves a subsystem projector file downloaded from PATRIC.  This file can
 * then be used to project subsystems during genome improvement.
 *
 * The positional parameter is the name of the output file.
 *
 * @author Bruce Parrello
 *
 */
public class P3SubsystemSaveProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(P3SubsystemSaveProcessor.class);

    // COMMAND-LINE OPTIONS

    /** name of the output file */
    @Argument(index = 0, metaVar = "outFile.ser", usage = "name of the output file", required = true)
    private File outFile;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        log.info("Output file will be {}.", this.outFile);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        log.info("Connecting to PATRIC.");
        P3Connection p3 = new P3Connection();
        log.info("Downloading subsystem definitions.");
        var projector = new P3SubsystemProjector(p3);
        log.info("Writing projector to {}.", this.outFile);
        projector.save(this.outFile);
    }

}
