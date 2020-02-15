/**
 *
 */
package org.theseed.genome;

import java.io.File;
import java.io.FilenameFilter;
import java.util.regex.Pattern;

/**
 * This is a file-name filter that only accepts bin GTO files.
 *
 * @author Bruce Parrello
 *
 */
public class BinFilter implements FilenameFilter {

    private static final Pattern BIN_GTO_PATTERN = Pattern.compile("bin\\d+\\.gto");

    @Override
    public boolean accept(File dir, String name) {
        return (BIN_GTO_PATTERN.matcher(name).matches());
    }

}
