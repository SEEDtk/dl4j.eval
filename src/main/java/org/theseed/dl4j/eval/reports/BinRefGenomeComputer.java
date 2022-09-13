/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;

import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;

/**
 * This object computes the reference genomes for a binning job.  For each binning directory, the
 * "setSampleDir" method should be called before producing reports.  The reference genome will be
 * loaded from the appropriate json file in the binning directory.
 *
 * @author Bruce Parrello
 *
 */
public class BinRefGenomeComputer extends RefGenomeComputer {

    // FIELDS
    /** binning sample directory */
    private File binDir;

    /**
     * Specify the binning directory containing the sample's files.
     *
     * @param binDir	current binning directory
     */
    public void setSampleDir(File binDir) {
        this.binDir = binDir;
    }
    @Override
    protected void initialize(GenomeStats[] reports) {
        // Loop through the binned genomes.
        for (GenomeStats gReport : reports) {
            if (gReport != null) {
                Genome genome = gReport.getGenome();
                // Get the reference genome ID.
                String refGenomeId = genome.getBinRefGenomeId();
                if (refGenomeId != null) {
                    // Compute the file containing the reference genome.
                    File refGenomeFile = new File(this.binDir, refGenomeId + ".json");
                    if (refGenomeFile.canRead()) {
                        try {
                            // Attach the reference genome to this genome.
                            Genome refGenome = new Genome(refGenomeFile);
                            this.put(genome.getId(), refGenome);
                        } catch (IOException e) {
                            throw new UncheckedIOException("Error loading reference genome " + refGenomeId + ".", e);
                        }
                    }
                }
            }
        }
    }


}
