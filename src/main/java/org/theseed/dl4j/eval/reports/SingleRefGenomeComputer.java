/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.IOException;
import java.io.UncheckedIOException;

import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.p3api.P3Genome;

/**
 * This is a simplified reference genome ID processor that always returns the same genome ID.
 *
 * @author Bruce Parrello
 *
 */
public class SingleRefGenomeComputer extends RefGenomeComputer {

    // FIELDS
    /** genome to always use */
    private Genome refGenome;

    /**
     * Construct a reference-genome computer that always returns the same genome ID.
     *
     * @param refGenomeId	ID to always return
     */
    public SingleRefGenomeComputer(String refGenomeId) {
        // Connect to PATRIC.
        P3CursorConnection p3 = new P3CursorConnection();
        try {
            this.refGenome = P3Genome.load(p3, refGenomeId, P3Genome.Details.PROTEINS);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    protected void initialize(GenomeStats[] reports) {
        for (GenomeStats report : reports) {
            if (report != null)
                this.put(report.getId(), this.refGenome);
        }
    }

}
