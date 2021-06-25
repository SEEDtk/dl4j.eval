/**
 *
 */
package org.theseed.dl4j.eval.reports;

import org.theseed.dl4j.eval.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.p3api.Connection;
import org.theseed.p3api.P3Genome;

/**
 * This is a simplified reference genoem ID processor that always returns the same genome ID.
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
        Connection p3 = new Connection();
        this.refGenome = P3Genome.load(p3, refGenomeId, P3Genome.Details.PROTEINS);
    }

    @Override
    protected void initialize(GenomeStats[] reports) {
        for (GenomeStats report : reports)
            this.put(report.getId(), this.refGenome);
    }

}
