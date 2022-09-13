/**
 *
 */
package org.theseed.dl4j.eval.reports;

import org.theseed.dl4j.eval.stats.GenomeStats;

/**
 * @author Bruce Parrello
 *
 */
public class NullRefGenomeComputer extends RefGenomeComputer {

    @Override
    protected void initialize(GenomeStats[] reports) {
        for (GenomeStats report : reports)
            this.put(report.getId(), null);
    }

}
