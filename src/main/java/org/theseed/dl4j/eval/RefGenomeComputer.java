/**
 *
 */
package org.theseed.dl4j.eval;

import java.util.HashMap;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;

/**
 * This is a base class for computation of reference genomes.  At the beginning of
 * each batch, it is called to build a map of reported genomes to reference
 * genome objects.  Then for each individual genome it is called to get
 * the relevant object.
 *
 * @author Bruce Parrello
 *
 */
public abstract class RefGenomeComputer {

    // FIELDS

    /** maximum acceptable distance for a reference genome */
    public static final double MAX_GENOME_DIST = 0.8;

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(EvalProcessor.class);

    /** master map of incoming genome IDs to reference genome objects */
    private Map<String, Genome>	master;

    public RefGenomeComputer() {
        this.master = new HashMap<String, Genome>();
    }

    /**
     * Set up the reference genomes for an evaluation batch.
     *
     * @param reports	array of evaluated genomes
     */
    public void setupReferences(GenomeStats[] reports) {
        // Erase our old map.
        this.master.clear();
        // Create the new map.
        this.initialize(reports);
    }

    /**
     * Initialize the reference genome map for a given batch of evaluations.
     *
     * @param reports	array of evaluated genomes
     */
    protected abstract void initialize(GenomeStats[] reports);

    /**
     * @return the reference genome for a specified genome to be displayed,
     * 		   or NULL if there is no useful reference
     *
     * @param genomeId		ID of the genome being displayed
     */
    public Genome ref(String genomeId) {
        Genome retVal = master.get(genomeId);
        if (retVal != null && genomeId.contentEquals(retVal.getId())) {
            log.info("{} is a reference for itself:  no useful information.", genomeId);
            retVal = null;
        }
        return retVal;
    }

    /**
     * Assign a reference genome to an incoming genome ID
     *
     * @param genomeId		ID of genome being evaluated
     * @param refGenome		reference genome to assign
     */
    protected void put(String genomeId, Genome refGenome) {
        this.master.put(genomeId, refGenome);
    }

}
