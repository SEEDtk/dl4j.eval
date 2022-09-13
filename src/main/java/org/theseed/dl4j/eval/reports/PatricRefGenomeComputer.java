/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.kmers.KmerCollectionGroup;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * This class uses representative genomes to compute the reference genome for an incoming
 * evaluation.  The reference genomes are downloaded from PATRIC.
 *
 * @author Bruce Parrello
 *
 */
public class PatricRefGenomeComputer extends RefGenomeComputer {

    // FIELDS

    /** kmer collection table for computing reference genomes */
    private KmerCollectionGroup referenceGenomes;
    /** buffer of reference genomes in memory */
    private Map<String, Genome> referenceBuffer;
    /** connection to PATRIC */
    private P3Connection p3;

    /**
     * Initialize the reference genome computer from the files in a specified
     * evaluation directory.
     *
     * @param modelDir	source evaluation directory
     * @throws IOException
     */
    public PatricRefGenomeComputer(File modelDir) throws IOException {
        super();
        // Connect to PATRIC.
        this.p3 = new P3Connection();
        // Create the reference buffer.
        this.referenceBuffer = new HashMap<String, Genome>();
        // Read in the reference-genome database.
        File refGenomeFile = new File(modelDir, "refGenomes.fa");
        try (FastaInputStream refStream = new FastaInputStream(refGenomeFile)) {
            // Create the kmer database object.
            this.referenceGenomes = new KmerCollectionGroup();
            log.info("Reading reference genomes from {}", refGenomeFile);
            for (Sequence seq : refStream) {
                this.referenceGenomes.addSequence(seq, seq.getLabel());
            }
            log.info("{} reference genomes read.", this.referenceGenomes.size());
        }

    }

    @Override
    public void initialize(GenomeStats[] reports) {
        // Clear the reference buffer.
        this.referenceBuffer.clear();
        // Loop through the incoming reports, mapping them to reference genomes.
        for (GenomeStats report : reports) {
            if (report != null) {
                Genome refGenome = this.computeRef(report);
                if (refGenome != null) {
                    this.put(report.getId(), refGenome);
                }
            }
        }
    }

    /**
     * Compute the reference genome for an incoming genome.
     *
     * @param gReport	evaluation report for the genome whose reference is desired
     *
     * @return a genome object for the reference genome, or NULL if there is none
     */
    public Genome computeRef(GenomeStats gReport) {
        Genome retVal = null;
        String seedProt = gReport.getSeed();
        // Find the closest genome in the reference genome database.
        log.info("Computing reference genome for {}: {}", gReport.getId(), gReport.getName());
        KmerCollectionGroup.Result refGenomeData = this.referenceGenomes.getBest(seedProt);
        double refGenomeDistance = refGenomeData.getDistance();
        String refGenomeId = refGenomeData.getGroup();
        if (refGenomeId == null) {
            log.info("No reference genome found for {}.", gReport.getId());
        } else if (refGenomeDistance > MAX_GENOME_DIST) {
            log.info("Reference genome {} found, but distance of {} exceeds the maximum of {}.", refGenomeId, refGenomeDistance, MAX_GENOME_DIST);
        } else if (refGenomeId.contentEquals(gReport.getId())) {
            log.info("Genome {} is a reference to itself.  No useful data is available.", refGenomeId);
        } else {
            // Read in the genome and buffer it in case we reuse it.
            retVal = this.referenceBuffer.computeIfAbsent(refGenomeId, k -> P3Genome.load(p3, k, P3Genome.Details.PROTEINS));
        }
        return retVal;
    }
}
