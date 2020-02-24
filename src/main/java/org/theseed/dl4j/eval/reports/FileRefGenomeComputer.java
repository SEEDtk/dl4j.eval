/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.theseed.dl4j.eval.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.genome.TaxItem;
import org.theseed.io.TabbedLineReader;

/**
 * This object creates a mapping from taxonomic IDs to GTO files for reference genomes.
 * When analyzing a genome, if one of its taxonomic lineage IDs is in the file,
 * the corresponding GTO will be used as the reference genome.  Otherwise, no reference
 * genome is used.
 *
 * @author Bruce Parrello
 *
 */
public class FileRefGenomeComputer extends RefGenomeComputer {

    /** map of taxon IDs to reference genomes */
    private Map<String, Genome> refMap;

    /**
     * Create a file-based reference-genome engine.
     *
     * @param refGenomeFile		tab-delimited file containing taxonomic IDs and reference genome file names
     *
     * @throws IOException
     */
    public FileRefGenomeComputer(File refGenomeFile) throws IOException {
        this.refMap = new HashMap<String, Genome>();
        log.info("Reading reference genomes from {}.", refGenomeFile);
        // Get the directory containing the reference genome file.  The file references are relative to this.
        File baseDir = refGenomeFile.getAbsoluteFile().getParentFile();
        // This is a headerless file with two columns.
        try (TabbedLineReader refFile = new TabbedLineReader(refGenomeFile, 2)) {
            for (TabbedLineReader.Line line : refFile) {
                String taxId = line.get(0);
                // Only proceed if this is NOT a blank line.
                if (! taxId.isEmpty()) {
                    File gtoFile = new File(baseDir, line.get(1));
                    Genome gto = new Genome(gtoFile);
                    this.refMap.put(taxId, gto);
                }
            }
        }
        log.info("Reference genomes provided for {} taxonomic groups.", this.refMap.size());
    }

    @Override
    protected void initialize(GenomeStats[] reports) {
        // Loop through the genomes, computing the best GTO.
        for (GenomeStats report : reports) {
            // Only process reports with data in them.
            if (report != null) {
                Genome genome = report.getGenome();
                Iterator<TaxItem> taxonomy = genome.taxonomy();
                Genome refGenome = null;
                // Find the reference genome for the smallest taxonomic group.
                while (taxonomy.hasNext() && refGenome == null) {
                    String taxId = taxonomy.next().getId();
                    refGenome = this.refMap.get(taxId);
                }
                if (refGenome != null)
                    this.put(genome.getId(), refGenome);
            }
        }
     }

}
