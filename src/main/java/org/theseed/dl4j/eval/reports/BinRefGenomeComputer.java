/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.QualityKeys;
import org.theseed.bins.Bin;
import org.theseed.bins.BinGroup;

import com.github.cliftonlabs.json_simple.JsonException;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This object computes the reference genomes for a binning job.  For each binning directory, the
 * "setSampleDir" method should be called before producing reports.  The reference genome will be
 * loaded from the appropriate json file in the binning directory.  The bins.json file is used to
 * associate the bins with the reference genomes.
 *
 * @author Bruce Parrello
 *
 */
public class BinRefGenomeComputer extends RefGenomeComputer {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(BinRefGenomeComputer.class);
    /** binning sample directory */
    private File binDir;
    /** genome directory for reference genomes */
    private GenomeDirectory refGenomes;

     /**
     * Specify the binning directory containing the sample's files.
     *
     * @param binDir	current binning directory
     */
    public void setSampleDir(File binDir) {
        this.binDir = binDir;
        File refDir = new File(binDir, "RefGenomes");
        try {
            this.refGenomes = new GenomeDirectory(refDir);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    protected void initialize(GenomeStats[] reports) {
        try {
            // Use the bins.json file to compute the reference genomes.
            File binsJsonFile = new File(this.binDir, "bins.json");
            BinGroup binGroup = new BinGroup(binsJsonFile);
            var realBins = binGroup.getSignificantBins();
            if (realBins.size() < reports.length)
                throw new IllegalArgumentException(String.format("%d bins presented, but only %d found in bins.json.",
                        reports.length, realBins.size()));
            // Loop through the binned genomes.
            for (GenomeStats gReport : reports) {
                if (gReport != null) {
                    Genome genome = gReport.getGenome();
                    // Get the descriptor for this bin.
                    int speciesId = genome.getTaxonomyId();
                    // Get the reference genome ID.
                    Bin sourceBin = realBins.stream().filter(x -> x.getTaxonID() == speciesId).findAny().orElse(null);
                    if (sourceBin == null)
                        log.warn("Cannot find bin for genome {}.", genome);
                    else {
                        String refGenomeId = sourceBin.getRefGenome();
                        // Store the reference genome data in the bin.
                        this.storeBinData(genome, sourceBin);
                        // Store the coverage in the genome report.
                        gReport.setBinCoverage(sourceBin.getCoverage());
                        // Get the reference genome from the genome source and attach it.
                        Genome refGenome = this.refGenomes.getGenome(refGenomeId);
                        this.put(genome.getId(), refGenome);
                    }
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        } catch (JsonException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Store the reference-genome ID in the genome quality data.  This insures it is available after the bin is saved
     * in the output GTO.
     *
     * @param genome		source genome
     * @param binData		bin descriptor
     */
    private void storeBinData(Genome genome, Bin binData) {
        log.info("Genome {} has {} contigs, compared to {} in the original bin.", genome, genome.getContigCount(),
                binData.getContigs().size());
        JsonObject quality = genome.getQuality();
        quality.put(QualityKeys.BIN_REF_GENOME.getKey(), binData.getRefGenome());
    }

    /**
     * Store the reference-genome ID and the coverage in the genome quality data.
     *
     * @param genome		genome to update
     * @param refId			reference genome ID
     * @param coverage		binning coverage
     */
    public static void storeBinData(Genome genome, String refId, double coverage) {
        JsonObject quality = genome.getQuality();
        if (refId != null)
            quality.put(QualityKeys.BIN_REF_GENOME.getKey(), refId);
        if (coverage != 0.0)
            quality.put(QualityKeys.BIN_COVERAGE.getKey(), coverage);
    }

    /**
     * @return the coverage data from a genome.
     *
     * @param genome	genome whose coverage is desired
     */
    public static double getBinCoverage(Genome genome) {
        JsonObject quality = genome.getQuality();
        double retVal = quality.getDoubleOrDefault(QualityKeys.BIN_COVERAGE);
        return retVal;
    }

}
