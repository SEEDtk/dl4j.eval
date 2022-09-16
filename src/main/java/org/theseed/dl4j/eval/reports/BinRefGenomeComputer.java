/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import org.apache.commons.lang3.StringUtils;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.io.LineReader;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonException;
import com.github.cliftonlabs.json_simple.JsonKey;
import com.github.cliftonlabs.json_simple.JsonObject;
import com.github.cliftonlabs.json_simple.Jsoner;

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
    /** binning sample directory */
    private File binDir;


    /** reference map */

    /**
     * This object contains information about a bin's reference genome relationship, including the genome ID, the
     * number of contigs in the bin, and its mean coverage.
     */
    private static class BinData {

        /** mean coverage of the bin */
        private double coverage;
        /** reference genome ID */
        private String refGenomeId;
        /** contig count */
        private int contigCount;
        /** null list for JSON values */
        private static final JsonArray noEntries = new JsonArray();

        public static enum BinKeys implements JsonKey {
            REF_GENOMES("refGenomes", noEntries),
            CONTIGS("contigs", noEntries);

            private String key;
            private Object defaultValue;

            /**
             * Construct a bin key, saving the key string and the default value for the key.
             *
             * @param key				key string
             * @param defaultValue		value to use if the key is absent
             */
            private BinKeys(String key, Object defaultValue) {
                this.key = key;
                this.defaultValue = defaultValue;
            }

            @Override
            public String getKey() {
                return this.key;
            }

            @Override
            public Object getValue() {
                return this.defaultValue;
            }

        }

        /**
         * Construct a bin descriptor from the JSON description of the bin.
         *
         * @param json		JSON description of the bin from the bins.json file
         */
        public BinData(JsonObject json) {
            // Start with the reference genome.  This is the first genome in the reference genome list.
            JsonArray refs = json.getCollectionOrDefault(BinKeys.REF_GENOMES);
            if (refs.size() > 0)
                this.refGenomeId = refs.getString(0);
            else
                this.refGenomeId = null;
            // Now loop through the contigs, computing the coverage and count.
            JsonArray contigs = json.getCollectionOrDefault(BinKeys.CONTIGS);
            final int n = contigs.size();
            double covg = 0.0;
            long length = 0;
            for (int i = 0; i < n; i++) {
                JsonArray contigData = contigs.getCollection(i);
                int contigLen = contigData.getInteger(1);
                double contigCovg = contigData.getDouble(2);
                length += contigLen;
                covg += contigLen * contigCovg;
            }
            // Save the contig count and the mean coverage per base pair.
            this.contigCount = n;
            this.coverage = covg / length;
        }

        /**
         * @return the mean coverage for the bin
         */
        public double getCoverage() {
            return this.coverage;
        }

        /**
         * @return the ID of the bin's reference genome
         */
        public String getRefGenomeId() {
            return this.refGenomeId;
        }

        /**
         * @return the number of contigs in the bin
         */
        public int getContigCount() {
            return this.contigCount;
        }

    }

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
        // Use the bins.json file to compute the reference genomes.
        File binsJsonFile = new File(this.binDir, "bins.json");
        BinData[] binMap;
        try {
            binMap = this.readBinData(binsJsonFile);
        } catch (IOException | JsonException e1) {
            // Convert the common exceptions to runtime.  They are fatal.
            throw new RuntimeException("Error reading bins.", e1);
        }
        if (binMap.length < reports.length)
            throw new IllegalArgumentException(String.format("%d bins presented, but only %d found in bins.json.",
                    reports.length, binMap.length));
        // Loop through the binned genomes.
        for (int i = 0; i < reports.length; i++) {
            GenomeStats gReport = reports[i];
            if (gReport != null) {
                Genome genome = gReport.getGenome();
                // Get the descriptor for this bin.
                BinData binData = binMap[i];
                // Get the reference genome ID.
                String refGenomeId = binData.getRefGenomeId();
                if (refGenomeId != null) {
                    // Store the reference genome data in the bin.
                    this.storeBinData(genome, binData);
                    // Store the coverage in the genome report.
                    gReport.setBinCoverage(binData.getCoverage());
                    // Compute the file containing the reference genome.
                    File refGenomeFile = new File(this.binDir, refGenomeId + ".json");
                    if (! refGenomeFile.canRead())
                        log.warn("Reference genome {} for {} not found in {}.", refGenomeId, genome, this.binDir);
                    else {
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

    /**
     * Read the bins.json file and create an array of bin descriptors.
     *
     * @param binsJsonFile	name of the bins.json file
     *
     * @return an array of the bin descriptors, in order
     *
     * @throws IOException
     * @throws JsonException
     */
    private BinData[] readBinData(File binsJsonFile) throws IOException, JsonException {
        var retVal = new ArrayList<BinData>();
        try (LineReader jsonReader = new LineReader(binsJsonFile)) {
            // Each bin is represented by a json string, with a "//" record separating the strings.
            while (jsonReader.hasNext()) {
                JsonObject json = this.getJsonSection(jsonReader);
                retVal.add(new BinData(json));
            }
        }
        return retVal.stream().toArray(BinData[]::new);
    }

    /**
     * Here we read the current section of the file and return it as a JSON object.
     *
     * @param jsonReader	input file stream
     *
     * @return a JSON object for the current section.
     *
     * @throws JsonException
     */
    private JsonObject getJsonSection(LineReader jsonReader) throws JsonException {
        boolean done = false;
        var lines = new ArrayList<String>();
        while (jsonReader.hasNext() && ! done) {
            String line = jsonReader.next();
            if (line.contentEquals("//"))
                done = true;
            else
                lines.add(line);
        }
        // Now convert the strings into a JsonObject.
        String serial = StringUtils.join(lines, " ");
        JsonObject retVal = (JsonObject) Jsoner.deserialize(serial);
        return retVal;
    }

    /**
     * Store the reference-genome ID in the genome quality data.  This insures it is available after the bin is saved
     * in the output GTO.
     *
     * @param genome		source genome
     * @param binData		bin descriptor
     */
    private void storeBinData(Genome genome, BinData binData) {
        if (genome.getContigCount() != binData.getContigCount())
            log.warn("Genome {} does not have the same number of contigs as described in the bins.json.");
        JsonObject quality = genome.getQuality();
        quality.put("bin_ref_genome", binData.getRefGenomeId());
    }


}
