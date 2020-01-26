/**
 *
 */
package org.theseed.dl4j.eval;

import java.util.Arrays;
import java.util.Collection;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;

/**
 * This is a simple utility class that is used to contain quality statistics for a genome.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeStats {

    // FIELDS

    /** genome ID */
    private String id;
    /** completeness group */
    private String group;
    /** number of extra role occurrences for completeness */
    private int contaminationCount;
    /** number of missing roles for completeness */
    private int missingCount;
    /** total number of completeness roles */
    private int completeCount;
    /** coarse-consistent roles */
    private int coarse;
    /** fine-consistent roles */
    private int fine;
    /** total consistency roles */
    private int consistentCount;
    /** total pegs */
    private int pegCount;
    /** total hypothetical pegs */
    private int hypoCount;
    /** number of features with local protein families */
    private int plfamCount;
    /** hash of problematic roles */
    private SortedMap<String, ProblematicRole> problematicRoles;
    /** total contig length */
    private int dnaSize;
    /** number of contigs */
    private int contigCount;
    /** smallest number of contigs that covers 50% of the total */
    private int l50;
    /** smallest contig in the set of longest contigs covering 50% of the total */
    private int n50;
    /** number of seed proteins found */
    private int seedCount;
    /** longest seed protein found */
    private String seedProt;
    /** domain of this genome */
    private String domain;

    /**
     * This object describes the statistics of a problematic role.
     */
    public class ProblematicRole {
        /** TRUE if this role is problematic because of completeness */
        private boolean universal;
        /** number of predicted occurrences */
        private int predicted;
        /** number of actual occurrences */
        private int actual;

        /**
         * Construct a problematic role record.
         *
         * @param universal		TRUE if this role relates to completeness/contamination
         * @param predicted		number of predicted occurrences
         * @param actual		number of actual occurrences
         */
        private ProblematicRole(boolean universal, int predicted, int actual) {
            this.universal = universal;
            this.predicted = predicted;
            this.actual = actual;
        }

        /**
         * @return TRUE if this role is problematic because of completeness
         */
        public boolean isUniversal() {
            return universal;
        }

        /**
         * @return the number of predicted occurrences
         */
        public int getPredicted() {
            return predicted;
        }

        /**
         * @return the number of actual occurrences
         */
        public int getActual() {
            return actual;
        }


    }
    /** Construct an empty genome stats object.
     *
     * @param id	ID of the genome
     */
    public GenomeStats(String id, String domain) {
        this.id = id;
        this.group = null;
        this.contaminationCount = 0;
        this.missingCount = 0;
        this.completeCount = 0;
        this.coarse = 0;
        this.fine = 0;
        this.consistentCount = 0;
        this.pegCount = 0;
        this.hypoCount = 0;
        this.plfamCount = 0;
        this.problematicRoles = new TreeMap<String, ProblematicRole>();
        this.seedCount = 0;
        this.seedProt = "";
        this.l50 = 0;
        this.n50 = 0;
        this.dnaSize = 0;
        this.contigCount = 0;
        this.domain = domain;
    }

    /**
     * @return the genome ID
     */
    public String getId() {
        return id;
    }

    /**
     * @return the completeness group
     */
    public String getGroup() {
        return group;
    }

    /**
     * @param completeness group for this genome
     */
    public void setGroup(String group) {
        this.group = group;
    }

    /**
     * Record a role result from the consistency check.
     *
     * @param role			ID of the role being checked
     * @param predicted		number of occurrences predicted
     * @param actual		actual number of occurrences
     */
    public void consistentRole(String role, int predicted, int actual) {
        this.consistentCount++;
        if (predicted == actual) {
            this.coarse++;
            this.fine++;
        } else {
            if ((predicted == 0) == (actual == 0))
                this.coarse++;
            ProblematicRole report = new ProblematicRole(false, predicted, actual);
            this.problematicRoles.put(role, report);
        }
    }

    /**
     * Record a role result from the completeness check.  Here the prediction is always 1.
     *
     *	@param role		ID of the role in question
     *  @param actual	number of actual occurrences of the role
     */
    public void completeRole(String role, int actual) {
        this.completeCount++;
        if (actual != 1) {
            ProblematicRole report = new ProblematicRole(true, 1, actual);
            this.problematicRoles.put(role, report);
            if (actual == 0) {
                this.missingCount++;
            } else {
                this.contaminationCount += actual - 1;
            }
        }
    }

    /**
     * Record a seed protein.
     *
     *  @param prot		seed protein string to record
     */
    public void countSeed(String prot) {
        this.seedCount++;
        if (prot.length() > this.seedProt.length())
            this.seedProt = prot;
    }

    /**
     * Compute the contig metrics for this genome.
     *
     * @param genome	genome whose contig metrics are desired
     */
    public void computeMetrics(Genome genome) {
        // Get the list of contig lengths and the total DNA size.
        Collection<Contig> contigs = genome.getContigs();
        this.contigCount = contigs.size();
        int[] lengths = new int[this.contigCount];
        int i = 0;
        this.dnaSize = 0;
        for (Contig contig : contigs) {
            lengths[i] = contig.length();
            this.dnaSize += lengths[i];
            i++;
        }
        // Sort in ascending order.
        Arrays.sort(lengths);
        // Handle the no-contigs case.
        if (this.dnaSize <= 0) {
            this.l50 = 0;
            this.n50 = 0;
        } else {
            // Work backwards until we find the split point.  First we do a cheap divide by 2 rounding up.
            int remaining = (this.dnaSize + 1) >> 1;
            i = lengths.length - 1;
            for (i = lengths.length - 1; remaining > 0; i--)
                remaining -= lengths[i];
            // Get back the correct array index for the split point.
            i++;
            // Store the results.
            this.l50 = lengths.length - i;
            this.n50 = lengths[i];
        }
    }

    /**
     * Record the non-role-related information about a protein-encoding gene.
     *
     * @param peg	feature to record
     */
    public void countPeg(Feature peg) {
        this.pegCount++;
        String function = peg.getFunction();
        if (function == null)
            this.hypoCount++;
        else {
            String normalized = StringUtils.substringBefore(function, "#").trim().toLowerCase();
            if (normalized.contains("hypothetical") || normalized.isEmpty())
                this.hypoCount++;
        }
        if (peg.getPlfam() != null)
            this.plfamCount++;
    }

    /**
     * @return the coarse consistency percentage
     */
    public double getCoarsePercent() {
        double retVal = 0;
        if (consistentCount > 0)
            retVal = coarse * 100.0 / consistentCount;
        return retVal;
    }

    /**
     * @return the fine consistency percentage
     */
    public double getFinePercent() {
        double retVal = 0;
        if (consistentCount > 0)
            retVal = fine * 100.0 / consistentCount;
        return retVal;
    }

    /**
     * @return the percent complete
     */
    public double getCompletePercent() {
        double retVal = 0;
        if (completeCount > 0)
            retVal = (completeCount - missingCount) * 100.0 / completeCount;
        return retVal;
    }

    /**
     * @return the percent contamination
     */
    public double getContaminationPercent() {
        double retVal = 100;
        if (completeCount > 0)
            retVal = contaminationCount * 100.0 / (contaminationCount + completeCount);
        return retVal;
    }

    /**
     * @return the percent hypothetical
     */
    public double getHypotheticalPercent() {
        double retVal = 100;
        if (pegCount > 0)
            retVal = hypoCount * 100.0 / pegCount;
        return retVal;
    }

    /**
     * @return the list of potentially problematic roles, in sorted order by ID
     */
    public Set<String> getProblematicRoles() {
        return this.problematicRoles.keySet();
    }

    /**
     * @return the problematic role report for the specified role (if any)
     */
    public ProblematicRole getReport(String role) {
        return this.problematicRoles.get(role);
    }

    /**
     * @return the percent of proteins in local families
     */
    public double getPlfamPercent() {
        double retVal = 0;
        if (pegCount > 0)
            retVal = plfamCount * 100.0 / pegCount;
        return retVal;
    }

    /**
     * @return the seed protein string
     */
    public String getSeed() {
        return this.seedProt;
    }

    /**
     * @return TRUE if the seed protein is good, else FALSE
     */
    public boolean isGoodSeed() {
        boolean retVal = false;
        if (this.seedCount == 1) {
            int len = this.seedProt.length();
            if (this.domain.charAt(0) == 'A') {
                retVal = (len >= 293 && len <= 652);
            } else {
                retVal = (len >= 209 && len <= 405);
            }
        }
        return retVal;
    }

    /**
     * @return the contig L50 metric
     */
    public int getL50() {
        return this.l50;
    }

    /**
     * @return the contig N50 metric
     */
    public int getN50() {
        return this.n50;
    }

    /**
     * @return the number of base pairs in the genome
     */
    public int getDnaSize() {
        return dnaSize;
    }

    /**
     * @return the number of contigs in the genome
     */
    public int getContigCount() {
        return contigCount;
    }

    /**
     * @return the percent of the genome covered by peg features
     */
    public double getCdsPercent() {
        double retVal = 0;
        if (this.dnaSize > 0)
            retVal = (this.pegCount * 100000.0) / this.dnaSize;
        return retVal;
    }

    /**
     * @return the number of protein-encoding features in the genome
     */
    public int getPegCount() {
        return pegCount;
    }

    /**
     * @return the number of hypothetical proteins in the genome
     */
    public int getHypoCount() {
        return hypoCount;
    }

    /**
     * @return the number of proteins with local families in the genome
     */
    public int getPlfamCount() {
        return plfamCount;
    }

    /**
     * @return the total quality score for this genome
     */
    public double getScore() {
        return getFinePercent() * 1.09 + getCompletePercent() - 5 * getContaminationPercent() + (100 - getHypotheticalPercent()) - getContigCount() / 100.0;
    }


}
