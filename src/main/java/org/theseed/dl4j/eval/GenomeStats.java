/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.UnsupportedEncodingException;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.UUID;

import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.MD5Hex;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This is a simple utility class that is used to contain quality statistics for a genome.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeStats {

    // FIELDS

    /** source genome */
    private Genome genome;
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
    /** map of problematic roles */
    private SortedMap<String, ProblematicRole> problematicRoles;
    /** set of useful roles */
    private Set<String> usefulRoles;
    /** total contig length */
    private int dnaSize;
    /** number of contigs */
    private int contigCount;
    /** smallest number of contigs that covers 50% of the total */
    private int l50;
    /** smallest contig in the set of longest contigs covering 50% of the total */
    private int n50;
    /** smallest number of contigs that covers 70% of the total */
    private int l70;
    /** smallest contig in the set of longest contigs covering 70% of the total */
    private int n70;
    /** number of seed proteins found */
    private int seedCount;
    /** longest seed protein found */
    private String seedProt;

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
        /** list of features containing this role */
        private List<Feature> features;

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
            this.features = new ArrayList<Feature>(5);
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

        /**
         * @return the list of features containing this role
         */
        public List<Feature> getFeatures() {
            return this.features;
        }

    }

    /**
     * Return structure for feature counts.
     */
    public class Counts {

        private int fidCount;
        private int pegCount;
        private int hypoCount;
        private int usefulCount;
        private int knownCount;

        /**
         * Compute the feature counts for a genome.
         *
         * @param genome			the genome of interest
         * @param roleDefinitions	role definition map
         */
        public Counts(Genome genome, RoleMap roleDefinitions) {
            this.fidCount = 0;
            this.pegCount = 0;
            this.hypoCount = 0;
            this.usefulCount = 0;
            this.knownCount = 0;
            // Run through all the features of the genome, counting.
            for (Feature feat : genome.getFeatures()) {
                this.fidCount++;
                if (feat.isProtein()) {
                    this.pegCount++;
                    if (Feature.isHypothetical(feat.getFunction())) {
                        this.hypoCount++;
                    } else {
                        // We are not hypothetical, so determine if one of our roles
                        // is known, and if one of them is useful with respect to this
                        // object's validation.
                        boolean useful = false;
                        boolean known = false;
                        for (Role role : feat.getUsefulRoles(roleDefinitions)) {
                            known = true;
                            if (GenomeStats.this.isUseful(role.getId())) {
                                useful = true;
                            }
                        }
                        // Do the last two counts.
                        if (useful) this.usefulCount++;
                        if (known) this.knownCount++;
                    }
                }
            }
        }

        /**
         * @return the number of features
         */
        public int getFidCount() {
            return fidCount;
        }

        /**
         * @return the number of protein features
         */
        public int getPegCount() {
            return pegCount;
        }

        /**
         * @return the number of hypothetical proteins
         */
        public int getHypoCount() {
            return hypoCount;
        }

        /**
         * @return the number of proteins with roles useful in evaluation
         */
        public int getUsefulCount() {
            return usefulCount;
        }

        /**
         * @return the number of proteins with roles used in subsystems
         */
        public int getKnownCount() {
            return knownCount;
        }



    }


    /**
     * Status of a feature:  GOOD, BAD, or UNKNOWN
     */
    public enum FeatureStatus {
        GOOD, BAD, UNKNOWN;
    }

    /** Construct an empty genome stats object.
     *
     * @param genome	source genome object
     */
    public GenomeStats(Genome genome) {
        this.genome = genome;
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
        this.usefulRoles = new HashSet<String>(2500);
        this.seedCount = 0;
        this.seedProt = "";
        this.l50 = 0;
        this.n50 = 0;
    }

    /**
     * @return the genome ID
     */
    public String getId() {
        return this.genome.getId();
    }

    /**
     * @return the genome name
     */
    public String getName() {
        return this.genome.getName();
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
        this.usefulRoles.add(role);
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
        this.usefulRoles.add(role);
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
            this.n70 = 0;
            this.l70 = 0;
        } else {
            // Work backwards until we find the split point for 50% of the genome.
            i = findSplitPoint(lengths, 50);
            // Store the results.
            this.l50 = lengths.length - i;
            this.n50 = lengths[i];
            // Do it again for N70.
            i = findSplitPoint(lengths, 70);
            this.l70 = lengths.length - i;
            this.n70 = lengths[i];
        }
    }

    /**
     * Given a sorted array of contig lengths, find the contig length at the split point for the
     * specified ratio.  A ratio of 70, for example, will return the contig length for the N70.
     *
     * @param lengths	array of contig lengths, sorted in ascending order
     * @param ratio		percentage for the split point
     *
     * @return the index of the contig length at the split point
     */
    private int findSplitPoint(int[] lengths, int ratio) {
        int retVal;
        // Compute the split point, rounding up.
        int remaining = (this.dnaSize * ratio + 99) / 100;
        retVal = lengths.length - 1;
        for (retVal = lengths.length - 1; remaining > 0; retVal--)
            remaining -= lengths[retVal];
        // Get back the correct array index for the split point.
        retVal++;
        return retVal;
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
            boolean hypo = Feature.isHypothetical(function);
            if (hypo)
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
            if (this.getDomain().charAt(0) == 'A') {
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

    /**
     * @return the genome object for these statistics
     */
    public Genome getGenome() {
        return this.genome;
    }

    /**
     * Store this quality information in the specified Genome's GTO.
     *
     * @param genome	genome into which the quality information should be stored
     * @param roles		role definition table
     * @param version	version string for the evaluation database
     * @param options	original command-line arguments
     *
     * @throws NoSuchAlgorithmException
     * @throws UnsupportedEncodingException
     */
    public void store(Genome genome, RoleMap roles, String version, String[] options) throws NoSuchAlgorithmException, UnsupportedEncodingException {
        // Record this as an analysis event.
        JsonObject gto = genome.toJson();
        JsonArray events = (JsonArray) gto.get("analysis_events");
        if (events == null) {
            events = new JsonArray();
            gto.put("analysis_events", events);
        }
        JsonObject thisEvent = new JsonObject().putChain("id", UUID.randomUUID().toString())
                .putChain("tool_name", "org.theseed.dl4j.eval")
                .putChain("execute_time", System.currentTimeMillis() / 1000.0);
        try {
            thisEvent.put("hostname", InetAddress.getLocalHost().getCanonicalHostName());
        } catch (UnknownHostException e) { }
        events.add(thisEvent);
        // Insure there is a quality member in the GTO and get access to it.
        JsonObject quality = (JsonObject) gto.get("quality");
        if (quality == null) {
            quality = new JsonObject();
            gto.put("quality", quality);
        }
        quality.put("genome_length", this.dnaSize);
        quality.put("contigs", this.contigCount);
        quality.put("plfam_cds_ratio", this.getPlfamPercent());
        quality.put("coarse_consistency", this.getCoarsePercent());
        quality.put("contamination", this.getContaminationPercent());
        quality.put("hypothetical_cds_ratio", this.getHypotheticalPercent());
        quality.put("completeness", this.getCompletePercent());
        quality.put("cds_ratio", this.getCdsPercent());
        quality.put("fine_consistency", this.getFinePercent());
        quality.put("completeness_group", this.getGroup());
        // Store the version string.
        quality.put("eval_version", version);
        // Genome metrics is a hash.
        JsonObject metrics = new JsonObject().putChain("N50", this.getN50()).putChain("L50", this.getL50())
                .putChain("N70", this.getN70()).putChain("L70", this.getL70())
                .putChain("totlen", this.dnaSize);
        quality.put("genome_metrics", metrics);
        // Finally, the problematic roles report.  We create two lists in parallel, one for consistency and one for completeness.
        // We also need to count the overs and unders, and we need to save the role definitions.
        JsonObject pprComplete = new JsonObject();
        JsonObject pprConsistent = new JsonObject();
        JsonObject roleMap = new JsonObject();
        int over = 0;
        int under = 0;
        for (String role : this.getProblematicRoles()) {
            ProblematicRole ppr = this.getReport(role);
            JsonObject target = (ppr.universal ? pprComplete : pprConsistent);
            JsonArray comparison = new JsonArray().addChain(ppr.predicted).addChain(ppr.actual);
            target.put(role, comparison);
            roleMap.put(role, roles.getName(role));
            if (ppr.actual > ppr.predicted)
                over++;
            else
                under++;
        }
        JsonObject problematicRoles = new JsonObject().putChain("consistency_roles", pprConsistent)
                .putChain("completeness_roles", pprComplete).putChain("role_map", roleMap);
        problematicRoles.put("consistency_checked", this.consistentCount);
        problematicRoles.put("completeness_checked", this.completeCount);
        problematicRoles.put("over_present", over);
        problematicRoles.put("under_present", under);
        quality.put("problematic_roles_report", problematicRoles);
        if (genome.hasContigs()) {
            MD5Hex md5Engine = new MD5Hex();
            quality.put("dna_md5", md5Engine.sequenceMD5(genome));
        }
    }

    /**
     * @return TRUE if this genome is clean, else FALSE
     */
    public boolean isClean() {
        return (getContaminationPercent() < 10.0);
    }

    /**
     * @return TRUE if this genome is good, else FALSE
     */
    public boolean isGood() {
        return isClean() && isUnderstood() && isComplete() && isConsistent() && isGoodSeed();
    }

    /**
     * @return TRUE if this genome is consistently annotated, else FALSE
     */
    public boolean isConsistent() {
        return (getFinePercent() >= 87.0);
    }

    /**
     * @return TRUE if this genome is mostly complete, else FALSE
     */
    public boolean isComplete() {
        return (getCompletePercent() >= 80.0);
    }

    /**
     * @return TRUE if this genome's proteins are understood, else FALSE
     */
    public boolean isUnderstood() {
        return (getHypotheticalPercent() <= 70.0);
    }

    /**
     * @return the number of roles used to check consistency
     */
    public int getConsistencyRoleCount() {
        return this.consistentCount;
    }

    /**
     * @return the number of roles used to check completeness
     */
    public int getCompletenessRoleCount() {
        return this.completeCount;
    }

    /**
     * @return the domain of this genome
     */
    public String getDomain() {
        return this.genome.getDomain();
    }

    /**
     * Process a feature to determine if it is problematic, and return a good/bad flag.
     *
     * @param feat		feature to check
     * @param roles		collection of roles in the feature
     *
     * @return the quality status of this feature
     */
    public FeatureStatus checkProblematicRoles(Feature feat, Collection<Role> roles) {
        // Loop through the roles in this genome.  We presume the feature is unknown unless
        // it has at least one good or bad role.  Good roles override the bad ones.
        FeatureStatus retVal = FeatureStatus.UNKNOWN;
        for (Role role : roles) {
            String roleId = role.getId();
            ProblematicRole ppr = this.getReport(roleId);
            if (ppr != null) {
                // Here the role is bad.
                ppr.features.add(feat);
                if (retVal != FeatureStatus.GOOD)
                    retVal = FeatureStatus.BAD;
            } else if (this.usefulRoles.contains(roleId)) {
                // Here we checked the role and it is not bad, so it is good.
                retVal = FeatureStatus.GOOD;
            }
        }
        return retVal;
    }

    /**
     * @return the contig L70 metric
     */
    public int getL70() {
        return l70;
    }

    /**
     * @return the contig N70 metric
     */
    public int getN70() {
        return n70;
    }

    /**
     * @return the number of roles checked by both methods
     */
    public int getCommonRoleCount() {
        return this.completeCount + this.consistentCount - this.usefulRoles.size();
    }

    /**
     * @return TRUE if this role was checked during evaluation
     */
    public boolean isUseful(String role) {
        return this.usefulRoles.contains(role);
    }

}
