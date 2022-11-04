/**
 *
 */
package org.theseed.dl4j.eval.stats;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.kmers.KmerCollectionGroup;
import org.theseed.reports.NaturalSort;

import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * A genome analysis contains information about a comparison between a genome and its reference genome.
 *
 * It is built from a GenomeStats object, and analyzes each problematic role to determine which features are
 * good and bad.  In addition, the SSU rRNA is considered a good feature, and any feature which is the sole
 * instance of a subsystem role is considered good.
 *
 * A feature is good if its role is useful but not problematic, it is close to a same-role protein in the
 * reference genome, or it is an SSU rRNA.  We track good and bad roles in contigs, and this information is
 * used to eliminate bad contigs.  In addition, the analysis itself is used to create explanatory text
 * about problematic roles.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeAnalysis {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeAnalysis.class);
    /** map of problematic roles to feature kmer databases */
    private Map<String, KmerCollectionGroup> featureFinder;
    /** reference genome ID */
    private String refGenomeId;
    /** reference genome object */
    private Genome refGenomeObj;
    /** feature ID to status map */
    private Map<String, FeatureAnalysis> featureMap;
    /** contig ID to analysis map */
    private Map<String, ContigAnalysis> contigMap;
    /** number of good contigs */
    private int goodContigs;
    /** number of subsystem-locked contigs */
    private int subContigs;
    /** list of bad contigs */
    private Set<Contig> badContigs;
    /** maximum distance between close features */
    private final double maxDist;
    /** genome quality report */
    private final GenomeStats gReport;
    /** role definition map */
    private final RoleMap roleMap;
    /** margin of error on split features; used to compute the maximum size of the sum of the pieces */
    private final static double SPLIT_RATIO = 1.3;


    /**
     * This class describes the characteristics of a feature.  The feature can either be good or bad.
     * A close feature performing the same role in the reference genome may be linked.  In addition,
     * it could be marked as starting or ending near the edge of a contig.
     */
    public class FeatureAnalysis {

        /** associated feature */
        private Feature feat;
        /** status (good, bad, unknown) */
        private FeatureStatus status;
        /** ID of the closest reference feature (if any) */
        private String refId;
        /** TRUE if it starts near the edge of the contig */
        private boolean startsEdge;
        /** TRUE if it ends near the edge of the contig */
        private boolean endsEdge;

        /**
         * Construct a default-status feature analysis object.
         *
         * @param inFeat	feature to analyze
         * @param status	proposed feature status
         * @param refId		ID of the closest feature in the reference genome (or NULL if there is none)
         * @param contig	contig containing the feature
         */
        protected FeatureAnalysis(Feature inFeat, FeatureStatus status, String refId, Contig contig) {
            this.feat = inFeat;
            this.status = status;
            this.refId = refId;
            // Check the contig edges.
            Location loc = inFeat.getLocation();
            this.startsEdge = loc.isEdgeStart(contig);
            this.endsEdge = loc.isEdgeEnd(contig);
        }

        /**
         * @return the ID of the analyzed feature
         */
        public String getFid() {
            return this.feat.getId();
        }

        /**
         * @return the status (good/bad/unknown)
         */
        public FeatureStatus getStatus() {
            return this.status;
        }

        /**
         * @return the closest reference-genome feature's ID (if any)
         */
        public String getRefId() {
            return this.refId;
        }

        /**
         * @return TRUE if this feature starts near the edge of the contig
         */
        public boolean isStartsEdge() {
            return this.startsEdge;
        }

        /**
         * @return TRUE if this feature ends near the edge of the contig
         */
        public boolean isEndsEdge() {
            return this.endsEdge;
        }

        /**
         * @return the location of this feature
         */
        public Location getLocation() {
            return this.feat.getLocation();
        }

    }

    /**
     * Construct a genome analysis object.
     *
     * @param gReport		genome report from quality run
     * @param refGenome		optional reference genome
     * @param roleMap		role definitions
     * @param maxDist		maximum distance for a close feature
     */
    public GenomeAnalysis(GenomeStats gReport, Genome refGenome, RoleMap roleMap, double maxDist) {
        // Save the genome report, the role map, and the distance.
        this.roleMap = roleMap;
        this.maxDist = maxDist;
        this.gReport = gReport;
        // Create the contig map.
        Genome sourceGenome = this.gReport.getGenome();
        Collection<Contig> contigs = sourceGenome.getContigs();
        this.contigMap = new HashMap<>(contigs.size() * 3 / 2 + 1);
        // Initialize the contig analysis objects.
        contigs.stream().forEach(x -> this.contigMap.put(x.getId(), new ContigAnalysis(x, gReport)));
        // Memorize the reference genome and create the feature finder.
        this.refGenomeObj = refGenome;
        if (refGenome == null) {
            // No reference genome.  Create an empty feature finder.
            this.featureFinder = new HashMap<>(5);
            this.refGenomeId = null;
        } else {
            // Here we have a reference genome and we must analyze.
            log.info("Analyzing reference genome {}.", refGenome);
            this.refGenomeId = refGenome.getId();
            // Create the problematic role directory.
            Set<String> targetRoles = gReport.getProblematicRoles();
            this.featureFinder = new HashMap<String, KmerCollectionGroup>(targetRoles.size() * 3 / 2 + 1);
            targetRoles.stream().forEach(x -> this.featureFinder.put(x, new KmerCollectionGroup()));
            // Loop through the features, putting them in the directory.  Note we only save the
            // features with problematic roles that have valid protein translations.  If our reference
            // genome does not have proteins in it, these means no close proteins will be found and
            // every problematic role's features will stay bad.
            for (Feature feat : this.refGenomeObj.getPegs()) {
                Collection<Role> roles = feat.getUsefulRoles(roleMap);
                for (Role role : roles) {
                    KmerCollectionGroup roleGroup = this.featureFinder.get(role.getId());
                    if (roleGroup != null) {
                        String prot = feat.getProteinTranslation();
                        if (! StringUtils.isEmpty(prot))
                            roleGroup.addSequence(feat.getProteinTranslation(), feat.getId());
                    }
                }
            }
        }
        // Now we create the feature map, which maps feature IDs to feature analysis objects.  We only
        // keep a feature if it has a problematic role.  Note that we mark the SSU rRNA as a good feature.
        // Good and bad features are counted in the contigs.
        this.featureMap = new HashMap<String, FeatureAnalysis>(sourceGenome.getFeatureCount() * 3 / 2 + 1);
        for (Feature feat : sourceGenome.getFeatures()) {
            // Get the contig analysis.
            ContigAnalysis contig = this.contigMap.get(feat.getLocation().getContigId());
            // Is this feature an rRNA? We need to check for the SSU rRNA, which is always a good feature.
            String fid = feat.getId();
            String type = feat.getType();
            switch (type) {
            case "rna" :
                // Here we have an RNA.  If it is an SSU rRNA or is subsystem-locked, record it in the contig
                // as a good feature.
                if (Genome.isSSURole(feat))
                    contig.countFeature(feat,  FeatureStatus.GOOD);
                else if (feat.isSubsystemLocked())
                    contig.countSubLock();
                break;
            case "CDS" :
            case "peg" :
                // Here we have a protein.  We need to find if it has a problematic role.  Note that if
                // it does, the checkProblematicRoles method adds it to the PPR list for the role.
                Collection<Role> roles = feat.getUsefulRoles(this.roleMap);
                FeatureStatus status = gReport.checkProblematicRoles(feat, roles);
                // The ID of the closest reference-genome protein will be stored in here.
                String refId = null;
                // If the feature has a problematic role, the status will have come back as bad.
                if (status == FeatureStatus.BAD) {
                    double bestDist = this.maxDist;
                    String prot = feat.getProteinTranslation();
                    for (Role role : roles) {
                        // Get the features in the reference genome with this role.
                        KmerCollectionGroup roleFeatures = this.featureFinder.get(role.getId());
                        if (roleFeatures != null) {
                            // Find the closest feature.
                            KmerCollectionGroup.Result match = roleFeatures.getBest(prot);
                            if (match.getDistance() <= bestDist) {
                                refId = match.getGroup();
                                bestDist = match.getDistance();
                                status = FeatureStatus.GOOD;
                            }
                        }
                    }
                }
                // Count the feature as a peg.
                contig.countPeg();
                // Add the feature to the feature map and count its status in the contig.
                this.featureMap.put(fid, new FeatureAnalysis(feat, status, refId, contig.getContig()));
                if (status != FeatureStatus.GOOD && feat.isSubsystemLocked())
                    contig.countSubLock();
                else
                    contig.countFeature(feat, status);
                break;
            default :
                // All other features, check for subsystem lock.
                if (feat.isSubsystemLocked())
                    contig.countSubLock();
            }
        }
        // With all the features processed, we now have counts for all the contigs.  Process them here.
        // A contig is bad if it has no good or subsystem-locked features.
        this.badContigs = new TreeSet<Contig>();
        this.goodContigs = 0;
        this.subContigs = 0;
        for (ContigAnalysis contig : this.contigMap.values()) {
            if (contig.getGoodCount() > 0)
                this.goodContigs++;
            else if (contig.getSubCount() > 0)
                this.subContigs++;
            else
                this.badContigs.add(contig.getContig());
        }
    }

     /**
     * @return the ID of the reference genome (or NULL if there is none)
     */
    public String getRefGenomeId() {
        return this.refGenomeId;
    }

    /**
     * @return the analysis for the specified contig
     *
     * @param contigId		ID of the target contig
     */
    public ContigAnalysis getContigData(String contigId) {
        return this.contigMap.get(contigId);
    }

    /**
     * @return the analysis for a source genome feature
     *
     * @param fid		ID of the target feature
     */
    public FeatureAnalysis getFeatureData(String fid) {
        return this.featureMap.get(fid);
    }

    /**
     * Compute quality-analysis comments about the role implemented using the features listed.
     *
     * @param roleId		ID of this role
     * @param features		list of features performing the specified role in the source genome
     * @param predicted		number of features predicted for the role
     *
     * @return a possibly-empty list of comments about the role implemented in this genome
     */
    public List<DomContent> computeRoleComment(String roleId, List<Feature> features, int predicted) {
        List<DomContent> retVal = new ArrayList<>();
        // The first job is to sort the implementing features by closest reference-genome feature,
        // which we do using a map.  If there is no reference genome, all the features will be
        // collected as orphans.
        var featuresByRef = new TreeMap<String, List<FeatureAnalysis>>(new NaturalSort());
        var orphans = new ArrayList<FeatureAnalysis>(features.size());
        for (Feature feat : features) {
            FeatureAnalysis fData = this.getFeatureData(feat.getId());
            String refId = fData.getRefId();
            if (refId == null)
                orphans.add(fData);
            else {
                List<FeatureAnalysis> list = featuresByRef.computeIfAbsent(refId, x -> new ArrayList<FeatureAnalysis>(features.size()));
                list.add(fData);
            }
        }
        // Now loop through the reference features found.
        for (var refData : featuresByRef.entrySet()) {
            String refId = refData.getKey();
            List<FeatureAnalysis> featList = refData.getValue();
            // First we check for a split.  We have a split if there are two matching features both are shorter than
            // the reference genome feature, at least one starts near an edge, and the other ends near an edge.
            boolean isSplit = (featList.size() == 2 ? this.checkSplit(refId, featList.get(0), featList.get(1)) : false);
            if (isSplit) {
                // Here we have a split protein, so we use a special message.
                retVal.add(this.splitMessage(refId, featList.get(0), featList.get(1)));
            } else {
                // Not a split, so use the usual messages.
                retVal.addAll(featList.stream().map(x -> this.refMessage(refId, x)).collect(Collectors.toList()));
            }
        }
        // Next, process the orphans.
        for (FeatureAnalysis feat : orphans)
            retVal.add(this.featMessage(feat));
        // If we have a reference genome, add a message about how the role is implemented.
        if (this.refGenomeObj != null)
            retVal.add(this.roleRefMessage(roleId, features.size(), predicted));
        return retVal;
    }

    /**
     * @return a basic HTML message about a feature in the source genome
     *
     * @param feat	analysis of the feature of interest
     */
    private DomContent featMessage(FeatureAnalysis feat) {
        Genome source = this.gReport.getGenome();
        // Get the location comment for the feature.
        Location loc = feat.getLocation();
        ContigAnalysis contig = this.getContigData(loc.getContigId());
        DomContent locComment = contig.locationComment(loc);
        // Form the full sentence using the feature link.
        DomContent retVal = join(source.featureRegionLink(feat.getFid()), locComment);
        return retVal;
    }

    /**
     * Add a message about how this role is implemented in the reference genome.
     *
     * @param roleId		ID of the role in question
     * @param actual		actual number of role implementations in the source genome
     * @param predicted		predicted number of role implementations in the reference genome.
     *
     * @return an HTML message about how often the role occurs compared to what is actually there
     */
    private DomContent roleRefMessage(String roleId, int actual, int predicted) {
        // We have two phrases:  the declaration and an optional codicil about the predicted/actual.
        DomContent retVal;
        PhraseBuilder codicil = new PhraseBuilder(", which", ".");
        // First, find out if this role occurs at all.
        var refKmers = this.featureFinder.get(roleId);
        int refCount;
        if (refKmers == null || refKmers.size() == 0) {
            // Here it is missing in the reference genome.
            refCount = 0;
            retVal = text("This role is missing in the reference genome");
        } else {
            // Here it is implemented in the reference genome.
            refCount = refKmers.size();
            String verb = (refCount == 1 ? "implements" : "implement");
            retVal = join(this.refGenomeObj.featureListLink(refKmers.getKeys()), verb + " this role in the reference genome");
        }
        // Check for a match to the actual or predicted count.
        if (refCount == actual)
            codicil.addPhrase("matches this genome");
        else if (refCount == predicted)
            codicil.addPhrase("is expected from the prediction");
        // Form the final comment.
        retVal = join(retVal, codicil.toString());
        return retVal;
    }

    /**
     * Create a message describing two features that appear to be a single split protein.
     *
     * @param refId		ID of the reference feature
     * @param feat1		first split feature
     * @param feat2		second split feature
     *
     * @return an HTML message describing the split
     */
    private DomContent splitMessage(String refId, FeatureAnalysis feat1, FeatureAnalysis feat2) {
        Genome source = this.gReport.getGenome();
        // Get the IDs of the first and second features in the split.
        String fid1;
        String fid2;
        if (feat1.isEndsEdge() && feat2.isStartsEdge()) {
            fid1 = feat2.getFid();
            fid2 = feat1.getFid();
        } else {
            fid1 = feat1.getFid();
            fid2 = feat2.getFid();
        }
        // Build the message.
        DomContent retVal = join(source.featureRegionLink(fid1), "and", source.featureRegionLink(fid2),
                "appear to be parts of", this.refGenomeObj.featureLink(refId), "split across a contig boundary.");
        return retVal;
    }

    /**
     * Create a message about a feature that is implemented in the reference genome.
     *
     * @param refId		ID of the reference feature
     * @param desc		feature descriptor for the source genome feature
     *
     * @return an HTML message about the source feature and its closest reference
     */
    private DomContent refMessage(String refId, FeatureAnalysis desc) {
        // Start with the base message.
        DomContent base = this.featMessage(desc);
        // Add the reference data.
        DomContent retVal = join(base, "It is close to", this.refGenomeObj.featureLink(refId),
                ", which performs this role in the reference genome.");
        return retVal;
    }

    /**
     * @return TRUE if the two features are likely a split implementation of the specified reference feature
     *
     * @param refId		ID of the reference feature matched
     * @param feat1		analysis of the first feature
     * @param feat2		analysis of the second feature
     */
    private boolean checkSplit(String refId, FeatureAnalysis feat1, FeatureAnalysis feat2) {
        boolean retVal = false;
        // Verify we are straddling a contig boundary.
        if (feat1.isStartsEdge() && feat2.isEndsEdge() || feat1.isEndsEdge() && feat2.isStartsEdge()) {
            // Yes.  Now we check lengths.
            Location loc1 = feat1.getLocation();
            Location loc2 = feat2.getLocation();
            Feature refFeat = this.refGenomeObj.getFeature(refId);
            int refLen = refFeat.getLocation().getLength();
            // If the pieces look lik
            if (loc1.getLength() + loc2.getLength() <= refLen * SPLIT_RATIO)
                retVal = true;
        }
        return retVal;
    }

    /**
     * @return the reference genome object.
     */
    public Genome getRefGenome() {
        return this.refGenomeObj;
    }

    /**
     * @return the number of good contigs
     */
    public int getGoodContigCount() {
        return this.goodContigs;
    }

    /**
     * @return the number of subsystem-locked contigs
     */
    public int getSubContigCount() {
        return this.subContigs;
    }

    /**
     * @return the set of bad contigs
     */
    public Set<Contig> getBadContigs() {
        return this.badContigs;
    }

    /**
     * @return TRUE if this analysis used a reference genome, else FALSE
     */
    public boolean hasRefGenome() {
        return this.refGenomeObj != null;
    }

    /**
     * @return TRUE if the genome is good, else FALSE
     */
    public boolean isGood() {
        return this.gReport.isGood();
    }

}
