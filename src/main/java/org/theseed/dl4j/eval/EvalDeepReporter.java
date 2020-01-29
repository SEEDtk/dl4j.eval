/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.theseed.dl4j.eval.GenomeStats.FeatureStatus;
import org.theseed.genome.Feature;
import org.theseed.p3api.Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.kmers.KmerCollectionGroup;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * This is an advanced version of the HTML reporting service that includes comparisons to a reference genome.
 *
 * @author Bruce Parrello
 *
 */
public class EvalDeepReporter extends EvalHtmlReporter {

    /** kmer database for finding the reference genome */
    private KmerCollectionGroup referenceGenomes;
    /** map of problematic roles to feature kmer databases */
    private Map<String, KmerCollectionGroup> featureFinder;
    /** reference genome ID */
    private String refGenomeId;
    /** reference genome name */
    private String refGenomeName;
    /** distance to the reference genome */
    private double refGenomeDistance;
    /** connection to PATRIC */
    private Connection p3;
    /** map of bad features to closest reference features */
    private Map<String,String> closeFeatureMap;
    /** reference genome ID override */
    private String refGenomeOverride;

    /** maximum acceptable distance for a reference genome */
    private static final double MAX_GENOME_DIST = 0.8;

    /** maximum acceptable distance for a close protein */
    private static final double MAX_PROTEIN_DIST = 0.9;

    /**
     * Construct a deep HTML reporting object.
     *
     * @param outDir	target output directory
     */
    public EvalDeepReporter(File outDir) {
        super(outDir);
        this.refGenomeOverride = null;
    }

    /**
     * @return the detail level needed in genomes read from PATRIC
     */
    public P3Genome.Details getDetailLevel() {
        // We need the proteins to do deep analysis.
        return P3Genome.Details.PROTEINS;
    }

    /**
     * Initialize the reference-genome database.
     *
     * @param modelDir	evaluation data directory
     */
    @Override
    protected void initialize(File modelDir) throws IOException {
        // Initialize the base class.
        super.initialize(modelDir);
        // Only bother to read in the reference-genome table if we are NOT in override mode.
        if (this.refGenomeOverride == null) {
            // Open the reference genome FASTA file.
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
        // Get a PATRIC connection.
        this.p3 = new Connection();
    }

    /**
     * Add the reference genome information to the detail table rows.
     *
     * @param detailRows	list of accumulating table rows
     */
    @Override
    protected void advancedDetailRows(ArrayList<DomContent> detailRows) {
        if (this.refGenomeId != null) {
            this.detailRow(detailRows, "Reference Genome", td(join(genomeLink(this.refGenomeId), this.refGenomeName)));
        }
    }

    /**
     * @return a possibly-modified status for a bad feature
     *
     * @param feat		feature of interest
     * @param roles		roles performed by the feature
     */
    @Override
    protected FeatureStatus advancedFeatureAnalysis(Feature feat, Collection<Role> roles) {
        // The default is not to modify the status.
        FeatureStatus retVal = FeatureStatus.BAD;
        // Only proceed if we have a reference genome.
        if (this.refGenomeId != null) {
            // Loop through the roles.
            for (Role role : roles) {
                String roleId = role.getId();
                // Get the features in the reference genome with this role.
                KmerCollectionGroup roleFeatures = this.featureFinder.get(roleId);
                if (roleFeatures != null) {
                    // Find the closest feature.
                    KmerCollectionGroup.Result match = roleFeatures.getBest(feat.getProteinTranslation());
                    if (match.getDistance() <= MAX_PROTEIN_DIST) {
                        // We found one.  Save it and upgrade the status.
                        this.closeFeatureMap.put(feat.getId(), match.getGroup());
                        retVal = FeatureStatus.GOOD;
                    }
                }
            }
        }
        return retVal;
    }

    /**
     * This is a stub that a subclass can use to create more advanced comments about a feature.
     *
     * @param feat		feature of interest
     * @param gReport	genome quality report
     * @param role		role of interest
     *
     * @return text describing the feature, or NULL if there is nothing of interest
     */
    @Override
    protected DomContent advancedFeatureComment(Feature feat, GenomeStats gReport, String role) {
        DomContent retVal = null;
        // Only proceed if we found a reference genome.
        if (this.refGenomeId != null) {
            // Do we have a close reference feature?
            String closeFid = this.closeFeatureMap.get(feat.getId());
            if (closeFid != null) {
                // We found one. Form our message.
                retVal = join(this.featureLink(closeFid), "in the reference genome is close, and performs the same role.");
            }
        }
        return retVal;
    }

    /**
     * This is a stub that a subclass can use to create more advanced comments about a missing role.
     *
     * @param listRows	an HTML list to which the advanced comment can be added
     * @param gReport	genome quality report
     * @param role		ID of the role of interest
     */
    @Override
    protected void advancedRoleComment(ContainerTag list, GenomeStats gReport, String role) {
        // Only proceed if a reference genome was found.
        if (this.refGenomeId != null) {
            // Get the features for this role from the kmer map.
            KmerCollectionGroup roleFeatures = this.featureFinder.get(role);
            if (roleFeatures == null || roleFeatures.size() == 0) {
                list.with(li("No features perform this role in the reference genome."));
            } else {
                Collection<String> fidList = roleFeatures.getKeys();
                list.with(li(join("This role is performed by ", featureListLink(fidList), "in the reference genome.")));
            }
        }
    }

    /**
     * Compute the reference genome and fill in the protein kmer database.
     *
     * @param gReport	quality report on the genome of interest
     */
    protected void advancedGenomeSetup(GenomeStats gReport) {
        if (this.refGenomeOverride != null) {
            // In override mode we use the same genome every time, but we need to make sure it is not the
            // genome we're evaluating.
            this.refGenomeId = this.refGenomeOverride;
            if (this.refGenomeId.contentEquals(gReport.getId())) {
                log.info("Genome {} is a reference to itself.  No useful data is available.", this.refGenomeId);
                this.refGenomeId = null;
            }
        } else {
            // Here we need to compute a reference. Get the seed protein.
            String seedProt = gReport.getSeed();
            // Find the closest genome in the reference genome database.
            KmerCollectionGroup.Result refGenomeData = this.referenceGenomes.getBest(seedProt);
            this.refGenomeDistance = refGenomeData.getDistance();
            this.refGenomeId = refGenomeData.getGroup();
            if (this.refGenomeId == null) {
                log.info("No reference genome found for {}.", gReport.getId());
            } else if (this.refGenomeDistance > MAX_GENOME_DIST) {
                log.info("Reference genome {} found, but distance of {} exceeds the maximum of {}.", this.refGenomeId, this.refGenomeDistance, MAX_GENOME_DIST);
                this.refGenomeId = null;
            } else if (this.refGenomeId.contentEquals(gReport.getId())) {
                log.info("Genome {} is a reference to itself.  No useful data is available.", this.refGenomeId);
                this.refGenomeId = null;
            } else {
                // Read in the genome.
                P3Genome refGenome = P3Genome.Load(p3, this.refGenomeId, P3Genome.Details.PROTEINS);
                if (refGenome == null) {
                    log.error("Reference genome {} not found in PATRIC.", this.refGenomeId);
                    this.refGenomeId = null;
                } else {
                    // Here we can use the reference genome.
                    this.refGenomeName = refGenome.getName();
                    log.info("Analyzing reference genome {}: {}.", this.refGenomeId, this.refGenomeName);
                    // Create the problematic role directory.
                    Set<String> targetRoles = gReport.getProblematicRoles();
                    this.featureFinder = new HashMap<String, KmerCollectionGroup>(targetRoles.size());
                    for (String role : targetRoles) {
                        featureFinder.put(role, new KmerCollectionGroup());
                    }
                    // Get the role definitions.
                    RoleMap roleMap = this.getRoleMap();
                    // Loop through the features, putting them in the directory.
                    for (Feature feat : refGenome.getPegs()) {
                        Collection<Role> roles = feat.getUsefulRoles(roleMap);
                        for (Role role : roles) {
                            KmerCollectionGroup roleGroup = this.featureFinder.get(role.getId());
                            if (roleGroup != null)
                                roleGroup.addSequence(feat.getProteinTranslation(), feat.getId());
                        }
                    }
                }
            }
        }
        // Finally, initialize the close-feature map.
        this.closeFeatureMap = new HashMap<String, String>();
    }

    /**
     * Specify an override reference genome ID.
     *
     * @param refGenomeId	ID of the reference genome to be used for all reports
     */
    public void setRefGenomeOverride(String refGenomeId) {
        this.refGenomeOverride = refGenomeId;
        // Read in the genome.
        this.p3 = new Connection();
        P3Genome refGenome = P3Genome.Load(p3, refGenomeId, P3Genome.Details.PROTEINS);
        if (refGenome == null) {
            throw new IllegalArgumentException("Reference genome " + refGenome + " not found in PATRIC.");
        } else {
            // Here we can use the reference genome.
            this.refGenomeName = refGenome.getName();
            log.info("Analyzing reference genome {}: {}.", this.refGenomeId, this.refGenomeName);
            // Create the problematic role directory.
            this.featureFinder = new HashMap<String, KmerCollectionGroup>();
            // Get the role definitions.  Note we are going to save all the roles, not just the problematic ones,
            // because we are re-using this genome.
            RoleMap roleMap = this.getRoleMap();
            // Loop through the features, putting them in the directory.
            for (Feature feat : refGenome.getPegs()) {
                Collection<Role> roles = feat.getUsefulRoles(roleMap);
                for (Role role : roles) {
                    KmerCollectionGroup roleGroup = this.featureFinder.get(role.getId());
                    if (roleGroup == null) {
                        roleGroup = new KmerCollectionGroup();
                        this.featureFinder.put(role.getId(), roleGroup);
                    }
                    roleGroup.addSequence(feat.getProteinTranslation(), feat.getId());
                }
            }
        }
    }


}
