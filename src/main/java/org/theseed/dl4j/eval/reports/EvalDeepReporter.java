/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.theseed.counters.GenomeEval;
import org.theseed.dl4j.eval.GenomeStats;
import org.theseed.dl4j.eval.GenomeStats.FeatureStatus;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.compare.CompareFeatures;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.kmers.KmerCollectionGroup;
import org.theseed.reports.Html;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * This is an advanced version of the HTML reporting service that includes comparisons to a reference genome.
 *
 * @author Bruce Parrello
 *
 */
public class EvalDeepReporter extends EvalHtmlReporter implements IRefReporter {

    /** computation object for finding the reference genome */
    private RefGenomeComputer refEngine;
    /** map of problematic roles to feature kmer databases */
    private Map<String, KmerCollectionGroup> featureFinder;
    /** reference genome ID */
    private String refGenomeId;
    /** reference genome object */
    private Genome refGenomeObj;
    /** map of bad features to closest reference features */
    private Map<String,String> closeFeatureMap;
    /** reference genome ID override */
    private double maxProtDist;
    /** genome comparator */
    private CompareFeatures compareObj;
    /** orf report */
    private List<DomContent> orfReportRows;

    /**
     * Construct a deep HTML reporting object.
     */
    public EvalDeepReporter() {
        this.maxProtDist = RefGenomeComputer.MAX_GENOME_DIST;
        try {
            this.compareObj = new CompareFeatures();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * @return the detail level needed in genomes read from PATRIC
     */
    public P3Genome.Details getDetailLevel() {
        // We need the proteins to do deep analysis.
        return P3Genome.Details.PROTEINS;
    }

    /* Add the reference genome information to the detail table rows.
     *
     * @param detailRows	list of accumulating table rows
     */
    @Override
    protected void advancedDetailRows(ArrayList<DomContent> detailRows) {
        if (this.refGenomeId != null) {
            Html.detailRow(detailRows, "Reference Genome", td(join(this.refGenomeObj.genomeLink(), refGenomeName())));
        }
    }

    /**
     * Here we set up the array to contain the ORF report rows.
     */
    @Override
    protected void initialize(File modelDir) throws IOException {
        super.initialize(modelDir);
        this.orfReportRows = new ArrayList<DomContent>();
        this.orfReportRows.add(tr(th("Genome ID"), th("Genome Name"), th("Ref Genome"), th("Subsystem Score").withClass("num"),
                th("Annotation Score").withClass("num"), th("Function Score").withClass("num"), th("False Positive Rate").withClass("num")));
    }

    /**
     * Create the comparison report for the reference genome.
     *
     * @param gReport	quality
     *
     * @return the HTML to insert in the output report
     */
    protected DomContent advancedGenomeAnalysis(GenomeStats gReport) {
        DomContent retVal = null;
        if (this.refGenomeId != null) {
            // We need to compute the feature counts for this genome and the reference.
            RoleMap roleDefinitions = this.getRoleMap();
            GenomeStats.Counts newCounts = gReport.new Counts(gReport.getGenome(), roleDefinitions);
            GenomeStats.Counts refCounts = gReport.new Counts(this.refGenomeObj, roleDefinitions);
            // Create a table of the counts.
            ArrayList<DomContent> tableRows = new ArrayList<DomContent>(5);
            tableRows.add(tr(th("Feature count"), th("This Genome").withClass("num"),
                    th("Ref Genome").withClass("num"), th("% This / Ref").withClass("num")));
            tableRows.add(compareTableRow("Total features in the genome", newCounts.getFidCount(), refCounts.getFidCount()));
            tableRows.add(compareTableRow("Features that are protein-coding genes", newCounts.getPegCount(), refCounts.getPegCount()));
            tableRows.add(compareTableRow("Features with hypothetical proteins", newCounts.getHypoCount(), refCounts.getHypoCount()));
            tableRows.add(compareTableRow("Features performing subsystem-related roles", newCounts.getKnownCount(), refCounts.getKnownCount()));
            tableRows.add(compareTableRow("Features performing one of the roles used in this evaluation",
                    newCounts.getUsefulCount(), refCounts.getUsefulCount()));
            retVal = Html.formatTable("Comparison of " + gReport.getId() + " with Reference Genome " + this.refGenomeId,
                    tableRows);
            // Check to see if we can do a comparison report.
            try {
                boolean comparable = this.compareObj.compare(gReport.getGenome(), this.refGenomeObj);
                if (comparable) {
                    // Create the ORF comparison report.
                    CompareFeatures comparison = this.compareObj;
                    DomContent report = compareReport(comparison);
                    retVal = join(retVal, report);
                    // Add a report row for this genome to the master ORF report.  We make a safety check in case the reference
                    // genome is terrible and would cause a divide-by-0 error.
                    double denominator = refCounts.getPegCount();
                    if (denominator > 0 && refCounts.getKnownCount() > 0) {
                        DomContent orfRow = tr(
                                td(gReport.getGenome().genomeLink()),
                                td(Html.gPageLink(gReport.getId(), gReport.getName())),
                                td(this.refGenomeObj.genomeLink()),
                                Html.numCell(newCounts.getKnownCount() * 100 / (double) refCounts.getKnownCount()),
                                Html.numCell(this.compareObj.getIdentical() * 100 / denominator),
                                Html.numCell((this.compareObj.getIdentical() + this.compareObj.getShorter() + this.compareObj.getLonger()) * 100 / denominator),
                                Html.numCell(this.compareObj.getNewOnlyCount() * 100 / denominator));
                        this.orfReportRows.add(orfRow);
                    }
                }
            } catch (UnsupportedEncodingException e) {
                throw new RuntimeException(e);
            }
        }
        return retVal;
    }

    /**
     * Create a table row for the comparison counts.
     *
     * @param label		description of the count
     * @param newValue	count for this genome
     * @param oldValue	count for reference genome
     */
    private DomContent compareTableRow(String label, int newValue, int oldValue) {
        double percent = (oldValue > 0 ? newValue * 100.0 / oldValue : 0);
        return tr(td(label), Html.numCell(newValue), Html.numCell(oldValue), Html.numCell(percent));
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
                    if (match.getDistance() <= maxProtDist) {
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
    protected DomContent advancedFeatureComment(Feature feat, GenomeEval gReport, String role) {
        DomContent retVal = null;
        // Only proceed if we found a reference genome.
        if (this.refGenomeId != null) {
            // Do we have a close reference feature?
            String closeFid = this.closeFeatureMap.get(feat.getId());
            if (closeFid != null) {
                // We found one. Form our message.
                retVal = join(this.refGenomeObj.featureLink(closeFid), "in the reference genome is close, and performs the same role.");
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
    protected void advancedRoleComment(ContainerTag list, GenomeEval gReport, String role) {
        // Only proceed if a reference genome was found.
        if (this.refGenomeId != null) {
            // Get the features for this role from the kmer map.
            KmerCollectionGroup roleFeatures = this.featureFinder.get(role);
            if (roleFeatures == null || roleFeatures.size() == 0) {
                list.with(li("No features perform this role in the reference genome."));
            } else {
                Collection<String> fidList = roleFeatures.getKeys();
                list.with(li(join("This role is performed by ", this.refGenomeObj.featureListLink(fidList), "in the reference genome.")));
            }
        }
    }

    /**
     * Compute the reference genome and fill in the protein kmer database.
     *
     * @param gReport	quality report on the genome of interest
     */
    protected void advancedGenomeSetup(GenomeStats gReport) {
        // Compute the reference genome.
        this.refGenomeObj = this.refEngine.ref(gReport.getGenome());
        if (this.refGenomeObj == null) {
            this.refGenomeId = null;
        } else {
            this.refGenomeId = this.refGenomeObj.getId();
            log.info("Analyzing reference genome {}: {}.", this.refGenomeId, refGenomeName());
            // Create the problematic role directory.
            Set<String> targetRoles = gReport.getProblematicRoles();
            this.featureFinder = new HashMap<String, KmerCollectionGroup>(targetRoles.size());
            for (String role : targetRoles) {
                featureFinder.put(role, new KmerCollectionGroup());
            }
            // Get the role definitions.
            RoleMap roleMap = this.getRoleMap();
            // Loop through the features, putting them in the directory.
            for (Feature feat : this.refGenomeObj.getPegs()) {
                Collection<Role> roles = feat.getUsefulRoles(roleMap);
                for (Role role : roles) {
                    KmerCollectionGroup roleGroup = this.featureFinder.get(role.getId());
                    if (roleGroup != null)
                        roleGroup.addSequence(feat.getProteinTranslation(), feat.getId());
                }
            }
        }
        // Finally, initialize the close-feature map.
        this.closeFeatureMap = new HashMap<String, String>();
    }

    /**
     * Add the ORF report to the summary page.
     */
    @Override
    protected DomContent advancedSummaryReport() {
        DomContent retVal = null;
        if (this.orfReportRows.size() > 1)
            retVal = Html.formatTable("ORF Comparison Summary", this.orfReportRows);
        return retVal;
    }


    /**
     * @return the name of the current reference genome
     */
    private String refGenomeName() {
        return this.refGenomeObj.getName();
    }

    /**
     * Set the maximum protein distance for the protein comparisons.
     *
     * @param newDistance	new maximum distance
     */
    public void setSensitivity(double newDistance) {
        this.maxProtDist = newDistance;
    }

    /**
     * Specify the engine for computing reference genomes.
     *
     * @param refEngine		object for computing reference genome IDs
     */
    public void setEngine(RefGenomeComputer refEngine) {
        this.refEngine = refEngine;
    }

    /**
     * Initialize the reference genome engine for this batch
     *
     * @param reports	array of evaluated genomes
     */
    @Override
    public void setupGenomes(GenomeStats[] reports) {
        this.refEngine.setupReferences(reports);
    }
}
