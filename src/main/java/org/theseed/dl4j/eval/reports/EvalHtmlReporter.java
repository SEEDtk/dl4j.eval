/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.commons.io.FileUtils;
import org.theseed.dl4j.eval.ContigAnalysis;
import org.theseed.dl4j.eval.GenomeStats;
import org.theseed.dl4j.eval.GenomeStats.FeatureStatus;
import org.theseed.dl4j.eval.GenomeStats.ProblematicRole;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.compare.CompareFeatures;
import org.theseed.locations.Location;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.reports.Html;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * This class produces HTML versions of the reports.  All of them have ".html" suffixes.
 *
 * @author Bruce Parrello
 *
 */
public class EvalHtmlReporter extends EvalReporter {

    /**
      * Subclass for sorting genomes by score and then ID.
      */
    protected static class EvalSorter implements Comparable<EvalSorter> {

        public static final EvalSorter HEADER = new EvalSorter();

        private double score;
        private String id;

        public EvalSorter(GenomeStats gReport) {
            this.score = gReport.getScore();
            this.id = gReport.getId();
        }

        public EvalSorter() {
            this.score = Double.POSITIVE_INFINITY;
            this.id = "";
        }

        @Override
        public int compareTo(EvalSorter other) {
            // Sort descending by score, then ascending by name.
            int retVal = Double.compare(other.score, this.score);
            if (retVal == 0)
                retVal = this.id.compareTo(other.id);
            return retVal;
        }

    }

     // FIELDS

    /** list of good genome table rows */
    SortedMap<EvalSorter, DomContent> goodRows;
    /** list of bad genome table rows */
    SortedMap<EvalSorter, DomContent> badRows;
    /** map of contig comments */
    Map<String, ContigAnalysis> contigMap;

    /** length of a short feature */
    private static final int SHORT_FEATURE = 180;


    @Override
    protected void initialize(File modelDir) throws IOException {
    }

    @Override
    protected void startSummary() throws IOException {
        // Create the two accumulating tables.
        this.goodRows = new TreeMap<EvalSorter, DomContent>();
        this.badRows = new TreeMap<EvalSorter, DomContent>();
        // Put a header on each.
        ContainerTag headerRow = tr(
                th("Score").withClass("num"),
                th("Genome ID"),
                th("Genome Name"),
                th("Coarse %").withClass("num"),
                th("Fine %").withClass("num"));
        if (this.hasCompleteness())
            headerRow.with(th("Complete %").withClass("num"))
                    .with(th("Contamination %").withClass("num"))
                    .with(th("Completeness Group"));
        headerRow.with(th("Contigs").withClass("num"))
                .with(th("Base Pairs").withClass("num"))
                .with(th("Hypothetical %").withClass("num"))
                .with(th("Good Phes").withClass("flag"))
                .with(th("SSU rRNA").withClass("flag"))
                .with(th("Quality").withClass("flag"));
        this.goodRows.put(EvalSorter.HEADER, headerRow);
        this.badRows.put(EvalSorter.HEADER, headerRow);
    }

    @Override
    protected void writeDetails(GenomeStats gReport) throws IOException {
        // Get the genome ID and compute the output file name.
        String genomeId = gReport.getId();
        File outFile = this.htmlFile(genomeId);
        // Determine if we're good or poor and get the genome name.
        String rating = (gReport.isGood() ? "good" : "poor");
        String seedRating = (gReport.isGoodSeed() ? "good" : "missing or invalid");
        String ssuRating = (gReport.hasSsuRRna() ? "known" : "unknown");
        String gName = gReport.getName();
        // Perform any special genome initialization.
        this.advancedGenomeSetup(gReport);
        // Analyze the contigs.
        this.analyzeContigs(gReport);
        // Analyze the features.
        this.analyzeFeatures(gReport);
        // Create the PPR rows.
        ArrayList<DomContent> roleRows = new ArrayList<DomContent>();
        DomContent headerRow = tr(th("Role description"), th("Predicted").withClass("num"),
                th("Actual").withClass("num"), th("Comment"));
        roleRows.add(headerRow);
        // We need to count the roles we process.
        int overCount = 0;
        int underCount = 0;
        // Loop through the problematic roles.
        Set<String> badRoles = gReport.getProblematicRoles();
        for (String role : badRoles) {
            String roleName = this.getRoleName(role);
            GenomeStats.ProblematicRole ppr = gReport.getReport(role);
            int predicted = ppr.getPredicted();
            int actual = ppr.getActual();
            DomContent comment = computeComment(gReport, ppr, role);
            if (actual < predicted) {
                underCount++;
            } else {
                overCount++;
            }
            roleRows.add(tr(td(roleName), Html.numCell(predicted), Html.numCell(actual), td(comment)));
        }
        // Create the labeling detail rows.
        ArrayList<DomContent> detailRows = new ArrayList<DomContent>();
        Html.detailRow(detailRows, "Genome ID", td(gReport.getGenome().genomeLink()));
        Html.detailRow(detailRows, "Genome Name", td(gName));
        // Ask the subclass for any additional rows.
        advancedDetailRows(detailRows);
        // Fill in all the quality-data statistic rows.
        qualityRows(gReport, detailRows, this.hasCompleteness());
        Html.detailRow(detailRows, "Overpresent Roles", Html.numCell(overCount));
        Html.detailRow(detailRows, "Underpresent Roles", Html.numCell(underCount));
        // Add the coverage if this came from a bin.
        double coverage = gReport.getGenome().getBinCoverage();
        if (coverage > 0.0) {
            Html.detailRow(detailRows, "Mean Contig Coverage", Html.numCell(coverage));
        }
        // Ask for extra tables.
        DomContent extraHtml = this.advancedGenomeAnalysis(gReport);
        // Format the page.
        String page = Html.page(gName + " Evaluation Report",
                    h1("Evaluation Report for " + genomeId),
                    p(String.format("This genome has an overall score of %4.2f using evaluator version %s and is of %s quality." +
                            "The PheS protein is %s.  The SSU rRNA sequence is %s.",
                            gReport.getScore(), this.getVersion(), rating, seedRating, ssuRating)),
                    div(table().with(detailRows.stream()).withClass(Html.TABLE_CLASS)).withClass("shrinker"),
                    extraHtml,
                    Html.formatTable("Potentially Problematic Roles", roleRows)
                );
        FileUtils.writeStringToFile(outFile, page, "UTF-8");
    }

    /**
     * Add the standard quality statistic rows to a table.
     *
     * @param gReport		evaluated genome
     * @param qualityRows	table being built
     * @param completeness	TRUE if completeness data is valid, else FALSE
     */
    public static void qualityRows(GenomeStats gReport, List<DomContent> qualityRows, boolean completeness) {
        Html.detailRow(qualityRows, "Completeness Group", td(gReport.getGroup()));
        Html.detailRow(qualityRows, "Domain", td(gReport.getDomain()));
        Html.detailRow(qualityRows, "DNA size (base pairs)", Html.numCell(gReport.getDnaSize()));
        Html.detailRow(qualityRows, "Number of Contigs", Html.numCell(gReport.getContigCount()));
        Html.detailRow(qualityRows, "Contig L50", Html.numCell(gReport.getL50()));
        Html.detailRow(qualityRows, "Contig N50", Html.numCell(gReport.getN50()));
        Html.detailRow(qualityRows, "Consistency Roles", Html.numCell(gReport.getConsistencyRoleCount()));
        if (completeness) {
            Html.detailRow(qualityRows, "Completeness Roles", Html.numCell(gReport.getCompletenessRoleCount()));
            Html.detailRow(qualityRows, "Shared Completeness / Consistency Roles", Html.numCell(gReport.getCommonRoleCount()));
        }
        Html.detailRow(qualityRows, "CDS Features", Html.numCell(gReport.getPegCount()));
        Html.detailRow(qualityRows, "CDS Features in Local Protein Families", Html.numCell(gReport.getPlfamCount()));
        Html.detailRow(qualityRows, "CDS Features without annotation", Html.numCell(gReport.getHypoCount()));
        Html.detailRow(qualityRows, "CDS Features with annotation", Html.numCell(gReport.getPegCount() - gReport.getHypoCount()));
        Html.detailRow(qualityRows, "Coarse Consistency %", Html.numCell(gReport.getCoarsePercent()));
        Html.detailRow(qualityRows, "Fine Consistency %", Html.colorCell(gReport.isConsistent(), gReport.getFinePercent()));
        if (completeness) {
            Html.detailRow(qualityRows, "Completeness %", Html.colorCell(gReport.isComplete(), gReport.getCompletePercent()));
            Html.detailRow(qualityRows, "Contamination %", Html.colorCell(gReport.isClean(), gReport.getContaminationPercent()));
        }
        Html.detailRow(qualityRows, "CDS Coverage %", Html.numCell(gReport.getCdsPercent()));
        Html.detailRow(qualityRows, "Hypothetical Protein %", Html.colorCell(gReport.isUnderstood(), gReport.getHypotheticalPercent()));
        Html.detailRow(qualityRows, "% CDS Features in Local Protein Families", Html.numCell(gReport.getPlfamPercent()));
    }

    /**
     * Subclasses can use this to create additional data in the output.
     *
     * @param gReport	quality
     *
     * @return the HTML to insert in the output report
     */
    protected DomContent advancedGenomeAnalysis(GenomeStats gReport) {
        return null;
    }

    /**
     * Subclasses can use this to add additional detail rows to the genome page.
     *
     * @param detailRows	list of accumulating table rows
     */
    protected void advancedDetailRows(ArrayList<DomContent> detailRows) {
    }

    /**
     * Subclasses can initialize special genome-related data structures here.
     *
     * @param gReport	quality report on the genome of interest
     */
    protected void advancedGenomeSetup(GenomeStats gReport) {
    }

    /**
     * Determine which features are good and bad.
     *
     * @param gReport	quality report for this genome
     */
    private void analyzeFeatures(GenomeStats gReport) {
        // Get the role map.
        RoleMap roleDefinitions = this.getRoleMap();
        // Loop through all the protein features in this genome.
        for (Feature feat : gReport.getGenome().getPegs()) {
            // Determine the status of the feature.  This also adds it to the problematic role list if it is bad.
            Collection<Role> roles = feat.getUsefulRoles(roleDefinitions);
            GenomeStats.FeatureStatus status = gReport.checkProblematicRoles(feat, roles);
            if (status == GenomeStats.FeatureStatus.BAD)
                status = this.advancedFeatureAnalysis(feat, roles);
            // Get the feature's contig, and count the feature on it.
            String contigId = feat.getLocation().getContigId();
            ContigAnalysis contigData = this.contigMap.get(contigId);
            contigData.countFeature(feat, status);
        }
    }

    /**
     * @return a possibly-modified status for a bad feature
     *
     * @param feat		feature of interest
     * @param roles		roles performed by the feature
     */
    protected FeatureStatus advancedFeatureAnalysis(Feature feat, Collection<Role> roles) {
        // The default is not to modify the status.
        return FeatureStatus.BAD;
    }

    /**
     * Run through all the contigs, creating the contig analysis map.
     *
     * @param gReport	quality report for this genome
     */
    private void analyzeContigs(GenomeStats gReport) {
        Genome genome = gReport.getGenome();
        this.contigMap = new HashMap<String, ContigAnalysis>(genome.getContigCount());
        for (Contig contig : genome.getContigs()) {
            ContigAnalysis analysis = new ContigAnalysis(contig, gReport);
            this.contigMap.put(contig.getId(), analysis);
        }
    }

    /**
     * Create the comments for this problematic role.  If the role occurs in a feature, the feature must be linked
     * in one of the comments.
     *
     * @param gReport	the evaluation report for the current genome
     * @param ppr		the problematic role object, with the features filled in
     * @param role		the ID of the role
     *
     * @return HTML describing the role.
     */
    private DomContent computeComment(GenomeStats gReport, ProblematicRole ppr, String role) {
        ContainerTag retVal = ul();
        if (ppr.getActual() == 0) {
            if (ppr.isUniversal()) {
                retVal.with(li("Missing universal role."));
            }
            // Allow the subclass to say more.
            this.advancedRoleComment(retVal, gReport, role);
        } else {
            // Here we need to comment about each feature containing the role.  These are put in a bullet list.
            if (ppr.isUniversal())
                retVal.with(li("Redundant universal role, indicating possible contamination."));
            for (Feature feat : ppr.getFeatures()) {
                // Figure out where the feature is.
                Location loc = feat.getLocation();
                ContigAnalysis contigObject = this.contigMap.get(loc.getContigId());
                // Get the feature's contig-related comment.
                DomContent contigComment = contigObject.locationComment(loc);
                // Form the full feature comment.
                DomContent featureComment = li(join(gReport.getGenome().featureRegionLink(feat.getId()),
                        iff(loc.getLength() < SHORT_FEATURE, text("is short and")), contigComment,
                        this.advancedFeatureComment(feat, gReport, role)));
                retVal.with(featureComment);
            }
            // Allow the subclass to say more.
            this.advancedRoleComment(retVal, gReport, role);
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
    protected DomContent advancedFeatureComment(Feature feat, GenomeStats gReport, String role) {
        return null;
    }

    /**
     * This is a stub that a subclass can use to create more advanced comments about a missing role.
     *
     * @param listRows	an HTML list to which the advanced comment can be added
     * @param gReport	genome quality report
     * @param role		role of interest
     */

    protected void advancedRoleComment(ContainerTag list, GenomeStats gReport, String role) {
    }

    /**
     * @return a table containing statistics about an ORF comparison
     *
     * @param comparison	ORF comparison results
     * @return
     */
    public DomContent compareReport(CompareFeatures comparison) {
        List<DomContent> tableRows = new ArrayList<DomContent>(10);
        Html.detailRow(tableRows, "Number of ORFs only annotated in the reference genome", Html.numCell(comparison.getOldOnlyCount()));
        Html.detailRow(tableRows, "Number of ORFs only annotated in this genome", Html.numCell(comparison.getNewOnlyCount()));
        Html.detailRow(tableRows, "Number of ORFs annotated in both genomes", Html.numCell(comparison.getCommon()));
        Html.detailRow(tableRows, "Number of ORFs annotated in both genomes with identical functions and lengths", Html.numCell(comparison.getIdentical()));
        Html.detailRow(tableRows, "Number of ORFs annotated in both genomes, but with different functions", Html.numCell(comparison.getDifferentFunctions()));
        Html.detailRow(tableRows, "Number of identically-annotated ORFs whose proteins are longer in the reference genome", Html.numCell(comparison.getShorter()));
        Html.detailRow(tableRows, "Number of identically-annotated ORFs whose proteins are longer in this genome", Html.numCell(comparison.getLonger()));
        DomContent retVal = join(h2("Comparison of Annotation Results by Open Reading Frame"),
                div(table().with(tableRows.stream()).withClass(Html.TABLE_CLASS)).withClass("shrinker"));
        return retVal;
    }


    @Override
    protected void writeSummary(GenomeStats gReport) throws IOException {
        ContainerTag detailRow = tr(
                td(Html.gPageLink(gReport.getId(), Html.num(gReport.getScore()))).withClass("num"),
                td(gReport.getGenome().genomeLink()),
                td(gReport.getName()),
                Html.numCell(gReport.getCoarsePercent()),
                Html.colorCell(gReport.isConsistent(), gReport.getFinePercent()));
        if (this.hasCompleteness()) {
            detailRow.with(Html.colorCell(gReport.isComplete(), gReport.getCompletePercent()))
                .with(Html.colorCell(gReport.isClean(), gReport.getContaminationPercent()))
                .with(td(gReport.getGroup()));
        }
        detailRow.with(Html.numCell(gReport.getContigCount()))
            .with(Html.numCell(gReport.getDnaSize()))
            .with(Html.colorCell(gReport.isUnderstood(), gReport.getHypotheticalPercent()))
            .with(Html.flagCell(gReport.isGoodSeed(), "Y", ""))
            .with(Html.flagCell(gReport.hasSsuRRna(), "Y", ""))
            .with(Html.flagCell(gReport.isGood(), "Good", "Poor"));
        if (gReport.isGood())
            this.goodRows.put(new EvalSorter(gReport), detailRow);
        else
            this.badRows.put(new EvalSorter(gReport), detailRow);
    }

    @Override
    protected void endSummary() throws IOException {
        // Compute the number of good and bad genomes.
        int good = this.getGoodCount();
        int total = this.getGenomeCount();
        int bad = total - good;
        DomContent countNotes = formatCounts();
        // Create genome tables.
        DomContent goodRegion = (good == 0 ? p() : Html.formatTable("Good Genomes", goodRows.values()));
        DomContent badRegion = (bad == 0 ? p() : Html.formatTable("Poor Genomes", badRows.values()));
        DomContent extraRegion = this.advancedSummaryReport();
        // Here we render the whole page.
        String page = Html.page("Evaluation Summary Report",
                    h1("Evaluation Summary Report"),
                    p(String.format("%d genomes processed using evaluator version %s. %d good and %d poor.",
                            this.getGenomeCount(), this.getVersion(), good, bad)),
                    countNotes,
                    goodRegion,
                    badRegion,
                    extraRegion
                );
        File summaryFile = new File(this.getOutDir(), "index.html");
        FileUtils.writeStringToFile(summaryFile, page, "UTF-8");
    }

    /**
     * Subclasses can use this to add extra content to the summary page.
     */
    protected DomContent advancedSummaryReport() {
        return null;
    }

    /**
     * @return a paragraph about the various quality counts
     */
    private DomContent formatCounts() {
        ArrayList<DomContent> notes = new ArrayList<DomContent>();
        formatCount(notes, this.getConsistentCount(), "inconsistent");
        if (this.hasCompleteness()) {
            formatCount(notes, this.getCompleteCount(), "incomplete");
            formatCount(notes, this.getCleanCount(), "contaminated");
        }
        formatCount(notes, this.getUnderstoodCount(), "insufficiently-characterized");
        int bad = this.getGenomeCount() - this.getGoodSeedCount();
        if (bad > 0) {
            notes.add(li(String.format("%d genomes have a missing or improper PheS protein", bad)));
        }
        bad = this.getGenomeCount() - this.getSsuFoundCount();
        if (bad > 0) {
            notes.add(li(String.format("%d genomes had an unknown SSU rRNA sequence", bad)));
        }
        DomContent retVal;
        if (notes.size() == 0) {
            retVal = p("All genomes are good.");
        } else {
            retVal = ul().with(notes.stream());
        }
        return retVal;
    }

    /**
     * Store a list item about this type of problem if any occurred.
     *
     * @param notes		buffer for list items
     * @param count		count of genomes without the problem
     * @param name		name of the problem
     */
    protected void formatCount(ArrayList<DomContent> notes, int count, String name) {
        int bad = this.getGenomeCount() - count;
        if (bad > 0) {
            notes.add(li(String.format("%d %s genomes.", bad, name)));
        }
    }

    @Override
    protected void finish() {
        // Insure we are not holding onto the table rows.
        this.goodRows = null;
        this.badRows = null;
    }

    @Override
    public void setupGenomes(GenomeStats[] reports) {
    }

}
