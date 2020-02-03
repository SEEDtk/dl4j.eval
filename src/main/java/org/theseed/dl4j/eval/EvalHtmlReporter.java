/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.theseed.dl4j.eval.GenomeStats.FeatureStatus;
import org.theseed.dl4j.eval.GenomeStats.ProblematicRole;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

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

    private static final String CSS_HREF = "https://patricbrc.org/js/3.6.2/p3/resources/p3.css";

    private static final String NUM_FORMAT = "%-8.2f";

    private static final String INT_FORMAT = "%d";

    protected static final String TABLE_CLASS = "p3basic";

    private static final String BODY_CLASS = "claro";

    private static final String EXTRA_STYLES = "	table." + TABLE_CLASS + " th {\n" +
            "		background-color: #c8c8c8;\n" +
            "	}\n" +
            "	table." + TABLE_CLASS + " th.num, table." + TABLE_CLASS + " td.num {\n" +
            "		text-align: right;\n" +
            "	}\n" +
            "	table." + TABLE_CLASS + " th.flag, table." + TABLE_CLASS + " td.flag {\n" +
            "		text-align: center;\n" +
            "	}\n" +
            "   h1, h2 {\n" +
            "		font-weight: bolder;\n" +
            "	}\n" +
            "   h1, h2, p, ul {\n" +
            "		margin: 12px 12px 0px 12px;\n" +
            "	}\n" +
            "   table.p3basic ul {\n" +
            "		margin: 0px;\n" +
            "		list-style: disc outside none;\n" +
            "		padding-left: 20px;\n" +
            "	}\n" +
            "   table.p3basic li {\n" +
            "		margin: 3px 0px;\n" +
            "	}\n" +
            "   div.wrapper {\n" +
            "       margin: 12px;\n" +
            "   }\n" +
            "   div.shrinker {\n" +
            "       margin: 12px;\n" +
            "       display: inline-block;\n" +
            "       min-width: 0;\n" +
            "       width: auto;\n" +
            "   }\n" +
            "   li {\n" +
            "		margin: 6px 12px 0px 12px;\n" +
            "	}\n" +
            "	table." + TABLE_CLASS + " {\n" +
            "	    display:table;\n" +
            "	}\n";

    /** background color for values indicating bad genomes */
    private static final String BAD_STYLE = "background-color: gold;";

    /** length of a short feature */
    private static final int SHORT_FEATURE = 180;


    // LINK CONSTANTS

    private static final String GENOME_LINK_FMT = "https://www.patricbrc.org/view/Genome/%s";

    private static final String PAGE_LINK_FMT = "%s.html";

    private static final String FEATURE_CR_LINK = "https://www.patricbrc.org/view/Feature/%s#view_tab=compareRegionViewer";

    private static final String FEATURE_VIEW_LINK = "https://www.patricbrc.org/view/Feature/%s";

    private static final String FEATURE_LIST_LINK = "https://www.patricbrc.org/view/FeatureList/?in(patric_id,(\"%s\"))";

     // FIELDS

    /** list of good genome table rows */
    ArrayList<DomContent> goodRows;
    /** list of bad genome table rows */
    ArrayList<DomContent> badRows;
    /** map of contig comments */
    Map<String, ContigAnalysis> contigMap;

    /**
     * @param outDir	output directory
     */
    public EvalHtmlReporter(File outDir) {
        super(outDir);
    }

    @Override
    protected void initialize(File modelDir) throws IOException {
    }

    @Override
    protected void startSummary() throws IOException {
        // Create the two accumulating tables.
        this.goodRows = new ArrayList<DomContent>();
        this.badRows = new ArrayList<DomContent>();
        // Put a header on each.
        DomContent headerRow = tr(
                th("Score").withClass("num"),
                th("Genome ID"),
                th("Genome Name"),
                th("Coarse %").withClass("num"),
                th("Fine %").withClass("num"),
                th("Complete %").withClass("num"),
                th("Contamination %").withClass("num"),
                th("Completeness Group"),
                th("Contigs").withClass("num"),
                th("Base Pairs").withClass("num"),
                th("Hypothetical %").withClass("num"),
                th("Good Phes").withClass("flag"),
                th("Quality").withClass("flag")
            );
        this.goodRows.add(headerRow);
        this.badRows.add(headerRow);
    }

    /**
     * Display the specified value in an alternate color if the flag is FALSE.
     *
     * @param flag	TRUE if the value is good, else FALSE
     * @param str	text to display
     */
    protected ContainerTag colorCell(boolean flag, String str) {
        DomContent cell;
        if (str.isEmpty())
            cell = rawHtml("&nbsp;");
        else
            cell = text(str);
        ContainerTag retVal = td(cell);
        if (! flag) {
            retVal = retVal.withStyle(BAD_STYLE);
        }
        return retVal;
    }

    /**
     * @return a link to the compare-regions viewer for a feature
     *
     * @param fid id of the target feature
     */
    protected DomContent featureRegionLink(String fid) {
        return a(fid).withHref(String.format(FEATURE_CR_LINK, fid))
                .withTarget("_blank");
    }

    /**
     * @return a link to the specified genome along with the genome ID
     *
     * @param genome_id		ID of the target genome
     */
    protected DomContent genomeLink(String genome_id) {
        return a(genome_id).withHref(String.format(GENOME_LINK_FMT, genome_id)).withTarget("_blank");
    }

    /**
     * @return a link to the specified genome's report page
     *
     * @param genome_id		ID of the target genome
     * @param text			text for the link
     */
    protected DomContent gPageLink(String genome_id, String text) {
        return a(text).withHref(String.format(PAGE_LINK_FMT, genome_id));
    }

    /**
     * Display the specified value in an alternate color if the flag is FALSE.
     *
     * @param flag	TRUE if the value is good, else FALSE
     * @param val	floating-point value to display
     */
    protected ContainerTag colorCell(boolean flag, double val) {
        return colorCell(flag, num(val)).withClass("num");
    }

    /**
     * Display the specified floating-point value in a table cell.
     *
     * @param val	floating-point value to display
     */
    protected ContainerTag numCell(double val) {
        return td(num(val)).withClass("num");
    }

    /**
     * Display the specified integer value in a table cell.
     *
     * @param val	floating-point value to display
     */
    protected ContainerTag numCell(int val) {
        return td(num(val)).withClass("num");
    }

    /**
     * Display a flag value in a cell.  A false value is in the alternate color.
     *
     * @param flag		TRUE if the relevant attribute is acceptable, else FALSE
     * @param trueText	text to display if the flag is TRUE
     * @param falseText	test to display if the flag is FALSE
     */
    protected ContainerTag flagCell(boolean flag, String trueText, String falseText) {
        String text = (flag ? trueText : falseText);
        ContainerTag retVal = this.colorCell(flag, text);
        return retVal.withClass("flag");
    }

    @Override
    protected void writeDetails(GenomeStats gReport) throws IOException {
        // Get the genome ID and compute the output file name.
        String genomeId = gReport.getId();
        File outFile = new File(this.getOutDir(), genomeId + ".html");
        // Determine if we're good or poor and get the genome name.
        String rating = (gReport.isGood() ? "good" : "poor");
        String seedRating = (gReport.isGoodSeed() ? "good" : "missing or invalid");
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
            roleRows.add(tr(td(roleName), numCell(predicted), numCell(actual), td(comment)));
        }
        // Create the labeling detail rows.
        ArrayList<DomContent> detailRows = new ArrayList<DomContent>();
        detailRow(detailRows, "Genome ID", td(genomeLink(genomeId)));
        detailRow(detailRows, "Genome Name", td(gName));
        // Ask the subclass for any additional rows.
        advancedDetailRows(detailRows);
        // Fill in all the quality-data statistic rows.
        detailRow(detailRows, "Completeness Group", td(gReport.getGroup()));
        detailRow(detailRows, "Domain", td(gReport.getDomain()));
        detailRow(detailRows, "DNA size (base pairs)", numCell(gReport.getDnaSize()));
        detailRow(detailRows, "Number of Contigs", numCell(gReport.getContigCount()));
        detailRow(detailRows, "Contig L50", numCell(gReport.getL50()));
        detailRow(detailRows, "Contig N50", numCell(gReport.getN50()));
        detailRow(detailRows, "Consistency Roles", numCell(gReport.getConsistencyRoleCount()));
        detailRow(detailRows, "Completeness Roles", numCell(gReport.getCompletenessRoleCount()));
        detailRow(detailRows, "Shared Completeness / Consistency Roles", numCell(gReport.getCommonRoleCount()));
        detailRow(detailRows, "Overpresent Roles", numCell(overCount));
        detailRow(detailRows, "Underpresent Roles", numCell(underCount));
        detailRow(detailRows, "CDS Features", numCell(gReport.getPegCount()));
        detailRow(detailRows, "CDS Features in Local Protein Families", numCell(gReport.getPlfamCount()));
        detailRow(detailRows, "CDS Features without annotation", numCell(gReport.getHypoCount()));
        detailRow(detailRows, "CDS Features with annotation", numCell(gReport.getPegCount() - gReport.getHypoCount()));
        detailRow(detailRows, "Coarse Consistency %", numCell(gReport.getCoarsePercent()));
        detailRow(detailRows, "Fine Consistency %", colorCell(gReport.isConsistent(), gReport.getFinePercent()));
        detailRow(detailRows, "Completeness %", colorCell(gReport.isComplete(), gReport.getCompletePercent()));
        detailRow(detailRows, "Contamination %", colorCell(gReport.isClean(), gReport.getContaminationPercent()));
        detailRow(detailRows, "CDS Coverage %", numCell(gReport.getCdsPercent()));
        detailRow(detailRows, "Hypothetical Protein %", colorCell(gReport.isUnderstood(), gReport.getHypotheticalPercent()));
        detailRow(detailRows, "% CDS Features in Local Protein Families", numCell(gReport.getPlfamPercent()));
        // Ask for extra tables.
        DomContent extraHtml = this.advancedGenomeAnalysis(gReport);
        // Format the page.
        String page = html(
                head(
                    title(gName + " Evaluation Report"),
                    link().withRel("stylesheet").withHref(CSS_HREF).withType("text/css"),
                    style(EXTRA_STYLES).withType("text/css")
                ),
                body(
                    h1("Evaluation Report for " + genomeId),
                    p(String.format("This genome has an overall score of %4.2f using evaluator version %s and is of %s quality." +
                            "The PheS protein is %s.",
                            gReport.getScore(), this.getVersion(), rating, seedRating)),
                    div(table().with(detailRows.stream()).withClass(TABLE_CLASS)).withClass("shrinker"),
                    extraHtml,
                    formatTable("Potentially Problematic Roles", roleRows)
                ).withClass(BODY_CLASS)
            ).render();
        FileUtils.writeStringToFile(outFile, page, "UTF-8");
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
                DomContent featureComment = li(join(this.featureRegionLink(feat.getId()),
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
     * Add a row to the statistical details table row collection.
     *
     * @param detailRows	statistical details table row collection
     * @param label		label for the row
     * @param cell			table cell with the data
     */
    protected void detailRow(List<DomContent> detailRows, String label, ContainerTag cell) {
        detailRows.add(tr(th(label), cell));
    }

    @Override
    protected void writeSummary(GenomeStats gReport) throws IOException {
        DomContent detailRow = tr(
                td(gPageLink(gReport.getId(), num(gReport.getScore()))).withClass("num"),
                td(genomeLink(gReport.getId())),
                td(gReport.getName()),
                numCell(gReport.getCoarsePercent()),
                colorCell(gReport.isConsistent(), gReport.getFinePercent()),
                colorCell(gReport.isComplete(), gReport.getCompletePercent()),
                colorCell(gReport.isClean(), gReport.getContaminationPercent()),
                td(gReport.getGroup()),
                numCell(gReport.getContigCount()),
                numCell(gReport.getDnaSize()),
                colorCell(gReport.isUnderstood(), gReport.getHypotheticalPercent()),
                flagCell(gReport.isGoodSeed(), "Y", ""),
                flagCell(gReport.isGood(), "Good", "Poor")
            );
        if (gReport.isGood())
            this.goodRows.add(detailRow);
        else
            this.badRows.add(detailRow);
    }

    @Override
    protected void endSummary() throws IOException {
        // Compute the number of good and bad genomes.
        int good = this.getGoodCount();
        int total = this.getGenomeCount();
        int bad = total - good;
        DomContent countNotes = formatCounts();
        // Create genome tables.
        DomContent goodRegion = (good == 0 ? p() : formatTable("Good Genomes", goodRows));
        DomContent badRegion = (bad == 0 ? p() : formatTable("Poor Genomes", badRows));
        // Here we render the whole page.
        String page = html(
                head(
                    title("Evaluation Summary Report"),
                    link().withRel("stylesheet").withHref(CSS_HREF).withType("text/css"),
                    style(EXTRA_STYLES).withType("text/css")
                ),
                body(
                    div(
                        h1("Evaluation Summary Report"),
                        p(String.format("%d genomes processed using evaluator version %s. %d good and %d poor.",
                                this.getGenomeCount(), this.getVersion(), good, bad)),
                        countNotes,
                        goodRegion,
                        badRegion
                    ).withClass(BODY_CLASS)
                )
            ).render();
        File summaryFile = new File(this.getOutDir(), "index.html");
        FileUtils.writeStringToFile(summaryFile, page, "UTF-8");
    }

    /**
     * @return a formatted integer
     *
     * @param val	integer to format
     */
    protected String num(int val) {
        return String.format(INT_FORMAT, val);
    }

    /**
     * @return a formatted floating-point number
     *
     * @param val	floating-point number to format
     */
    protected String num(double val) {
        return String.format(NUM_FORMAT, val);
    }

    /**
     * @return a table created from the specified detail rows
     *
     * @param tableRows		rows to put in the table
     */
    protected DomContent formatTable(String header, ArrayList<DomContent> tableRows) {
        return join(
                h2(header),
                div(table().with(tableRows.stream()).withClass(TABLE_CLASS + " striped")).withClass("wrapper")
            );
    }

    /**
     * @return a paragraph about the various quality counts
     */
    private DomContent formatCounts() {
        ArrayList<DomContent> notes = new ArrayList<DomContent>();
        formatCount(notes, this.getConsistentCount(), "inconsistent");
        formatCount(notes, this.getCompleteCount(), "incomplete");
        formatCount(notes, this.getCleanCount(), "contaminated");
        formatCount(notes, this.getUnderstoodCount(), "insufficiently-characterized");
        int bad = this.getGenomeCount() - this.getGoodSeedCount();
        if (bad > 0) {
            notes.add(li(String.format("%d genomes have a missing or improper PheS protein", bad)));
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

    /**
     * Create a link to list all the features in the collection.
     *
     * @param fidList	a collection of feature IDs
     *
     * @return a hyperlink that describes the features and can access them
     */
    protected DomContent featureListLink(Collection<String> fidList) {
        DomContent retVal = null;
        if (fidList.size() == 1) {
            // Only one feature.  We go to the feature landing page and display the feature ID.
            String fid = fidList.iterator().next();
            retVal = featureLink(fid);
        } else {
            // Multiple features.  We go to a feature list view.  This requires the feature IDs to be enclosed in quotes.
            String rawUrl = String.format(FEATURE_LIST_LINK, StringUtils.join(fidList, "\",\""));
            // We also have to URLEncode the vertical bars.
            String linkUrl = StringUtils.replace(rawUrl, "|", "%7c");
            // Apply the URL to the text.
            String linkText = String.format("%d features", fidList.size());
            retVal = a(linkText).withHref(linkUrl).withTarget("_blank");
        }
        return retVal;
    }

    /** Create a link to view a single feature.
     *
     * @param fid	ID of the feature to view
     *
     * @return a hyperlink to the feature's view page
     */
    protected DomContent featureLink(String fid) {
        return a(fid).withHref(String.format(FEATURE_VIEW_LINK, fid)).withTarget("_blank");
    }

}
