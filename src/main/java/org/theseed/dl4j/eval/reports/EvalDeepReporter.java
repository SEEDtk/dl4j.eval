/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;

import org.theseed.dl4j.eval.stats.GenomeAnalysis;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.genome.compare.CompareFeatures;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.RoleMap;
import org.theseed.reports.Html;

import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * This is an advanced version of the HTML reporting service that includes comparisons to a reference genome.
 *
 * @author Bruce Parrello
 *
 */
public class EvalDeepReporter extends EvalHtmlReporter implements IRefReporter {

    /** reference genome object */
    private Genome refGenomeObj;
    /** genome comparator */
    private CompareFeatures compareObj;
    /** orf report */
    private List<DomContent> orfReportRows;

    /**
     * Construct a deep HTML reporting object.
     */
    public EvalDeepReporter() {
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
        if (this.refGenomeObj != null) {
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
        if (this.refGenomeObj != null) {
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
            tableRows.add(compareTableRow("Total DNA length of the genome", gReport.getDnaSize(), this.refGenomeObj.getLength()));
            retVal = Html.formatTable("Comparison of " + gReport.getId() + " with Reference Genome " + this.refGenomeObj.getId(),
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
     * Compute the reference genome.
     *
     * @param gReport	quality report on the genome of interest
     * @param analysis 	analysis of genome of interest
     */
    @Override
    protected void advancedGenomeSetup(GenomeStats gReport, GenomeAnalysis analysis) {
        this.refGenomeObj = analysis.getRefGenome();
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

}
