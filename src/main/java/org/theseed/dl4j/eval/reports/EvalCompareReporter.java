/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static j2html.TagCreator.*;

import org.apache.commons.io.FileUtils;
import org.theseed.dl4j.eval.GenomeStats;
import org.theseed.genome.Annotation;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.compare.CompareFeatures;
import org.theseed.reports.Html;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

import j2html.tags.DomContent;

/**
 * This is an evaluation reporter that produces a detailed comparison of two genomes with the same contigs.
 * The target genome should have been annotated using kmers.anno.  The reference genome is presumed to
 * contain the ground truth.
 *
 * @author Bruce Parrello
 *
 */
public class EvalCompareReporter extends EvalReporter implements IRefReporter {

    // FIELDS

    /** reference-genome computation engine */
    private RefGenomeComputer refEngine;
    /** table rows for summary report */
    List<DomContent> summaryRows;
    /** genome comparator */
    private CompareFeatures compareObj;

    private static final Pattern EVIDENCE_PATTERN = Pattern.compile("Annotated with evidence (\\d+) and strength (\\d+\\.\\d+)");

    @Override
    public void setEngine(RefGenomeComputer refEngine) {
        this.refEngine = refEngine;
    }

    @Override
    protected void initialize(File modelDir) throws IOException {
        try {
            this.compareObj = new CompareFeatures();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    protected void startSummary() throws IOException {
        // Create the array for the summary lines.
        this.summaryRows = new ArrayList<DomContent>(50);
        // Add the header line.
        summaryRows.add(tr(th("Genome ID"), th("Genome Name"), th("Ref Genome"), th("Ref Genome Name"), th("Identical").withClass("num"), th("Same Function").withClass("num"),
                th("Different Function").withClass("num"), th("False Positive").withClass("num"), th("False Negative").withClass("num")));
    }

    @Override
    protected void writeDetails(GenomeStats gReport) throws IOException {
        // Get the genome ID and the output file.
        Genome newGenome = gReport.getGenome();
        String genomeId = newGenome.getId();
        File outFile = this.htmlFile(genomeId);
        // Get the reference genome.  Only proceed if we have one.
        Genome refGenome = this.refEngine.ref(newGenome);
        if (refGenome != null) {
            // Compare the genomes.  Only proceed if the compare works.
            boolean compared = this.compareObj.compare(newGenome, refGenome);
            if (compared) {
                // Create a FASTA output file for the false positives.
                try (FastaOutputStream fpOut = new FastaOutputStream(new File(this.getOutDir(), genomeId + ".faa"))) {
                    // Build the evaluation summary table.
                    List<DomContent> qualityRows = new ArrayList<DomContent>(25);
                    EvalHtmlReporter.qualityRows(gReport, qualityRows, this.hasCompleteness());
                    // Build the two role tables.  First, the false positives.
                    List<DomContent> tableRows = new ArrayList<DomContent>(this.compareObj.getNewOnlyCount());
                    tableRows.add(tr(th("Peg ID"), th("function"), th("Prot Len").withClass("num"), th("Evidence").withClass("num"),
                            th("Strength").withClass("num")));
                    // These track the minimum evidence and strength for good ORFs.
                    int minEvidence = Integer.MAX_VALUE;
                    double minStrength = 1.0;
                    // These track the number of good ORFs that fall within certain ranges.
                    int[] evidenceBy50Count = new int[10];
                    int[] strengthByP1Count = new int[10];
                    // Get the false positive set.
                    Set<Feature> newOnly = this.compareObj.getNewOnly();
                    // Get the features sorted by function.
                    SortedSet<Feature> sortedFeatures = new TreeSet<Feature>(new Feature.ByFunction());
                    sortedFeatures.addAll(newGenome.getPegs());
                    for (Feature peg : sortedFeatures) {
                        // Extract the strength and evidence from the annotations.
                        int evidence = Integer.MAX_VALUE;
                        double strength = 1.0;
                        for (Annotation anno : peg.getAnnotations()) {
                            Matcher match = EVIDENCE_PATTERN.matcher(anno.getComment());
                            if (match.matches()) {
                                evidence = Integer.parseInt(match.group(1));
                                strength = Double.parseDouble(match.group(2));
                            }
                        }
                        // If this is a false positive, put it in the FP table and the output FASTA. Otherwise, merge it into the trackers.
                        if (newOnly.contains(peg)) {
                            tableRows.add(tr(td(newGenome.featureRegionLink(peg.getId())), td(peg.getFunction()),
                                    Html.numCell(peg.getProteinLength()), Html.numCell(evidence),
                                    Html.numCell(strength * 3)));
                            Sequence fpSeq = new Sequence(peg.getId(), peg.getFunction(), peg.getProteinTranslation());
                            fpOut.write(fpSeq);
                        } else {
                            minEvidence = Math.min(evidence, minEvidence);
                            minStrength = Math.min(strength, minStrength);
                            int evidenceI = Math.min(evidence / 50, 9);
                            int strengthI = Math.min((int) (strength*30), 9);
                            evidenceBy50Count[evidenceI]++;
                            strengthByP1Count[strengthI]++;
                        }
                    }
                    DomContent fpTable = Html.formatTable("Proteins Only in This Genome", tableRows);
                    fpOut.flush();
                    // Add some more comparison data to the quality rows.
                    Html.detailRow(qualityRows, "Number of correctly-called functions",
                            Html.numCell(compareObj.getIdentical()+compareObj.getLonger()+compareObj.getShorter()));
                    Html.detailRow(qualityRows, "Number of differing functions", Html.numCell(compareObj.getDifferentFunctions()));
                    Html.detailRow(qualityRows, "Number of false positives", Html.numCell(compareObj.getNewOnlyCount()));
                    Html.detailRow(qualityRows, "Number of false negatives", Html.numCell(compareObj.getOldOnlyCount()));
                    Html.detailRow(qualityRows, "Minimum evidence for found ORFs", Html.numCell(minEvidence));
                    Html.detailRow(qualityRows, "Minimum strength for found ORFs", Html.numCell(minStrength * 3));
                    // Now, the false negatives, again sorted by function.
                    tableRows.clear();
                    tableRows.add(tr(th("Peg ID"), th("function"), th("Prot Len").withClass("num")));
                    sortedFeatures.clear();
                    sortedFeatures.addAll(this.compareObj.getOldOnly());
                    for (Feature peg : sortedFeatures) {
                        tableRows.add(tr(td(refGenome.featureRegionLink(peg.getId())), td(peg.getFunction()),
                                Html.numCell(peg.getProteinLength())));
                    }
                    DomContent fnTable = Html.formatTable("Proteins Not Found in This Genome", tableRows);
                    // Build the quality table.
                    DomContent qualityTable = div(table().with(qualityRows.stream()).withClass(Html.TABLE_CLASS)).withClass("shrinker");
                    // Build the tracker table.
                    List<DomContent> headerRow = new ArrayList<DomContent>(11);
                    List<DomContent> dataRow = new ArrayList<DomContent>(11);
                    headerRow.add(th(rawHtml("&nbsp;")));
                    dataRow.add(th("Evidence Counts"));
                    for (int i = 0; i < 9; i++) {
                        headerRow.add(th(String.format("%d to %d", i*50, i*50 + 49)).withClass("num"));
                        dataRow.add(Html.numCell(evidenceBy50Count[i]));
                    }
                    headerRow.add(th(">= 450"));
                    dataRow.add(Html.numCell(evidenceBy50Count[9]));
                    DomContent trackers = div(table().with(tr(each(headerRow.stream()))).with(tr(each(dataRow.stream())))
                            .withClass(Html.TABLE_CLASS)).withClass("wrapper");
                    headerRow.clear();
                    dataRow.clear();
                    headerRow.add(th(rawHtml("&nbsp;")));
                    dataRow.add(th("Strength Counts"));
                    for (int i = 0; i < 10; i++) {
                        headerRow.add(th(String.format("%3.1f to < %3.1f", i / 10.0, (i+1)/10.0)).withClass("num"));
                        dataRow.add(Html.numCell(strengthByP1Count[i]));
                    }
                    trackers = join(trackers, div(table().with(tr(each(headerRow.stream()))).with(tr(each(dataRow.stream())))
                            .withClass(Html.TABLE_CLASS)).withClass("wrapper"));

                    // Assemble all of this into a page.
                    String page = Html.page(newGenome.getName() + " Evaluation Report",
                                h1("Evaluation Report for " + genomeId + " " + newGenome.getName()),
                                qualityTable,
                                h2("Distribution of Found ORFs by Evidence and Strength"),
                                trackers,
                                fpTable,
                                fnTable
                            );
                    FileUtils.writeStringToFile(outFile, page, "UTF-8");
                    // Finally, create our summary row.  All numbers are expressed as a percent of the reference genome peg count.
                    double denom = refGenome.getPegs().size() / 100.0;
                    DomContent sumRow = tr(Html.colorCell(gReport.isGood(), newGenome.genomeLink().render()),
                            td(Html.gPageLink(genomeId, newGenome.getName())),
                            td(refGenome.genomeLink()),
                            td(refGenome.getName()),
                            Html.numCell(this.compareObj.getIdentical() / denom),
                            Html.numCell((this.compareObj.getIdentical() + this.compareObj.getLonger() + this.compareObj.getShorter()) / denom),
                            Html.numCell(this.compareObj.getDifferentFunctions() / denom),
                            Html.numCell(this.compareObj.getNewOnlyCount() / denom),
                            Html.numCell(this.compareObj.getOldOnlyCount() / denom)
                            );
                    this.summaryRows.add(sumRow);
                }
            }
        }
    }

    @Override
    protected void writeSummary(GenomeStats gReport) throws IOException {
    }

    @Override
    protected void endSummary() throws IOException {
        String page = Html.page("Evaluation Summary Report",
                    h1("Comparison Summary Report"),
                    p(String.format("%d genomes processed using evaluator version %s.",
                            this.getGenomeCount(), this.getVersion())),
                    p("All numbers shown in the ORF Comparison Results table are percentages relative to the number of PEGs in the reference genome."),
                    Html.formatTable("ORF Comparison Results", this.summaryRows)
                );
        File summaryFile = new File(this.getOutDir(), "index.html");
        FileUtils.writeStringToFile(summaryFile, page, "UTF-8");
    }

    @Override
    protected void finish() {
    }

    @Override
    public void setupGenomes(GenomeStats[] reports) {
        this.refEngine.setupReferences(reports);
    }

}
