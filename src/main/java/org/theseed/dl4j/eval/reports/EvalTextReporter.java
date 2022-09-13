/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;

import org.theseed.dl4j.eval.stats.GenomeAnalysis;
import org.theseed.dl4j.eval.stats.GenomeStats;

/**
 * Here the reports are being produced in tab-delimited files.  The summary report name is "summary" and all files have an extension of ".tsv"
 * so they can be opened in Excel on the Mac.
 *
 * @author Bruce Parrello
 *
 */
public class EvalTextReporter extends EvalReporter {

    private static final String SUMMARY_FORMAT = "%s\t%s\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8d\t%-8s\t%-8s\t%8.2f\t%s%n";
    private static final String END_FORMAT     = "%s\t%s\t%s\t%8d\t%8d\t%8d\t%8d\t%s\t%8d\t%8d\t%8d\t%8d%n";

    // FIELDS
    /** output stream for the summary report */
    PrintWriter summaryStream;

    @Override
    protected void initialize(File modelDir) throws IOException {
        // Denote we have no summary output file.
        this.summaryStream = null;
    }

    @Override
    protected void startSummary() throws IOException {
        // Open the summary file.
        File outFile = new File(this.getOutDir(), "summary.tsv");
        this.summaryStream = new PrintWriter(outFile);
        // Produce the summary header.
        this.summaryStream.println("Genome\tName\tCoarse\tFine\tCompleteness\tContamination\tHypothetical\tContigs\tGood_Seed\tSSU rRNA\tScore\tGood");
    }

    @Override
    protected void writeDetails(GenomeStats gReport, GenomeAnalysis analysis) throws IOException {
        String genome = gReport.getId();
        File outFile = new File(this.getOutDir(), genome + ".tsv");
        try (PrintWriter genomeStream = new PrintWriter(outFile)) {
            genomeStream.println("Role\tactual\tpredicted\tuniversal");
            Collection<String> pprs = gReport.getProblematicRoles();
            for (String role : pprs) {
                GenomeStats.ProblematicRole ppr = gReport.getReport(role);
                String uFlag = (ppr.isUniversal() ? "Y" : "");
                genomeStream.format("%s\t%d\t%d\t%s%n", role, ppr.getActual(), ppr.getPredicted(), uFlag);
            }
            genomeStream.println();
            genomeStream.format("*\tCompleteness Group\t%s\t%n", gReport.getGroup());
            genomeStream.format("*\tCompleteness\t%2.2f\t%n", gReport.getCompletePercent());
            genomeStream.format("*\tContamination\t%2.2f\t%n", gReport.getContaminationPercent());
            genomeStream.format("*\tCoarse Consistency\t%2.2f\t%n", gReport.getCoarsePercent());
            genomeStream.format("*\tFine Consistency\t%2.2f\t%n", gReport.getFinePercent());
            genomeStream.format("*\tHypothetical Rate\t%2.2f\t%n", gReport.getHypotheticalPercent());
            genomeStream.format("*\tPLFAM Coverage\t%2.2f\t%n", gReport.getPlfamPercent());
            genomeStream.format("*\tCDS Coverage\t%2.2f\t%n", gReport.getCdsPercent());
            genomeStream.format("*\tDNA Length\t%d\t%n", gReport.getDnaSize());
            genomeStream.format("*\tContig Count\t%d\t%n", gReport.getContigCount());
            genomeStream.format("*\tHypothetical Roles\t%d\t%n", gReport.getHypoCount());
            genomeStream.format("*\tContig Length L50\t%d\t%n", gReport.getL50());
            genomeStream.format("*\tContig Length N50\t%d\t%n", gReport.getN50());
            genomeStream.format("*\tCDS Count\t%d\t%n", gReport.getPegCount());
            genomeStream.format("*\tPLFAM Protein Count\t%d\t%n", gReport.getPlfamCount());
            genomeStream.format("*\tGood PheS found\t%s\t%n", (gReport.isGoodSeed() ? "Yes" : "No"));
            genomeStream.format("*\tSSU rRNA Sequence Found\t%s\t%n", (gReport.hasSsuRRna() ? "Yes" : "No"));
            genomeStream.format("*\tGenome Rating\t%s\t%n", (gReport.isGood() ? "Good" : "Poor"));
        }
    }

    @Override
    protected void writeSummary(GenomeStats gReport) throws IOException {
        String genome = gReport.getId();
        String gName = gReport.getName();
        double coarsePct = gReport.getCoarsePercent();
        double finePct = gReport.getFinePercent();
        double completePct = gReport.getCompletePercent();
        double contamPct = gReport.getContaminationPercent();
        double hypoPct = gReport.getHypotheticalPercent();
        int contigs = gReport.getContigCount();
        double score = gReport.getScore();
        String goodSeed = (gReport.isGoodSeed() ? "Y" : "");
        String goodSSU = (gReport.hasSsuRRna() ? "Y" : "");
        String goodGenome = (gReport.isGood() ? "Y" : "");
        this.summaryStream.format(SUMMARY_FORMAT,
                genome, gName, coarsePct, finePct, completePct, contamPct, hypoPct, contigs, goodSeed, goodSSU, score, goodGenome);
    }

    @Override
    protected void endSummary() throws IOException {
        // Print a spacer and format a total line.
        this.summaryStream.println();
        this.summaryStream.format(END_FORMAT, "Total", "Good in each Category", "", this.getConsistentCount(),
                this.getCompleteCount(), this.getCleanCount(), this.getUnderstoodCount(), "",
                this.getGoodSeedCount(), this.getSsuFoundCount(), this.getGenomeCount(), this.getGoodCount());
    }

    @Override
    protected void finish() {
        if (this.summaryStream != null)
            this.summaryStream.close();
    }


}
