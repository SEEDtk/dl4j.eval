/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.dl4j.eval.stats.QualityKeys;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.reports.Html;
import org.theseed.reports.TextBuilder;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

import com.github.cliftonlabs.json_simple.JsonObject;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;

import static j2html.TagCreator.*;

/**
 * This command produces the binning summary report.  It takes as input a control file describing the the bin GTO files and
 * the URL of each bin's genome report.  It outputs an HTML summary report with links to the bin genome reports.
 *
 * The control file is on the standard input.  It should be tab-delimited with headers, the first column containing
 * file names and the second column the corresponding URLs.  The report will be printed on the standard output in HTML format.
 *
 * There are no positional parameters.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input control file (if not STDIN)
 * -o	output HTML file (if not STDOUT)
 *
 * --job		binning job ID (if known)
 * --group		name of the output genome group (if any)
 *
 * @author Bruce Parrello
 *
 */
public class BinReportProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinReportProcessor.class);

    // COMMAND-LINE OPTIONS

    /** job ID */
    @Option(name = "--job", metaVar = "123456", usage = "binning job ID number (if known)")
    private String jobID;

    /** group ID */
    @Option(name = "--group", metaVar = "BinGroup", usage = "output group for bin genomes")
    private String groupName;

    /**
     * Genome status enum.
     */
    public static enum QualityStatus {
        GOOD("good", "Good Quality Genomes"),
        MOSTLY_GOOD("mostly good", "Generally Good-Quality Genomes with Missing SSU rRNA"),
        POOR("poor", "Poor Quality Genomes");

        /** status name (short) */
        private String name;
        /** status description (long) */
        private String description;

        /**
         * Create a quality status.
         *
         * @param name			short and readable status name
         * @param description	long status description
         */
        private QualityStatus(String name, String description) {
            this.name = name;
            this.description = description;
        }

        /**
         * @return the status description
         */
        public String getDescription() {
            return this.description;
        }

        /**
         * @return a short and readable status name
         */
        public String getShortName() {
            return this.name;
        }

    }

    /**
     * Genome descriptor from the input file.  It compares by status first (best to worst), then score.
     */
    private static class GenomeData implements Comparable<GenomeData> {

        /** ID of the genome */
        private String genomeId;
        /** genome quality score */
        private double score;
        /** genome quality status */
        private QualityStatus status;
        /** genome output table row */
        private ContainerTag tableRow;
        /** table header, matching the rows */
        protected static final ContainerTag TABLE_HEADER = tr(
                th("Score").withClass("num"),
                th("Genome ID"),
                th("Genome Name"),
                th("Ref Genome"),
                th("Coverage").withClass("num"),
                th("Coarse %").withClass("num"),
                th("Fine %").withClass("num"),
                th("Complete %").withClass("num"),
                th("Contamination %").withClass("num"),
                th("Hypothetical %").withClass("num"),
                th("Contigs").withClass("num"),
                th("Base Pairs").withClass("num"),
                th("Good Phes").withClass("flag"),
                th("SSU rRNA").withClass("flag")
            );

        /**
         * Construct the genome data object from the genome and the URL.
         *
         * @param genome	genome to be included in the report
         * @param url		URL for the genome report
         */
        protected GenomeData(Genome genome, String url) {
            this.genomeId = genome.getId();
            // Get the quality structure.
            JsonObject quality = genome.getQuality();
            // Compute the status.
            if (! quality.getBooleanOrDefault(QualityKeys.MOSTLY_GOOD))
                this.status = QualityStatus.POOR;
            else if (quality.getBooleanOrDefault(QualityKeys.HAS_SSU_RNA))
                this.status = QualityStatus.GOOD;
            else
                this.status = QualityStatus.MOSTLY_GOOD;
            // Save the score.
            this.score = quality.getDoubleOrDefault(QualityKeys.SCORE);
            // Now we must build the summary table row.  First comes the score.
            DomContent score = a(Html.num(this.score)).withHref(url);
            // For the genome ID, we need to build a link.
            DomContent genomeId = genome.genomeLink().withTarget("_blank");
            // The reference genome ID may not exist, so it is extra tricky.
            var linker = genome.getLinker();
            String refGenomeId = quality.getStringOrDefault(QualityKeys.BIN_REF_GENOME);
            DomContent refCell;
            if (refGenomeId == null)
                refCell = Html.emptyCell();
            else
                refCell = td(linker.genomeLink(refGenomeId));
            // Get the scores that decide the quality.
            double fineConsistency = quality.getDoubleOrDefault(QualityKeys.FINE_CONSISTENCY);
            double completeness = quality.getDoubleOrDefault(QualityKeys.COMPLETENESS);
            double contamination = quality.getDoubleOrDefault(QualityKeys.CONTAMINATION);
            double hypoPercent = quality.getDoubleOrDefault(QualityKeys.HYPOTHETICAL_CDS_RATIO);
            // Now we assemble the table row.
            this.tableRow = tr(td(score),
                    td(genomeId), td(genome.getName()), refCell,
                    Html.numCell(quality.getDoubleOrDefault(QualityKeys.BIN_COVERAGE)),
                    Html.numCell(quality.getDoubleOrDefault(QualityKeys.COARSE_CONSISTENCY)),
                    Html.colorCell(GenomeStats.indicatesConsistent(fineConsistency), fineConsistency),
                    Html.colorCell(GenomeStats.indicatesComplete(completeness), completeness),
                    Html.colorCell(GenomeStats.indicatesClean(contamination), contamination),
                    Html.colorCell(GenomeStats.indicatesUnderstood(hypoPercent), hypoPercent),
                    Html.numCell(genome.getContigCount()),
                    Html.numCell(genome.getLength()),
                    Html.flagCell(quality.getBooleanOrDefault(QualityKeys.HAS_SEED), "Y", ""),
                    Html.flagCell(quality.getBooleanOrDefault(QualityKeys.HAS_SSU_RNA), "Y", "")
                );
        }

        /**
         * @return the heading row for an output table
         */
        protected static ContainerTag getHeader() {
            return TABLE_HEADER;
        }

        @Override
        public int compareTo(GenomeData o) {
            int retVal = this.status.compareTo(o.status);
            if (retVal == 0) {
                retVal = Double.compare(o.score, this.score);
                if (retVal == 0)
                    retVal = this.genomeId.compareTo(o.genomeId);
            }
            return retVal;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.genomeId == null) ? 0 : this.genomeId.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof GenomeData)) {
                return false;
            }
            GenomeData other = (GenomeData) obj;
            if (this.genomeId == null) {
                if (other.genomeId != null) {
                    return false;
                }
            } else if (!this.genomeId.equals(other.genomeId)) {
                return false;
            }
            return true;
        }

        /**
         * @return the status
         */
        public QualityStatus getStatus() {
            return this.status;
        }

        /**
         * @return the HTML table row
         */
        public ContainerTag getTableRow() {
            return this.tableRow;
        }


    }

    @Override
    protected void setPipeDefaults() {
        this.jobID = null;
        this.groupName = null;
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Set up some counters.
        int linesIn = 0;
        int badFiles = 0;
        int bins = 0;
        // Create table row storage for each type of bin.
        var rowSetMap = new EnumMap<QualityStatus, Set<GenomeData>>(QualityStatus.class);
        // Loop through the input stream.  We process each genome and put its descriptor in the appropriate row set.
        for (var line : inputStream) {
            File gFile = new File(line.get(0));
            String url = line.get(1);
            linesIn++;
            if (! gFile.canRead()) {
                log.error("Genome file {} is not found or unreadable.", gFile);
                badFiles++;
            } else {
                Genome genome = new Genome(gFile);
                log.info("Genome {} with url {} loaded.", genome, url);
                var genomeData = new GenomeData(genome, url);
                // Store it in the row-set map.
                Set<GenomeData> rowSet = rowSetMap.computeIfAbsent(genomeData.getStatus(), x -> new TreeSet<GenomeData>());
                rowSet.add(genomeData);
                bins++;
            }
        }
        log.info("{} lines in, {} processed, {} bad files.", linesIn, bins, badFiles);
        // We now need to form the output page, which will then be unspooled to the report writer.
        // The output page has a section for each quality status and a header phrase containing the
        // counts.
        List<DomContent> sections = new ArrayList<DomContent>(QualityStatus.values().length + 2);
        // Compute the page title.
        String title;
        if (this.jobID != null)
            title = "Binning Job " + this.jobID;
        else
            title = "Binning Results";
        sections.add(h1(title));
        // Check for a failure.
        if (bins == 0)
            sections.add(p("No bins found."));
        else {
            // Next are the counts.  We add an empty placeholder tag so we can build the other sections in
            // the same loop.
            ContainerTag counts = ul();
            sections.add(counts);
            // This buffer will hold the count sentences.
            var buffer = new TextBuilder(40);
            // Put in the total count.
            buffer.append(bins, "bin", "bins").append(" found.");
            counts.with(li(buffer.toString()));
            // Add some text about the quality.
            for (var sectionData : rowSetMap.entrySet()) {
                var status = sectionData.getKey();
                var rows = sectionData.getValue();
                final int count = rows.size();
                log.info("{} genomes in {} section.", count, status);
                // Build the table rows.
                var htmlRows = new ArrayList<DomContent>(count + 1);
                htmlRows.add(GenomeData.getHeader());
                sectionData.getValue().stream().forEach(x -> htmlRows.add(x.getTableRow()));
                // Create the count sentence.  This will link to the appropriate section.
                buffer.clear();
                buffer.append(count, "is", "are").append(" ").append(status.getShortName()).append(".");
                String name = status.name();
                counts.with(li(a(buffer.toString()).withHref("#" + name)));
                // Build the table for this status.
                DomContent table = Html.formatTable(status.getDescription(), htmlRows);
                // Bracket the table with the name anchor.
                sections.add(a().withName(name));
                sections.add(table);
            }
            // Add a help sentence about the quality.
            counts.with(li(join(i("Mostly-good genomes"),
                    String.format("have completeness >= %2.0f, fine consistency >= %2.0f, contamination <= %2.0f, hypothetical percent <= %2.0f",
                    GenomeStats.MIN_COMPLETENESS, GenomeStats.MIN_CONSISTENCY, GenomeStats.MAX_CONTAMINATION,
                    GenomeStats.MAX_HYPOTHETICAL), "and a single PheS protein of reasonable size.")));
            counts.with(li(join("A mostly-good genome with an identifiable SSU rRNA is classified as ", i("good"), ".")));
            // Display the genome group, if any.
            if (this.groupName != null)
                counts.with(li(join("Genomes have been added to group", b(this.groupName), ".")));
        }
        // Format the output page.
        DomContent[] sectionArray = sections.stream().toArray(DomContent[]::new);
        String page = Html.page(title, sectionArray);
        // Write it to the report output.
        writer.print(page);
    }

}
