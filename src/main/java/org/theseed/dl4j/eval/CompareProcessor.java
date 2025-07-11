/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.compare.CompareORFs;
import org.theseed.genome.compare.CompareFeatures;
import org.theseed.locations.Location;
import org.theseed.reports.Html;
import org.theseed.sequence.MD5Hex;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * This class produces a web site showing an ORF-by-ORF comparison of genomes with identical contigs.
 * The positional parameters are the name of an input directory containing GTOs, a reference directory
 * containing the GTOs for comparison, and the name of an output directory.  For each input genome,
 * we will find a matching reference genome and then produce a report comparing the ORFs.
 *
 * The command-line options are as follows.
 *
 * -v	produce more detailed status output
 *
 * --clear	erase output directory before starting
 *
 * @author Bruce Parrello
 *
 */
public class CompareProcessor extends BaseProcessor {

    // FIELDS

    /** map of reference genomes, MD5 to file name */
    private Map<String, File> refGenomes;
    /** comparison processor */
    private CompareORFs comparator;
    /** map of scores to summary table rows with that score, used to sort the output */
    private TreeMap<Double, Collection<DomContent>> summaryData;
    /** list of detail table rows */
    private List<DomContent> tableRows;
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(Evaluator.class);
    /** first genome for comparison */
    private Genome genome1;
    /** second genome for comparison */
    private Genome genome2;
    /** genome identity computer */
    private MD5Hex md5Computer;


    // COMMAND LINE OPTIONS

    /** clear-output flag */
    @Option(name = "--clear", usage = "clear output directory before starting")
    private boolean clearOutDir;

    /** input genome directory */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory", required = true)
    private File inDir;

    /** reference genome directory */
    @Argument(index = 1, metaVar = "refDir", usage = "comparison genome directory (if not input directory)")
    private File refDir;

    /** output directory */
    @Argument(index = 2, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    @Override
    public boolean validateParms() throws IOException {
        if (! this.inDir.isDirectory()) {
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        } else if (! this.refDir.isDirectory()) {
            throw new FileNotFoundException("Reference input directory " + this.refDir + " is not found or invalid.");
        } else {
            if (! this.outDir.exists()) {
                log.info("Creating output directory {}.", this.outDir);
                if (! this.outDir.mkdir())
                    throw new IOException("Could not create output directory.");
            } else if (! this.outDir.isDirectory()) {
                throw new IOException("Output directory " + this.outDir + " is invalid.");
            } else if (this.clearOutDir) {
                log.info("Erasing output directory {}.", this.outDir);
                FileUtils.cleanDirectory(this.outDir);
            } else {
                log.info("Output will be in {}.", this.outDir);
            }
        }
        return true;
    }

    @Override
    public void setDefaults() {
        this.clearOutDir = false;
        this.refDir = null;
    }

    @Override
    public void runCommand() throws Exception {
        // Create the comparator and the MD5 calculator.
        this.comparator = new CompareFeatures();
        this.md5Computer = new MD5Hex();
        // Initialize the summary row table.  Note the header is given a low score to sort it to the top.
        this.summaryData = new TreeMap<Double, Collection<DomContent>>();
        this.summaryData.put(-1.0, Collections.singleton(tr(th("Match Score").withClass("num"), th("Genome 1 ID"), th("Genome 1 Name"), th("Genome 2 ID"), th("Genome 2 Name"))));
        // Get all the reference genomes.  We would like to hold them in memory, but it is too much.
        // Instead, we map the genome MD5 to its file name.
        log.info("Analyzing reference genomes from {}.", this.refDir);
        this.refGenomes = this.comparator.getMd5GenomeMap(this.refDir);
        // Get all the genomes.
        log.info("Reading genomes from {}.", this.inDir);
        GenomeDirectory genomesIn = new GenomeDirectory(this.inDir);
        for (Genome genome : genomesIn) {
            log.info("Processing genome {}.", genome);
            // Search for a match.
            String key = this.md5Computer.sequenceMD5(genome);
            File refGenomeFile = this.refGenomes.get(key);
            if (refGenomeFile == null) {
                log.warn("No match found for {}", genome);
            } else {
                // Read the reference genome.
                this.genome1 = genome;
                this.genome2 = new Genome(refGenomeFile);
                // Map the contigs.
                mapContigs();
                // Do a compare.
                this.compareReport();
            }
        }
        log.info("Writing summary page.");
        Collection<DomContent> summaryRows = new ArrayList<DomContent>(genomesIn.size());
        for (Collection<DomContent> summaryEntry : summaryData.values())
            summaryRows.addAll(summaryEntry);
        String page = Html.page("Summary of ORF Comparisons", Html.formatTable("Genome Comparisons", summaryRows));
        File outFile = new File(this.outDir, "index.html");
        FileUtils.writeStringToFile(outFile, page, "UTF-8");
    }

    /**
     * Change the name of the reference genome (genome1) contig IDs in the features to match the IDs of the contigs
     * in the target genome (genome2).
     *
     * @throws UnsupportedEncodingException
     */
    private void mapContigs() throws UnsupportedEncodingException {
        // Build a map of MD5s to IDs for each contig in the target genome.
        Map<String,String> md5Map = new HashMap<String, String>(this.genome2.getContigCount());
        for (Contig contig2 : genome2.getContigs()) {
            String contigKey = this.md5Computer.sequenceMD5(contig2.getSequence());
            md5Map.put(contigKey, contig2.getId());
        }
        // Now map each contig ID in the reference genome to the appropriate target genome contig ID.
        for (Contig contig1 : genome1.getContigs()) {
            String contigKey = this.md5Computer.sequenceMD5(contig1.getSequence());
            String contig2Id = md5Map.get(contigKey);
            // The two genomes are DNA-identical, so there should always be a matching contig ID.
            assert(contig2Id != null);
            genome1.updateContigId(contig1, contig2Id);
        }
    }

    /**
     * Produce a report comparing the two current genomes.
     *
     * @throws IOException
     */
    private void compareReport() throws IOException {
        // Start the table.
        this.tableRows = new ArrayList<DomContent>(3000);
        this.tableRows.add(tr(th("Location"), th("Peg 1"), th("Function 1"), th("Peg 2"), th("Function 2"), th("Diff").withClass("num")));
        // Get both sets of features in ORF order.
        log.info("Sorting genome features into ORFs for {} and {}.", genome1.getId(), genome2.getId());
        SortedSet<Feature> f1List = this.comparator.sortFeatures(genome1);
        SortedSet<Feature> f2List = this.comparator.sortFeatures(genome2);
        // These are used to compute the match score.
        int matches = 0;
        int total = 0;
        int extraLeft = 0;
        int extraRight = 0;
        // Get iterators through the sets.
        Iterator<Feature> iter1 = f1List.iterator();
        Iterator<Feature> iter2 = f2List.iterator();
        Feature f1 = this.comparator.next(iter1);
        Feature f2 = this.comparator.next(iter2);
        while (f1 != null && f2 != null) {
            // Compare these features.
            int comp = this.comparator.orfCompare(f1, f2);
            if (comp < 0) {
                // F1 is first.  It is an orphan.
                this.tableRow(f1, null);
                f1 = this.comparator.next(iter1);
                extraLeft++;
                total++;
            } else if (comp > 0) {
                // F2 is first.  It is an orphan.
                this.tableRow(null, f2);
                f2 = this.comparator.next(iter2);
                extraRight++;
                total++;
            } else {
                // Both are in the same ORF.
                this.tableRow(f1, f2);
                f1 = this.comparator.next(iter1);
                f2 = this.comparator.next(iter2);
                matches++;
                total++;
            }
        }
        // Run out the residuals.
        while (f1 != null) {
            this.tableRow(f1, null);
            f1 = this.comparator.next(iter1);
            extraLeft++;
            total++;
        }
        while (f2 != null) {
            this.tableRow(null, f2);
            f2 = this.comparator.next(iter2);
            extraRight++;
            total++;
        }
        // Create the web page.
        log.info("Creating web page.");
        double score = (double) (matches * 100) / total;
        String title = "ORF Comparison of " + genome1.getId() + " and " + genome2.getId();
        String page = Html.page(title, h1(title),
                ul(li(join("Genome 1 is", genome1.genomeLink(), genome1.getName(),
                        String.format("and has %d extra proteins.", extraLeft))),
                li(join("Genome 2 is ", genome2.genomeLink(), genome2.getName(),
                        String.format("and has %d extra proteins.", extraRight))),
                li(String.format("%4.2f percent of the ORFs matched out of %d total.", score, total))),
                Html.formatTable("Content of each ORF", tableRows));
        String fileName = genome1.getId() + ".html";
        File outFile = new File(this.outDir, fileName);
        FileUtils.writeStringToFile(outFile, page, "UTF-8");
        // Create the summary line and add it to the summary rows.
        DomContent scoreCell = td(a(String.format("%8.2f", score)).withHref(fileName)).withClass("num");
        DomContent summaryRow = tr(scoreCell, td(genome1.genomeLink()), td(genome1.getName()), td(genome2.genomeLink()), td(genome2.getName()));
        this.summaryData.computeIfAbsent(score, k -> new ArrayList<DomContent>()).add(summaryRow);
    }

    /**
     * Add a detail row to the detail table.
     *
     * @param f1	first genome feature, or NULL if the second genome feature is an orphan
     * @param f2	second genome feature, or NULL if the first genome feature is an orphan
     */
    private void tableRow(Feature f1, Feature f2) {
        if (f1 == null && f2 == null)
            throw new IllegalArgumentException("Cannot create table row with two null features.");
        ContainerTag cell1 = (f1 != null ? td(genome1.featureRegionLink(f1.getId())) : Html.emptyCell());
        ContainerTag cell2 = (f1 != null ? td(f1.getFunction()) : Html.emptyCell());
        ContainerTag cell3 = (f2 != null ? td(genome2.featureRegionLink(f2.getId())) : Html.emptyCell());
        ContainerTag cell4 = (f2 != null ? td(f2.getFunction()) : Html.emptyCell());
        if (f1 == null) {
            cell4 = cell4.withStyle(Html.BAD_STYLE);
        } else if (f2 == null) {
            cell2 = cell2.withStyle(Html.BAD_STYLE);
        } else if (! StringUtils.equalsIgnoreCase(f1.getFunction(), f2.getFunction())) {
            cell2 = cell2.withStyle(Html.BAD_STYLE);
            cell4 = cell4.withStyle(Html.BAD_STYLE);
        }
        // Compute the length difference.
        int len = (f2 != null ? f2.getProteinLength() : 0) - (f1 != null ? f1.getProteinLength() : 0);
        ContainerTag cell5 = Html.colorCell(len == 0, len);
        @SuppressWarnings("null") // because we checked for both null above
        Location loc = (f1 != null ? f1.getLocation() : f2.getLocation());
        ContainerTag row = tr(td(loc.toString()), cell1, cell2, cell3, cell4, cell5);
        this.tableRows.add(row);
    }

}
