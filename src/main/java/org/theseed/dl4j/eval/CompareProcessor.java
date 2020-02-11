/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Compare;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.locations.Location;
import org.theseed.reports.Html;
import org.theseed.utils.ICommand;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.LoggerContext;
import j2html.tags.ContainerTag;
import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * This class produces a web site showing an ORF-by-ORF comparison of genomes with identical contigs.
 * The positional parameters are the name of an input directory containing GTOs and the name of an
 * output directory.  The GTOs will be loaded into memory and then processed individually.  For each,
 * we will look for a genome that has identical contigs and then produce a report comparing the
 * ORFs.
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
public class CompareProcessor implements ICommand {

    // FIELDS

    /** list of genomes in memory */
    private List<Genome> genomes;
    /** comparison processor */
    private Compare comparator;
    /** list of summary table rows */
    private List<DomContent> summaryRows;
    /** list of detail table rows */
    private List<DomContent> tableRows;
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(Evaluator.class);
    /** first genome for comparison */
    private Genome genome1;
    /** second genome for comparison */
    private Genome genome2;


    // COMMAND LINE OPTIONS

    /** help option */
    @Option(name = "-h", aliases = { "--help" }, help = true)
    private boolean help;

    /** clear-output flag */
    @Option(name = "--clear", usage = "clear output directory before starting")
    private boolean clearOutDir;

    /** debug-message flag */
    @Option(name = "-v", aliases = { "--verbose", "--debug" }, usage = "show more detailed progress messages")
    private boolean debug;

    @Argument(index = 0, usage = "input directory", required = true)
    private File inDir;

    @Argument(index = 1, usage = "output directory", required = true)
    private File outDir;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.debug = false;
        this.clearOutDir = false;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                if (! this.inDir.isDirectory()) {
                    throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
                } else {
                    if (! this.outDir.exists()) {
                        log.info("Creating output directory {}.", this.outDir);
                        if (! this.outDir.mkdir())
                            throw new IOException("Could not create output directory.");
                        else
                            retVal = true;
                    } else if (! this.outDir.isDirectory()) {
                        throw new IOException("Output directory " + this.outDir + " is invalid.");
                    } else if (this.clearOutDir) {
                        log.info("Erasing output directory {}.", this.outDir);
                        FileUtils.cleanDirectory(this.outDir);
                        retVal = true;
                    } else {
                        log.info("Output will be in {}.", this.outDir);
                        retVal = true;
                    }
                }
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
        return retVal;
    }

    @Override
    public void run() {
        try {
            if (this.debug) {
                // To get more progress messages, we set the log level in logback.
                LoggerContext loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory();
                ch.qos.logback.classic.Logger logger = loggerContext.getLogger("org.theseed");
                logger.setLevel(Level.toLevel("TRACE"));
            }
            // Create the comparator.
            this.comparator = new Compare();
            // Initialize the summary row table.
            this.summaryRows = new ArrayList<DomContent>();
            this.summaryRows.add(tr(th("Match Score").withClass("num"), th("Genome 1 ID"), th("Genome 1 Name"), th("Genome 2 ID"), th("Genome 2 Name")));
            // Get all the genomes.
            this.genomes = new ArrayList<Genome>();
            log.info("Reading genomes from {}.", this.inDir);
            GenomeDirectory genomesIn = new GenomeDirectory(this.inDir);
            for (Genome genome : genomesIn) {
                String genomeId = genome.getId();
                String genomeName = genome.getName();
                log.info("Processing genome {} {}.", genomeId, genomeName);
                // Search for a match.
                int i = 0;
                while (i < this.genomes.size() && ! this.comparator.compare(genome, this.genomes.get(i))) i++;
                if (i >= this.genomes.size()) {
                    // No match found.  Add this genome to the array.
                    this.genomes.add(genome);
                } else {
                    // Delete the match from the array and do a compare.
                    this.genome1 = this.genomes.remove(i);
                    this.genome2 = genome;
                    this.compareReport();
                }
            }
            log.info("Writing summary page.");
            String page = Html.page("Summary of ORF Comparisons", Html.formatTable("Genome Comparisons", summaryRows));
            File outFile = new File(this.outDir, "index.html");
            FileUtils.writeStringToFile(outFile, page, "UTF-8");
        } catch (Exception e) {
            e.printStackTrace();
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
                total++;
            } else if (comp > 0) {
                // F2 is first.  It is an orphan.
                this.tableRow(null, f2);
                f2 = this.comparator.next(iter2);
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
            total++;
        }
        while (f2 != null) {
            this.tableRow(null, f2);
            f2 = this.comparator.next(iter2);
            total++;
        }
        // Create the web page.
        log.info("Creating web page.");
        double score = (double) (matches * 100) / total;
        String title = "ORF Comparison of " + genome1.getId() + " and " + genome2.getId();
        String page = Html.page(title, h1(title),
                ul(li(join("Genome 1 is", genome1.genomeLink(), genome1.getName())),
                li(join("Genome 2 is ", genome2.genomeLink(), genome2.getName())),
                li(String.format("%4.2f percent of the ORFs matched.", score))),
                Html.formatTable("Content of each ORF", tableRows));
        String fileName = genome1.getId() + ".html";
        File outFile = new File(this.outDir, fileName);
        FileUtils.writeStringToFile(outFile, page, "UTF-8");
        // Create the summary line and add it to the summary rows.
        DomContent scoreCell = td(a(String.format("%8.2f", score)).withHref(fileName)).withClass("num");
        DomContent summaryRow = tr(scoreCell, td(genome1.genomeLink()), td(genome1.getName()), td(genome2.genomeLink()), td(genome2.getName()));
        this.summaryRows.add(summaryRow);
    }

    /**
     * Add a detail row to the detail table.
     *
     * @param f1	first genome feature, or NULL if the second genome feature is an orphan
     * @param f2	second genome feature, or NULL if the first genome feature is an orphan
     */
    private void tableRow(Feature f1, Feature f2) {
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
        Location loc = (f1 != null ? f1.getLocation() : f2.getLocation());
        ContainerTag row = tr(td(loc.toString()), cell1, cell2, cell3, cell4, cell5);
        this.tableRows.add(row);
    }

}
