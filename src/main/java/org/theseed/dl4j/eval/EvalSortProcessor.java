/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command sorts the output of a mass evaluation so it can be fed to the various update
 * commands.  The positional parameters are the output file prefix and the name of the output
 * directory.  The output directory will contain three files.
 *
 * 	XXXXXX.sort.tbl		which is the sorted output
 *  XXXXXX.good.tbl		which contains only the good genomes
 *  XXXXXX.stat.tbl		which contains counts related to the genomes
 *
 * where XXXXXX is the output prefix.  The unsorted report should come in via the standard
 * input.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	if specified, the name of the input file (default STDIN)
 *
 * --prio	if specified, the name of a tab-delimited file (with headers) containing genome IDs
 * 			in the first column; the identified genomes will be sorted to the top of the
 * 			good-genome section
 *
 * @author Bruce Parrello
 *
 */
public class EvalSortProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(EvalSortProcessor.class);
    /** input stream */
    private TabbedLineReader inStream;
    /** priority genome ID set */
    private Set<String> priorityGenomes;
    /** genome counts */
    private CountMap<String> statMap;
    /** list of stat types */
    private static final Set<String> STAT_TYPES = Arrays.asList("HasRNA", "HasPheS", "Clean", "Complete",
            "Consistent", "Understood").stream().collect(Collectors.toSet());

    // COMMAND-LINE OPTONS

    /** input file (if not STDIN) */
    @Option(name = "-i", aliases = { "--input" }, usage = "input file (if not STDIN")
    private File inFile;

    /** priority genome file */
    @Option(name = "--prio", metaVar = "genomes.tbl", usage = "file containg priority genome IDs in first column")
    private File prioFile;

    /** output file name prefix */
    @Argument(index = 0, metaVar = "prefix", usage = "prefix to be put on output file names", required = true)
    private String prefix;

    /** output directory */
    @Argument(index = 1, metaVar = "outDir", usage = "output file directory")
    private File outDir;

    /**
     * This class provides an ordering for the input lines.  We sort by good before poor, then priority
     * before non-priority, then descending score, then genome ID.
     */
    private class LineSorter implements Comparator<TabbedLineReader.Line> {

        @Override
        public int compare(TabbedLineReader.Line o1, TabbedLineReader.Line o2) {
            // Prioritize good over poor.
            boolean g1 = o1.getFancyFlag(MassEvalProcessor.GOOD_COL);
            boolean g2 = o2.getFancyFlag(MassEvalProcessor.GOOD_COL);
            int retVal = Boolean.compare(g2, g1);
            if (retVal == 0) {
                // Get the genome IDs.
                String id1 = o1.get(MassEvalProcessor.GENOME_COL);
                String id2 = o2.get(MassEvalProcessor.GENOME_COL);
                // Determine the priorities. Priority genomes go first.
                boolean p1 = EvalSortProcessor.this.priorityGenomes.contains(id1);
                boolean p2 = EvalSortProcessor.this.priorityGenomes.contains(id2);
                retVal = Boolean.compare(p2, p1);
                if (retVal == 0) {
                    // Sort by descending score.
                    retVal = Double.compare(o2.getDouble(MassEvalProcessor.SCORE_COL),
                            o1.getDouble(MassEvalProcessor.SCORE_COL));
                    // Break ties on genome ID.
                    if (retVal == 0)
                        retVal = id1.compareTo(id2);
                }
            }
            return retVal;
        }

    }

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.prioFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else
            log.info("Output will be in {}.", this.outDir);
        // Set up the priority genome set.
        if (this.prioFile == null)
            this.priorityGenomes = Collections.emptySet();
        else {
            this.priorityGenomes = TabbedLineReader.readSet(this.prioFile, "1");
            log.info("{} genomes found in priority file.", this.priorityGenomes.size());
        }
        // Connect the input file.
        if (this.inFile == null) {
            log.info("Genome data will be taken from the standard input.");
            this.inStream = new TabbedLineReader(System.in);
        } else {
            log.info("Genome data will be read from {}.", this.inFile);
            this.inStream = new TabbedLineReader(this.inFile);
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        try {
            // We will accumulate the input lines in here in order to sort them.
            SortedSet<TabbedLineReader.Line> lineSorter = new TreeSet<TabbedLineReader.Line>(this.new LineSorter());
            // We will count the various genome attributes in here.
            this.statMap = new CountMap<String>();
            // Loop through the input file.
            log.info("Scanning input file.");
            for (TabbedLineReader.Line line : this.inStream) {
                lineSorter.add(line);
                // Here we need to analyze the line to compute the statistics.
                Set<String> stats = new TreeSet<String>();
                if (line.getFancyFlag(MassEvalProcessor.GOOD_COL))
                    this.statMap.count("Good");
                else
                    this.statMap.count("Bad");
                if (line.getFancyFlag(MassEvalProcessor.SSU_RNA_COL))
                    stats.add("HasRNA");
                else
                    this.statMap.count("MissingRNA");
                if (line.getFancyFlag(MassEvalProcessor.SEED_COL))
                    stats.add("HasPheS");
                else
                    this.statMap.count("MissingPheS");
                if (GenomeStats.indicatesClean(line.getDouble(MassEvalProcessor.CONTAM_COL)))
                    stats.add("Clean");
                else
                    this.statMap.count("Contaminated");
                if (GenomeStats.indicatesComplete(line.getDouble(MassEvalProcessor.COMPLETE_COL)))
                    stats.add("Complete");
                else
                    this.statMap.count("Incomplete");
                if (GenomeStats.indicatesConsistent(line.getDouble(MassEvalProcessor.FINE_COL)))
                    stats.add("Consistent");
                else
                    this.statMap.count("Inconsistent");
                if (GenomeStats.indicatesUnderstood(line.getDouble(MassEvalProcessor.HYPO_COL)))
                    stats.add("Understood");
                else
                    this.statMap.count("Misunderstood");
                // If there is a single failure point, note that.
                if (STAT_TYPES.size() - stats.size() == 1) {
                    for (String reason : STAT_TYPES) {
                        if (! stats.contains(reason))
                            this.statMap.count("Not " + reason);
                    }
                }
            }
            // Now we have processed all the lines, so we produce the output files.  Start with the line
            // output files.
            log.info("{} lines read.  Producing output.", lineSorter.size());
            File allFile = new File(this.outDir, this.prefix + ".sort.tbl");
            File goodFile = new File(this.outDir, this.prefix + ".good.tbl");
            try (PrintWriter allWriter = new PrintWriter(allFile);
                    PrintWriter goodWriter = new PrintWriter(goodFile)) {
                String header = StringUtils.join(MassEvalProcessor.DEFAULT_HEADERS, '\t');
                allWriter.println(header);
                goodWriter.println(header);
                for (TabbedLineReader.Line line : lineSorter) {
                    allWriter.println(line);
                    if (line.getFancyFlag(MassEvalProcessor.GOOD_COL))
                        goodWriter.println(line);
                }
            }
            File statFile = new File(this.outDir, this.prefix + ".stats.tbl");
            try (PrintWriter statWriter = new PrintWriter(statFile)) {
                statWriter.println("status\tcount");
                for (CountMap<String>.Count counter : this.statMap.sortedCounts())
                    statWriter.format("%s\t%d%n", counter.getKey(), counter.getCount());
            }
        } finally {
            this.inStream.close();
        }
    }

}
