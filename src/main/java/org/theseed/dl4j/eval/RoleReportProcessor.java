/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

/**
 * This command will create a report of the role occurrences in the genomes found in a genome source.
 * This report can be used by other commands to compute co-occurrence or build a machine learning input
 * matrix.
 *
 * The positional parameters are the name of the role definition file and the name of the genome source.
 * The report will be produced on the standard output.  Each record will contain a genome ID, a role ID,
 * the canonical role name, and the number of occurrences.
 *
 * The command-line options are as follows:
 *
 * -h	command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 * -t	genome source type (default DIR)
 *
 * --filter		name of a tab-delimited file containing the IDs of genomes to include in the first column;
 * 				if omitted, all genomes in the source are included
 *
 * @author Bruce Parrello
 *
 */
public class RoleReportProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RoleReportProcessor.class);
    /** genome source to traverse */
    private GenomeSource genomes;
    /** role definitions */
    private RoleMap roleMap;
    /** set of genome IDs to load */
    private SortedSet<String> genomeIds;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--type", aliases = { "--source", "-t" }, usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** name of optional filter file */
    @Option(name = "--filter", metaVar = "filter.tbl", usage = "list of genome IDs for genomes to use (default is to use all)")
    private File filterFile;

    /** name of the role definition file */
    @Argument(index = 0, metaVar = "roles.in.subsystems", usage = "role definition file", required = true)
    private File roleFile;

    /** name of the input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "input genome source (file or directory)")
    private File inDir;

    @Override
    protected void setReporterDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.filterFile = null;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Process the role file.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " is not found or unreadable.");
        this.roleMap = RoleMap.load(this.roleFile);
        // Validate the filter file.
        if (this.filterFile != null && ! this.filterFile.canRead())
            throw new FileNotFoundException("Filter file " + this.filterFile + " is not found or unreadable.");
        // Load the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Genome source " + this.inDir + " does not exist.");
        this.genomes = this.sourceType.create(this.inDir);
        // Compute the list of genomes to read.  This depends on whether or not we have a filter.
        this.genomeIds = new TreeSet<String>();
        if (this.filterFile != null)
            this.genomeIds.addAll(TabbedLineReader.readSet(this.filterFile, "1"));
        else
            this.genomeIds.addAll(this.genomes.getIDs());
        log.info("{} genomes selected for input,", this.genomeIds.size());
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write the header line.
        writer.println("genome_id\trole_id\tcount\trole_description");
        // Loop through the genomes.
        int genomeCount = 0;
        int skipCount = 0;
        int genomeTotal = this.genomeIds.size();
        for (String genomeId : this.genomeIds) {
            genomeCount++;
            Genome genome = this.genomes.getGenome(genomeId);
            if (genome == null) {
                log.info("Skipping genome {} of {}:  {} not found.", genomeCount, genomeTotal, genomeId);
                skipCount++;
            } else {
                log.info("Processing genome {} of {}:  {}.", genomeCount, genomeTotal, genome);
                // We will count the roles in here.
                CountMap<Role> counts = new CountMap<Role>();
                // Loop through the pegs, collecting roles.
                for (Feature feat : genome.getPegs()) {
                    List<Role> roles = feat.getUsefulRoles(this.roleMap);
                    for (Role role : roles)
                        counts.count(role);
                }
                // Output the role counts.
                for (CountMap<Role>.Count count : counts.sortedCounts()) {
                    Role role = count.getKey();
                    writer.format("%s\t%s\t%d\t%s%n", genomeId, role.getId(), count.getCount(), role.getName());
                }
                log.info("{} roles found in genome.", counts.size());
            }
        }
        log.info("{} genomes skipped out of {}.", skipCount, genomeCount);
    }
}
