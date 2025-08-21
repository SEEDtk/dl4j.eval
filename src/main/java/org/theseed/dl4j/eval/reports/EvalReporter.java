/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;

import org.theseed.dl4j.eval.stats.GenomeAnalysis;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.RoleMap;

/**
 * This is the base class for evaluation report writers.  Each such writer must support methods to
 * open and close the main output file and to write detail reports.
 *
 * @author Bruce Parrello
 *
 */
public abstract class EvalReporter implements AutoCloseable {

    // FIELDS
    /** output directory */
    private File outDir;
    /** TRUE to output summary */
    private boolean summary;
    /** if specified, the name of the output HTML file */
    private String htmlFileName;
    /** TRUE to output details */
    private boolean details;
    /** current version of the evaluation database */
    private String version;
    /** number of clean genomes */
    private int cleanCount;
    /** number of consistent genomes */
    private int consistentCount;
    /** number of complete genomes */
    private int completeCount;
    /** number of understood genomes */
    private int understoodCount;
    /** number of good seeds */
    private int goodSeedCount;
    /** number of SSU rRNAs found */
    private int ssuFoundCount;
    /** number of good genomes */
    private int goodCount;
    /** total number of genomes */
    private int genomeCount;
    /** role definition map */
    private RoleMap roleDefinitions;
    /** TRUE if completeness data is available, else FALSE */
    private boolean haveCompleteness;

    /**
     * output report format
     */
    public enum Type {
        /** tab-delimited text */
        TEXT,
        /** HTML web pages */
        HTML,
        /** HTML web pages with comparison to a reference genome */
        DEEP,
        /** HTML annotation comparison */
        COMPARE
    }

    /**
     * option flags
     */
    public enum Option {
        /** suppress summary report */
        NOSUMMARY,
        /** suppress detail reports */
        NODETAILS,
        /** PATRIC genome report output file */
        P3REPORT
    }

    /**
     * Construct a reporter for a specified output directory.
     *
     * @param outDir	output directory for reports
     */
    protected EvalReporter() {
        this.outDir = null;
        this.summary = true;
        this.details = true;
        this.version = null;
        this.haveCompleteness = false;
        this.htmlFileName = null;
    }

    /**
     * Construct a report of the specified type with the specified options.
     *
     * @param outDir	output directory
     * @param type		output format
     *
     * @return a reporter object for the specified output
     */
    public static EvalReporter create(Type type) {
        EvalReporter retVal;
        switch (type) {
        case TEXT -> retVal = new EvalTextReporter();
        case HTML -> retVal = new EvalHtmlReporter();
        case DEEP -> retVal = new EvalDeepReporter();
        case COMPARE -> retVal = new EvalCompareReporter();
        default -> throw new RuntimeException("Unsupported output format.");
        }
        return retVal;
    }

    /**
     * Specify an option.  Options must be set prior to initialization.
     *
     * @param option	option to set
     */
    public void setOption(Option option) {
        switch (option) {
        case NOSUMMARY:
            this.summary = false;
            break;
        case NODETAILS:
            this.details = false;
            break;
        case P3REPORT:
            this.htmlFileName = "GenomeReport.html";
            break;
        default:
            break;
        }
    }

    /**
     * Start the reports.
     *
     * @param version	evaluation database version
     * @param roleMap	role definition file
     * @param modeDir	evaluation directory
     *
     * @throws IOException
     */
    public void open(String version, RoleMap roleMap, File modelDir) throws IOException {
        // Store the caller's descriptive information.
        this.version = version;
        this.roleDefinitions = roleMap;
        // Initalize the subclass.
        this.initialize(modelDir);
        // Clear the counters.
        this.cleanCount = 0;
        this.completeCount = 0;
        this.consistentCount = 0;
        this.genomeCount = 0;
        this.goodCount = 0;
        this.goodSeedCount = 0;
        this.understoodCount = 0;
        this.ssuFoundCount = 0;
        // Start the summary report.
        if (this.summary)
            startSummary();
    }

    /**
     * Initialize the reporting data structures.  The version is known at this
     * point.
     *
     * @param modeDir	evaluation directory
     *
     * @throws IOException
     */
    protected abstract void initialize(File modelDir) throws IOException;

    /**
     * Start the summary report.
     *
     * @throws IOException
     */
    protected abstract void startSummary() throws IOException;

    /**
     * Write the detail report for a genome.
     *
     * @param gReport	stats for the genome of interest
     * @param analysis	analysis of the genome to write
     *
     * @throws IOException
     */
    protected abstract void writeDetails(GenomeStats gReport, GenomeAnalysis analysis) throws IOException;

    /**
     * Write the summary data for a genome.
     *
     * @param gReport	stats for the genome of interest
     *
     * @throws IOException
     */
    protected abstract void writeSummary(GenomeStats gReport) throws IOException;

    /**
     * Finish the summary report.
     *
     * @throws IOException
     */
    protected abstract void endSummary() throws IOException;

    /**
     * Output the details for a genome.
     *
     * @param gReport	stats for the genome to write
     * @param analysis	analysis of genome to write
     *
     * @throws IOException
     */
    public void writeGenome(GenomeStats gReport, GenomeAnalysis analysis) throws IOException {
        // Count this genome.
        this.genomeCount++;
        if (gReport.isGood()) this.goodCount++;
        if (gReport.isClean()) this.cleanCount++;
        if (gReport.isComplete()) this.completeCount++;
        if (gReport.isConsistent()) this.consistentCount++;
        if (gReport.isGoodSeed()) this.goodSeedCount++;
        if (gReport.isUnderstood()) this.understoodCount++;
        if (gReport.hasSsuRRna()) this.ssuFoundCount++;
        // Write this genome's data on the summary report.
        if (this.summary)
            writeSummary(gReport);
        // Write this genome's detail report.
        if (this.details)
            writeDetails(gReport, analysis);
    }

    /**
     * Release all resources held by this object.
     */
    protected abstract void finish();

    /**
     * Finish the reports.
     */
    @Override
    public void close() {
        // Close off the summary report.
        try {
            if (this.summary)
                this.endSummary();
        } catch (IOException e) {
            throw new UncheckedIOException("Error in summary report.", e);
        } finally {
            this.finish();
        }
    }

    /**
     * @return the output directory
     */
    protected File getOutDir() {
        return outDir;
    }

    /**
     * @return the version string
     */
    protected String getVersion() {
        return version;
    }

    /**
     * @return the number of clean genomes
     */
    protected int getCleanCount() {
        return cleanCount;
    }

    /**
     * @return the number of consistent genomes
     */
    protected int getConsistentCount() {
        return consistentCount;
    }

    /**
     * @return the number of mostly-complete genomes
     */
    protected int getCompleteCount() {
        return completeCount;
    }

    /**
     * @return the number of genomes with good seed proteins
     */
    protected int getGoodSeedCount() {
        return goodSeedCount;
    }

    /**
     * @return the number of good genomes
     */
    protected int getGoodCount() {
        return goodCount;
    }

    /**
     * @return the total number of genomes
     */
    protected int getGenomeCount() {
        return genomeCount;
    }

    /**
     * @return the number of reasonably understood genomes
     */
    protected int getUnderstoodCount() {
        return understoodCount;
    }

    /**
     * @return the number of SSU rRNAs found
     */
    protected int getSsuFoundCount() {
        return ssuFoundCount;
    }

    /**
     * @return the name of a specified role
     *
     * @param ID of the role whose name is desired
     */
    protected String getRoleName(String role) {
        String retVal = this.roleDefinitions.getName(role);
        if (retVal == null) {
            retVal = "Unknown role " + role;
        }
        return retVal;
    }

    /**
     * @return the full role map
     */
    protected RoleMap getRoleMap() {
        return this.roleDefinitions;
    }

    /**
     * @return the detail level needed in genomes read from PATRIC (can be overridden by subclasses)
     */
    public P3Genome.Details getDetailLevel() {
        return P3Genome.Details.STRUCTURE_ONLY;
    }

    /**
     * @return the name of the HTML file for a genome.
     *
     * @param genomeId	ID of the genome of interest
     */
    protected File htmlFile(String genomeId) {
        String fName;
        if (this.htmlFileName != null)
            fName = this.htmlFileName;
        else
            fName = genomeId + ".html";
        return new File(this.getOutDir(), fName);
    }


   /**
     * Specify the output directory.
     *
     * @param outDir	proposed output directory
     */
    public void setOutDir(File outDir) {
        this.outDir = outDir;

    }

    /**
     * @return TRUE if we have completeness data
     */
    public boolean hasCompleteness() {
        return this.haveCompleteness;
    }

    /**
     * @param haveCompleteness 	TRUE if we have completeness data, else FALSE
     */
    public void setHaveCompleteness(boolean haveCompleteness) {
        this.haveCompleteness = haveCompleteness;
    }

}
