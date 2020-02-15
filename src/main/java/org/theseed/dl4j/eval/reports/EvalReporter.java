/**
 *
 */
package org.theseed.dl4j.eval.reports;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.eval.GenomeStats;
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
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(EvalReporter.class);
    /** output directory */
    private File outDir;
    /** TRUE to output summary */
    private boolean summary;
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
    /** number of good genomes */
    private int goodCount;
    /** total number of genomes */
    private int genomeCount;
    /** role definition map */
    private RoleMap roleDefinitions;
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
        NODETAILS
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
        case TEXT :
            retVal = new EvalTextReporter();
            break;
        case HTML :
            retVal = new EvalHtmlReporter();
            break;
        case DEEP :
            retVal = new EvalDeepReporter();
            break;
        case COMPARE :
            retVal = new EvalCompareReporter();
            break;
        default :
            throw new RuntimeException("Unsupported output format.");
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
     *
     * @throws IOException
     */
    protected abstract void writeDetails(GenomeStats gReport) throws IOException;

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
     *
     * @throws IOException
     */
    public void writeGenome(GenomeStats gReport) throws IOException {
        // Count this genome.
        this.genomeCount++;
        if (gReport.isGood()) this.goodCount++;
        if (gReport.isClean()) this.cleanCount++;
        if (gReport.isComplete()) this.completeCount++;
        if (gReport.isConsistent()) this.consistentCount++;
        if (gReport.isGoodSeed()) this.goodSeedCount++;
        if (gReport.isUnderstood()) this.understoodCount++;
        // Write this genome's data on the summary report.
        if (this.summary)
            writeSummary(gReport);
        // Write this genome's detail report.
        if (this.details)
            writeDetails(gReport);
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
        return new File(this.getOutDir(), genomeId + ".html");
    }


    /**
     * Perform special setup for the specified batch of evaluated genomes.
     *
     * @param reports	array of genome evaluations
     */
    public abstract void setupGenomes(GenomeStats[] reports);

    /**
     * Specify the output directory.
     *
     * @param outDir	proposed output directory
     */
    public void setOutDir(File outDir) {
        this.outDir = outDir;

    }

}
