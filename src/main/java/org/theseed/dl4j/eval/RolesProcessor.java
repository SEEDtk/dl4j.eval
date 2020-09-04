/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseProcessor;

/**
 * This command produces a report on the predictable roles for the consistency checker.  The "training.tbl" file is
 * used to build a role-count matrix, and it is then run against all the roles in the "Roles" subdirectory.  This
 * information is then used to compute the accuracy for each role.  For each role, the coarse accuracy is the
 * percent of genomes for which the role count was correct.  The overcount is the amount by which the role was
 * over-counted, and the undercount is the amount by which the role was undercounted.
 *
 * The positional parameter is the evaluation directory.  The role report will be written to the standard output.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 * -m	minimum acceptable accuracy
 *
 * @author Bruce Parrello
 *
 */
public class RolesProcessor extends BaseProcessor implements IConsistencyChecker {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RolesProcessor.class);
    /** full list of role IDs */
    private List<String> roles;
    /** map of role IDs to role names */
    private Map<String, String> roleMap;
    /** matrix of role counts-- [genome][role] */
    private int[][] rolesActual;
    /** directory containing the role models */
    private File roleDir;
    /** roles-to-use file */
    private File roleFile;
    /** training file */
    private File trainFile;
    /** number of correct genomes for each role */
    private int[] goodCounts;
    /** total number of genomes */
    private int nGenomes;
    /** total overcount for each role */
    private int[] overCounts;
    /** total undercount for each role */
    private int[] underCounts;

    // COMMAND-LINE OPTIONS

    /** minimum acceptable accuracy */
    @Option(name = "-m", aliases = { "--min" }, metaVar = "95.0", usage = "minimum acceptable accuracy")
    private double minAccuracy;

    /** evaluation directory */
    @Argument(index = 0, metaVar = "evalDir", usage = "evaluation directory")
    private File evalDir;

    @Override
    protected void setDefaults() {
        this.minAccuracy = 93.0;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Verify the evaluation directory.
        if (! this.evalDir.isDirectory())
            throw new FileNotFoundException("Evaluation directory " + this.evalDir + " not found or invalid.");
        // Locate the role directory.
        this.roleDir = new File(this.evalDir, "Roles");
        if (! this.roleDir.isDirectory())
            throw new FileNotFoundException("No Roles subdirectory found in " + this.evalDir + ".");
        // Locate the training file.
        this.trainFile = new File(this.evalDir, "training.tbl");
        if (! this.trainFile.canRead())
            throw new FileNotFoundException("Training file " + this.trainFile + " not found or unreadable.");
        // Locate the roles-to-use file.
        this.roleFile = new File(this.evalDir, "roles.to.use");
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role file " + this.roleFile + " not found or unreadable.");
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // First, we read in roles-to-use.
        log.info("Reading roles.to.use.");
        this.roles = new ArrayList<String>(2500);
        this.roleMap = new HashMap<String, String>(2500);
        try (TabbedLineReader roleStream = new TabbedLineReader(this.roleFile, 4)) {
            for (TabbedLineReader.Line line : roleStream) {
                String roleId = line.get(0);
                this.roles.add(roleId);
                this.roleMap.put(roleId, line.get(3));
            }
        }
        log.info("{} roles found in {}.", this.roles.size(), this.roleFile);
        // Now we build the count matrix from the training table.
        try (TabbedLineReader trainStream = new TabbedLineReader(this.trainFile)) {
            log.info("Reading training.tbl");
            List<int[]> counts = new ArrayList<int[]>(1500);
            for (TabbedLineReader.Line line : trainStream) {
                int[] countRow = new int[this.roles.size()];
                for (int i = 0; i < countRow.length; i++)
                    countRow[i] = line.getInt(i+1);
                counts.add(countRow);
            }
            this.nGenomes = counts.size();
            this.rolesActual = new int[this.nGenomes][];
            this.rolesActual = counts.toArray(this.rolesActual);
            log.info("{} genomes in training table.", this.nGenomes);
        }
        // Initialize the result arrays.
        this.goodCounts = new int[this.roles.size()];
        this.overCounts = new int[this.roles.size()];
        this.underCounts = new int[this.roles.size()];
        // Now run the consistency checker.
        boolean[] rolesUsed = Evaluator.runConsistency(this, this.roles, this.rolesActual);
        // Compute the conversion factor from a good-count to percent accuracy.
        double factor = 100.0 / this.nGenomes;
        // Produce the output.
        int roleCount = 0;
        log.info("Pruning bad roles and writing report.");
        System.out.println("role_id\taccuracy\tover_count\tunder_count\tdescription");
        for (int i = 0; i < this.roles.size(); i++) {
            if (rolesUsed[i]) {
                // Here we have a predictable role.
                String roleId = this.roles.get(i);
                String roleName = this.roleMap.get(roleId);
                double accuracy = this.goodCounts[i] * factor;
                if (accuracy < this.minAccuracy) {
                    log.info("Deleting role {}: {}", roleId, roleName);
                    FileUtils.forceDelete(new File(this.roleDir, roleId + ".ser"));
                } else {
                    System.out.format("%s\t%4.2f\t%d\t%d\t%s%n", roleId, accuracy, this.overCounts[i],
                            this.underCounts[i], roleName);
                    roleCount++;
                }
            }
        }
        log.info("{} predictable roles output.",roleCount);
    }

    @Override
    public File getRoleDir() {
        return this.roleDir;
    }

    @Override
    public void storeActual(int iRole, String role, int iGenome, int count) {
        int oldCount = this.rolesActual[iGenome][iRole];
        if (count == oldCount)
            this.goodCounts[iRole]++;
        else if (count < oldCount)
            this.underCounts[iRole] += oldCount - count;
        else
            this.overCounts[iRole] += count - oldCount;
    }

}
