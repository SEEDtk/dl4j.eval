/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command produces a list of the singly-occurring roles for each genome in a set of input genomes.
 * Only the roles found in a specified role map are used.
 *
 * The report is written to the standard output.  Each output line consists of a genome ID, a feature ID, the role ID,
 * and the feature's functional assignment.
 *
 * The positional parameters are the name of the role map file and the name of the input file or directory.  The input
 * is a genome source-- either a master directory, a GTO directory (the default), or a file of PATRIC genome IDs.  The
 * role map is a tab-delimited file with role IDs in the first column and role descriptions in the third column.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output for report (if not STDOUT)
 *
 * --source		specify the type of input (PATRIC, DIR, or MASTER); default DIR
 *
 * @author Bruce Parrello
 */
public class UniRoleProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(UniRoleProcessor.class);
    /** role definitions */
    private RoleMap roleMap;
    /** genome source */
    private GenomeSource genomes;
    /** output stream */
    private OutputStream output;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome input (master genome directory, GTO directory, patric ID file)")
    private GenomeSource.Type inType;

    /** output file (if not STDOUT) */
    @Option(name = "--output", aliases = { "-o" }, usage = "output file (if not STDOUT)")
    private File outFile;

    /** role definition file */
    @Argument(index = 0, metaVar = "roles.in.subsystems", usage = "role definition file", required = true)
    private File roleFile;

    /** input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "input genome file or directory", required = true)
    private File inDir;


    @Override
    protected void setDefaults() {
        this.inType = GenomeSource.Type.DIR;
        this.outFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the role file.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " not found or unreadable.");
        this.roleMap = RoleMap.load(this.roleFile);
        log.info("{} role definitions loaded from {}.", this.roleFile);
        // Set up the output.
        if (this.outFile == null) {
            log.info("Report will be written to standard output.");
            this.output = System.out;
        } else {
            log.info("Report will be written to {}.", this.outFile);
            this.output = new FileOutputStream(this.outFile);
        }
        // Verify the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Genome source " + this.inDir + " not found.");
        this.genomes = this.inType.create(this.inDir);
        log.info("{} genomes found in {}.", this.genomes.size(), this.inDir);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the role list map.
        Map<String, List<Feature>> roleLists = new TreeMap<String, List<Feature>>();
        try (PrintWriter writer = new PrintWriter(this.output)) {
            // Write the output header.
            writer.println("genome_id\tfid\trole_id\tfunction");
            // Loop through the genomes.
            int count = 0;
            for (Genome genome : genomes) {
                count++;
                log.info("Processing genome {} of {}: {}.", count, genomes.size(), genome);
                // Sort the genome features into lists by role.
                roleLists.clear();
                for (Feature feat : genome.getPegs()) {
                    for (Role role : feat.getUsefulRoles(this.roleMap)) {
                        String roleId = role.getId();
                        List<Feature> flist = roleLists.computeIfAbsent(roleId, x -> new ArrayList<Feature>());
                        flist.add(feat);
                    }
                }
                // Now output the roles with only one feature.
                int rCount = 0;
                for (Map.Entry<String, List<Feature>> fEntry : roleLists.entrySet()) {
                    List<Feature> feats = fEntry.getValue();
                    if (feats.size() == 1) {
                        Feature feat = feats.get(0);
                        writer.format("%s\t%s\t%s\t%s%n", genome.getId(), feat.getId(),
                                fEntry.getKey(), feat.getPegFunction());
                        rCount++;
                    }
                }
                log.info("{} singleton roles found for {}.", rCount, genome);
            }
        } finally {
            // Insure we close the output stream.  This really only matters if we are
            // being called as a service by a client program, since in that case we
            // want the file closed immediately, rather than waiting for the process to
            // finish.
            if (this.outFile != null)
                this.output.close();
        }
    }

}
