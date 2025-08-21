/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.dl4j.decision.RandomForest;
import org.theseed.dl4j.eval.reports.RoleInfluence;
import org.theseed.io.TabbedLineReader;

/**
 * This command looks at a directory built by BuildProcessor and determines which input roles had the most impact on the outputs.
 * The coefficients from the input layer are summed and compared.  The output report is on the standard output.
 *
 * The positional parameter is the evaluation directory created by BuildProcessor.  We use "roles.to.use" to get the input
 * columns, then look in the "Roles" subdirectory for the models created.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 *
 * @author Bruce Parrello
 *
 */
public class AnalyzeProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(AnalyzeProcessor.class);
    /** list of roles */
    private RoleInfluence[] roles;
    /** input role file */
    private File rolesToUse;
    /** role model directory */
    private File roleDir;

    // COMMAND-LINE OPTIONS

    /** input evaluation directory */
    @Argument(index = 0, metaVar = "evalDir", usage = "input evaluation directory")
    private File evalDir;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected void validateParms() throws IOException {
        // Insure we have all the files we need.
        if (! this.evalDir.isDirectory())
            throw new FileNotFoundException("Evaluator directory" + this.evalDir + " not found or invalid.");
        this.rolesToUse = new File(this.evalDir, "roles.to.use");
        if (! this.rolesToUse.canRead())
            throw new FileNotFoundException("Role file " + this.rolesToUse + " not found or unreadable.");
        this.roleDir = new File(this.evalDir, "Roles");
        if (! this.roleDir.isDirectory())
            throw new FileNotFoundException("Model directory " + this.roleDir + " not found or invalid.");
    }

    @Override
    protected void runCommand() throws Exception {
        // Read in the role IDs and names.
        this.roles = this.readRoles();
        log.info("{} roles read from {}.", this.roles.length, this.rolesToUse);
        // Now we loop through the models, analyzing them.
        for (int i = 0; i < this.roles.length; i++) {
            RoleInfluence role = this.roles[i];
            log.info("Processing role {}: {}.", role.getId(), role.getName());
            File modelFile = new File(this.roleDir, role.getId() + ".ser");
            if (! modelFile.exists())
                log.warn("No model file for {}.", role.getId());
            else {
                RandomForest model = RandomForest.load(modelFile);
                // Get the impact.
                INDArray weights = model.computeImpact();
                log.info("Weight matrix has {} elements.", weights.length());
                for (int j = 0; j < weights.length(); j++) {
                    // Determine the role to which this weight belongs and add it it.
                    int jReal = (j >= i ? j + 1 : j);
                    this.roles[jReal].increment(Math.abs(weights.getDouble(j)));
                }
            }
        }
        // Sort and output the results.
        log.info("Sorting and writing results.");
        List<RoleInfluence> sortedRoles = Arrays.stream(this.roles).sorted().collect(Collectors.toList());
        System.out.println("role_id\trole_name\tinfluence");
        for (RoleInfluence role : sortedRoles)
            System.out.format("%s\t%s\t%8.4f%n", role.getId(), role.getName(), role.getRating());
    }

    private RoleInfluence[] readRoles() throws IOException {
        List<RoleInfluence> roleList = new ArrayList<RoleInfluence>(2500);
        try (TabbedLineReader roleStream = new TabbedLineReader(this.rolesToUse, 4)) {
            for (TabbedLineReader.Line line : roleStream) {
                RoleInfluence roleData = new RoleInfluence(line.get(0), line.get(3));
                roleList.add(roleData);
            }
        }
        RoleInfluence[] retVal = new RoleInfluence[roleList.size()];
        return roleList.toArray(retVal);
    }
}
