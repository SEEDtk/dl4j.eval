/**
 *
 */
package org.theseed.dl4j.eval;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.dl4j.train.ClassPredictError;
import org.theseed.io.TabbedLineReader;

/**
 * @author Bruce Parrello
 *
 */
class JackKnifeTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(JackKnifeTest.class);


    @Test
    void test() throws IOException, ParseFailureException {
        RoleTrainProcessor processor = new RoleTrainProcessor();
        File inDir = new File("data", "100226.15");
        // Insure the missing-roles file is not leftover from a previous test.
        File missingRolesFile = new File(inDir, "missingRoles.tbl");
        if (missingRolesFile.exists())
        	FileUtils.forceDelete(missingRolesFile);
        processor.setInDir(inDir);
        processor.setDefaults();
        processor.validateParms();
        processor.setupDataTable();
        DataSet master = processor.getMasterTable();
        File inFile = new File(inDir, "data.tbl");
        var roleNames = processor.getRoleNames();
        var uselessRoles = processor.getUselessRoles();
        // Verify each column.
        for (int c = 0; c < master.numInputs(); c++) {
            if (uselessRoles.get(c)) {
                log.info("Column {} is useless.", c);
            } else {
                try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
                    int inCol = c + 1;
                    String roleId = inStream.getLabels()[inCol];
                    assertThat(roleNames.get(c), equalTo(roleId));
                    log.info("Processing role {}.", roleId);
                    INDArray features = master.getFeatures();
                    int r = 0;
                    var gMap = new HashMap<String, Integer>(700);
                    for (TabbedLineReader.Line line : inStream) {
                        var genome = line.get(0);
                        var count = line.getDouble(inCol);
                        gMap.put(genome, (int) count);
                        assertThat(String.format("%s[%d,%d]", roleId, r, c), features.getDouble(r, c), equalTo(count));
                        r++;
                    }
                    // Check that we can jacknife this column correctly.
                    processor.buildDataSets(c);
                    log.info("Checking testing set: {} genomes.", processor.getTestGenomes().length);
                    this.checkActuals("Testing", c, gMap, processor.getTestGenomes(), processor.getTestingSet());
                    log.info("Checking training set: {} genomes.", processor.getTrainGenomes().length);
                    this.checkActuals("Training", c, gMap, processor.getTrainGenomes(), processor.getTrainingSet());
                }
            }
        }

    }


    private void checkActuals(String type, int col, HashMap<String, Integer> gMap,
            String[] genomes, DataSet set) {
        INDArray actuals = set.getLabels();
        for (int r = 0; r < set.numExamples(); r++) {
            int actualCount = ClassPredictError.computeBest(actuals, r);
            int inputCount = gMap.get(genomes[r]);
            assertThat(String.format("%s genome %s, row %d, col %d", type, genomes[r], r, col),
                    actualCount, equalTo(inputCount));
        }

    }

}
