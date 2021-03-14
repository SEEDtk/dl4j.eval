/**
 *
 */
package org.theseed.dl4j.eval;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.theseed.dl4j.TabbedDataSetReader;
import org.theseed.dl4j.decision.RandomForest;

/**
 * @author Bruce Parrello
 *
 */
public class TestUtils {

    @Test
    public void testSlicing() throws IOException {
        File testFile = new File("data", "testing.tbl");
        List<String> metaCols = Arrays.asList(new String[] { "genome" });
        TabbedDataSetReader reader = new TabbedDataSetReader(testFile, metaCols);
        reader.setBatchSize(200);
        DataSet inputSet = reader.next();
        RandomForest.flattenDataSet(inputSet);
        reader.close();
        DataSet testingSet = EvalUtilities.prepareData(inputSet, 453, "287,648,2,0,0,0");
        INDArray inFeatures = inputSet.getFeatures();
        INDArray outFeatures = testingSet.getFeatures();
        INDArray outLabels = testingSet.getLabels();
        assertThat(outLabels.columns(), equalTo(3));
        assertThat(outLabels.rows(), equalTo(inFeatures.rows()));
        assertThat(outFeatures.rows(), equalTo(inFeatures.rows()));
        for (int r = 0; r < inFeatures.rows(); r++) {
            for (int i = 0; i < outFeatures.columns(); i++) {
                int i1 = (i >= 453 ? i+1 : i);
                assertThat(String.format("Row %d col %d", r, i), outFeatures.getDouble(r, i), equalTo(inFeatures.getDouble(r, i1)));
            }
            int label = inFeatures.getInt(r, 453);
            for (int i = 0; i < 3; i++) {
                double expected = (i == label ? 1.0 : 0.0);
                assertThat(String.format("Row %d label %d", r, i), outLabels.getDouble(r, i), equalTo(expected));
            }
        }
    }

}
