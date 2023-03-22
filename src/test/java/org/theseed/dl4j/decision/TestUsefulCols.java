/**
 *
 */
package org.theseed.dl4j.decision;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.theseed.dl4j.TabbedDataSetReader;

/**
 * @author Bruce Parrello
 *
 */
class TestUsefulCols {

    @Test
    void testUsefulCols() throws IOException {
        List<String> outcomes = Arrays.asList("None", "Low", "High");
        List<String> meta = Arrays.asList("sample_id", "density", "production");
        File partFile = new File("data", "thr2.tbl");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "prod_level",
                outcomes, meta)) {
            reader.setBatchSize(3000);
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            int[] idxes = RandomForest.getUsefulFeatures(readSet);
            int expected = 0;
            for (int i = 0; i < idxes.length; i++) {
                assertThat(Integer.toString(i), idxes[i], equalTo(expected));
                expected++;
                while (expected == 7 || expected == 9)
                    expected++;
            }
        }
        partFile = new File("data", "thr.tbl");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "prod_level",
                outcomes, meta)) {
            reader.setBatchSize(3000);
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            int[] idxes = RandomForest.getUsefulFeatures(readSet);
            int expected = 0;
            for (int i = 0; i < idxes.length; i++) {
                assertThat(Integer.toString(i), idxes[i], equalTo(expected));
                expected++;
            }
        }

    }

}
