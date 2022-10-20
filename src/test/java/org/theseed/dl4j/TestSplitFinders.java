/**
 *
 */
package org.theseed.dl4j;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.theseed.dl4j.decision.DecisionTree;
import org.theseed.dl4j.decision.SequentialSplitPointFinder;
import org.theseed.dl4j.decision.SplitPointFinder;
import org.theseed.dl4j.decision.Splitter;

/**
 * @author Bruce Parrello
 *
 */
public class TestSplitFinders {

    @Test
    public void test() throws IOException {
        SplitPointFinder finder = new SplitPointFinder.Mean();
        File partFile = new File("src/test/data", "partial_iris.tbl");
        List<String> outcomes = Arrays.asList("virginica", "versicolor", "setosa");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "species", outcomes,
                Collections.emptyList())) {
            reader.setBatchSize(200);
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            double oldEntropy = DecisionTree.entropy(readSet);
            List<DataSet> rows = readSet.asList();
            for (int i = 0; i < 3; i++) {
                Splitter meanSplitter = finder.computeSplit(i, 3, rows, oldEntropy);
                finder = new SequentialSplitPointFinder();
                Splitter bestSplitter = finder.computeSplit(i, 3, rows, oldEntropy);
                assertThat(bestSplitter.compareTo(meanSplitter), lessThanOrEqualTo(0));
            }
        }
    }

}
