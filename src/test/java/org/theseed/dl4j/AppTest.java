package org.theseed.dl4j;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;

/**
 * Simple tests
 */
public class AppTest
{
    // Convert double[] to List<Double>
    public List<Double> collect(double[] parm) {
        return Arrays.stream(parm).boxed().collect(Collectors.toList());
    }

    /**
     * Test the training set reader.
     * @throws IOException
     */
    @Test
    public void testReader() throws IOException {
        File inFile = new File("data", "iris.tbl");
        List<String> labels = Arrays.asList("setosa", "versicolor", "virginica");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(inFile, "species", labels)) {
            reader.setBatchSize(11);
            assertThat("End of file too soon.", reader.hasNext(), equalTo(true));
            assertThat("Wrong input width in file.", reader.getWidth(), equalTo(4));
            DataSet set1 = reader.next();
            assertThat("Wrong number of classes.", set1.numOutcomes(), equalTo(3));
            assertThat("Wrong batch size.", set1.numExamples(), equalTo(11));
            INDArray features = set1.getFeatures();
            INDArray outputs = set1.getLabels();
            assertThat("Wrong sl for ex1.", features.getDouble(0, 1, 1, 0), closeTo(6.0, 0.001));
            assertThat("Wrong sw for ex1.", features.getDouble(0, 1, 1, 1), closeTo(3.4, 0.001));
            assertThat("Wrong pl for ex1.", features.getDouble(0, 1, 1, 2), closeTo(4.5, 0.001));
            assertThat("Wrong pw for ex1.", features.getDouble(0, 1, 1, 3), closeTo(1.6, 0.001));
            assertThat("Wrong label for ex1.", outputs.getDouble(0, 0), equalTo(0.0));
            assertThat("Wrong label aura 1 for ex1.", outputs.getDouble(0, 1), equalTo(1.0));
            assertThat("Wrong label aura 2 for ex1.", outputs.getDouble(0, 2), equalTo(0.0));
            assertThat("Wrong sl for ex9.", features.getDouble(8, 1, 1, 0), closeTo(5.2, 0.001));
            assertThat("Wrong sw for ex9.", features.getDouble(8, 1, 1, 1), closeTo(4.1, 0.001));
            assertThat("Wrong pl for ex9.", features.getDouble(8, 1, 1, 2), closeTo(1.5, 0.001));
            assertThat("Wrong pw for ex9.", features.getDouble(8, 1, 1, 3), closeTo(0.1, 0.001));
            assertThat("Wrong label for ex9.", outputs.getDouble(8, 0), equalTo(1.0));
            assertThat("Wrong label aura 1 for ex9.", outputs.getDouble(8, 1), equalTo(0.0));
            assertThat("Wrong label aura 2 for ex9.", outputs.getDouble(8, 2), equalTo(0.0));
            int count = 1;
            for (DataSet seti : reader) {
                count++;
                if (count <= 13) {
                    assertThat("Wrong batch size for " + count, seti.numExamples(), equalTo(11));
                } else {
                    assertThat("Wrong tail batch size", seti.numExamples(), equalTo(7));
                    assertThat("Wrong label for last example.", seti.getLabels().getDouble(6, 2),
                            equalTo(1.0));
                }
            }
            List<String> names = reader.getFeatureNames();
            assertThat(names, contains("sepal_length", "sepal_width", "petal_length", "petal_width"));
        }
        inFile = new File("data", "samples.tbl");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(inFile, "lcol", labels,
                Arrays.asList("mcol", "mcol2"))) {
            var names = reader.getFeatureNames();
            assertThat(names, contains("col1", "col2", "col3"));
            var set1 = reader.next();
            assertThat("Wrong number of classes.", set1.numOutcomes(), equalTo(3));
            assertThat("Wrong batch size.", set1.numExamples(), equalTo(2));
            var features = set1.getFeatures();
            var outputs = set1.getLabels();
            assertThat("Wrong col1 for ex1.", features.getDouble(0, 1, 1, 0), closeTo(1.0, 0.001));
            assertThat("Wrong col2 for ex1.", features.getDouble(0, 1, 1, 1), closeTo(3.0, 0.001));
            assertThat("Wrong col3 for ex1.", features.getDouble(0, 1, 1, 2), closeTo(6.0, 0.001));
            assertThat("Wrong label for ex1.", outputs.getDouble(0, 0), equalTo(1.0));
            assertThat("Wrong label aura 1 for ex1.", outputs.getDouble(0, 1), equalTo(0.0));
            assertThat("Wrong label aura 2 for ex1.", outputs.getDouble(0, 2), equalTo(0.0));
            assertThat("Wrong col1 for ex2.", features.getDouble(1, 1, 1, 0), closeTo(1.1, 0.001));
            assertThat("Wrong col2 for ex2.", features.getDouble(1, 1, 1, 1), closeTo(3.1, 0.001));
            assertThat("Wrong col3 for ex2.", features.getDouble(1, 1, 1, 2), closeTo(6.1, 0.001));
            assertThat("Wrong label for ex2.", outputs.getDouble(1, 2), equalTo(1.0));
            assertThat("Wrong label aura 1 for ex2.", outputs.getDouble(1, 0), equalTo(0.0));
            assertThat("Wrong label aura 2 for ex2.", outputs.getDouble(1, 1), equalTo(0.0));
        }
        inFile = new File("data", "regression.tbl");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(inFile, null, labels,
                Arrays.asList("mcol", "mcol2"))) {
            reader.setRegressionColumns();
            var names = reader.getFeatureNames();
            assertThat(names, contains("col1", "col2", "col3"));
            var set1 = reader.next();
            assertThat("Wrong number of classes.", set1.numOutcomes(), equalTo(3));
            assertThat("Wrong batch size.", set1.numExamples(), equalTo(2));
            var features = set1.getFeatures();
            var outputs = set1.getLabels();
            assertThat("Wrong col1 for ex1.", features.getDouble(0, 1, 1, 0), closeTo(1.0, 0.001));
            assertThat("Wrong col2 for ex1.", features.getDouble(0, 1, 1, 1), closeTo(3.0, 0.001));
            assertThat("Wrong col3 for ex1.", features.getDouble(0, 1, 1, 2), closeTo(6.0, 0.001));
            assertThat("Wrong setosa for ex1.", outputs.getDouble(0, 0), closeTo(4.0, 0.001));
            assertThat("Wrong versicolor for ex1.", outputs.getDouble(0, 1), closeTo(8.0, 0.001));
            assertThat("Wrong virginica for ex1.", outputs.getDouble(0, 2), closeTo(7.0, 0.001));
            assertThat("Wrong col1 for ex2.", features.getDouble(1, 1, 1, 0), closeTo(1.1, 0.001));
            assertThat("Wrong col2 for ex2.", features.getDouble(1, 1, 1, 1), closeTo(3.1, 0.001));
            assertThat("Wrong col3 for ex2.", features.getDouble(1, 1, 1, 2), closeTo(6.1, 0.001));
            assertThat("Wrong setosa for ex1.", outputs.getDouble(1, 0), closeTo(4.1, 0.001));
            assertThat("Wrong versicolor for ex1.", outputs.getDouble(1, 1), closeTo(8.1, 0.001));
            assertThat("Wrong virginica for ex1.", outputs.getDouble(1, 2), closeTo(7.1, 0.001));
        }
        inFile = new File("data", "predicts.tbl");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(inFile, Arrays.asList("mcol", "mcol2"))) {
            var set1 = reader.next();
            var features = set1.getFeatures();
            List<String> metaData = set1.getExampleMetaData(String.class);
            assertThat("Wrong col1 for ex1.", features.getDouble(0, 1, 1, 0), closeTo(1.0, 0.001));
            assertThat("Wrong col2 for ex1.", features.getDouble(0, 1, 1, 1), closeTo(3.0, 0.001));
            assertThat("Wrong col3 for ex1.", features.getDouble(0, 1, 1, 2), closeTo(6.0, 0.001));
            assertThat("Wrong meta for ex1.", metaData.get(0), equalTo("2.0\t5.0"));
            assertThat("Wrong col1 for ex2.", features.getDouble(1, 1, 1, 0), closeTo(1.1, 0.001));
            assertThat("Wrong col2 for ex2.", features.getDouble(1, 1, 1, 1), closeTo(3.1, 0.001));
            assertThat("Wrong col3 for ex2.", features.getDouble(1, 1, 1, 2), closeTo(6.1, 0.001));
            assertThat("Wrong label for ex2.", metaData.get(1), equalTo("2.1\t5.1"));
        }
    }

    /**
     * Test the channel dataset reader.
     * @throws IOException
     */
    @Test
    public void testChannelReader() throws IOException {
        File channelFile = new File("data", "channels.tbl");
        Map<String, double[]> channelMap = ChannelDataSetReader.readChannelFile(channelFile);
        assertThat(collect(channelMap.get("a")), contains(1.0, 0.0, 0.0, 0.0));
        assertThat(collect(channelMap.get("X")), contains(0.25, 0.25, 0.25, 0.25));
        assertThat(collect(channelMap.get("y")), contains(0.0, 0.5, 0.0, 0.5));
        List<String> labels = TabbedDataSetReader.readLabels(new File("data", "labels.tbl"));
        ChannelDataSetReader reader = new ChannelDataSetReader(new File("data", "training.tbl"),
                "protein", labels, channelMap);
        assertThat(reader.getChannels(), equalTo(4));
        DataSet set1 = reader.next();
        INDArray features = set1.getFeatures();
        INDArray outputs = set1.getLabels();
        assertThat(set1.numExamples(), equalTo(8));
        assertThat(features.getDouble(0, 0, 0, 0), equalTo(0.0));  // c in example 1 position 1
        assertThat(features.getDouble(0, 1, 0, 0), equalTo(1.0));
        assertThat(features.getDouble(0, 2, 0, 0), equalTo(0.0));
        assertThat(features.getDouble(0, 3, 0, 0), equalTo(0.0));
        assertThat(features.getDouble(0, 0, 0, 1), equalTo(0.0));	// g in example 1 position 2
        assertThat(features.getDouble(0, 1, 0, 1), equalTo(0.0));
        assertThat(features.getDouble(0, 2, 0, 1), equalTo(1.0));
        assertThat(features.getDouble(0, 3, 0, 1), equalTo(0.0));
        assertThat(features.getDouble(1, 0, 0, 0), equalTo(0.0));	// t in example 2 position 1
        assertThat(features.getDouble(1, 1, 0, 0), equalTo(0.0));
        assertThat(features.getDouble(1, 2, 0, 0), equalTo(0.0));
        assertThat(features.getDouble(1, 3, 0, 0), equalTo(1.0));
        assertThat(features.getDouble(1, 0, 0, 2), equalTo(1.0));	// a in example 2 position 3
        assertThat(features.getDouble(1, 1, 0, 2), equalTo(0.0));
        assertThat(features.getDouble(1, 2, 0, 2), equalTo(0.0));
        assertThat(features.getDouble(1, 3, 0, 2), equalTo(0.0));
        assertThat(outputs.getDouble(1, 0), equalTo(1.0));	// protein in example 2 is '*'
        assertThat(outputs.getDouble(0, 15), equalTo(1.0)); // protein in example 1 is 'R'
    }

}
