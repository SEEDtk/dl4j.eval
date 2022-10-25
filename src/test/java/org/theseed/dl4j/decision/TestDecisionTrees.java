/**
 *
 */
package org.theseed.dl4j.decision;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.TabbedDataSetReader;
import org.theseed.dl4j.decision.RandomForest.Method;
import org.theseed.dl4j.train.ClassPredictError;
import org.theseed.io.TabbedLineReader;

/**
 * Use TESTPAR environment variable to set parallelism.
 *
 * @author Bruce Parrello
 *
 */
public class TestDecisionTrees {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TestDecisionTrees.class);

    public TestDecisionTrees() {
        String width = System.getenv("TESTPAR");
        if (width == null)
            width = "7";
        if (! width.contentEquals("all"))
            System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", width);
    }

    @Test
    public void testDataset() throws IOException {
        log.info("Running testDataset().");
        File partFile = new File("data", "partial_iris.tbl");
        List<String> labels = Arrays.asList("virginica", "versicolor", "setosa");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "species", labels,
                Collections.emptyList())) {
            reader.setBatchSize(200);
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            assertThat(readSet.numExamples(), equalTo(102));
            INDArray row1 = readSet.getFeatures().getRow(0);
            assertThat(row1.getDouble(0), closeTo(6.0, 0.0001));
            assertThat(row1.getDouble(1), closeTo(3.4, 0.0001));
            assertThat(row1.getDouble(2), closeTo(4.5, 0.0001));
            assertThat(row1.getDouble(3), closeTo(1.6, 0.0001));
            assertThat(readSet.numInputs(), equalTo(4));
            assertThat(readSet.numOutcomes(), equalTo(3));
            assertThat(DecisionTree.bestLabel(readSet), equalTo(2));
        }
    }

    @Test
    public void testTree() throws IOException, ClassNotFoundException {
        log.info("Running testTree().");
        File partFile = new File("data", "partial_iris.tbl");
        List<String> outcomes = Arrays.asList("virginica", "versicolor", "setosa");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "species", outcomes,
                Collections.emptyList())) {
            reader.setBatchSize(200);
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            RandomForest.Parms parms = new RandomForest.Parms();
            Iterator<TreeFeatureSelectorFactory> factoryIter = new NormalTreeFeatureSelectorFactory(142857,
                    reader.getWidth(), 4, parms.getNumTrees());
            DecisionTree tree = new DecisionTree(readSet, parms, factoryIter.next());
            // Create a label array for output.
            INDArray predictions = Nd4j.zeros(readSet.numExamples(), readSet.numOutcomes());
            tree.vote(readSet.getFeatures(), predictions);
            // Get the actual labels.
            INDArray expectations = readSet.getLabels();
            // Compare output to actual.
            double good = 0.0;
            double total = 0.0;
            for (int i = 0; i < readSet.numExamples(); i++) {
                int actual = ClassPredictError.computeBest(expectations, i);
                int predicted = ClassPredictError.computeBest(predictions, i);
                total++;
                if (predicted == actual) good++;
            }
            log.info("Good = {}, total = {}.", good, total);
            assertThat(good / total, greaterThan(0.9));
            INDArray impact = tree.computeImpact();
            for (int i = 0; i < impact.columns(); i++) {
                log.info("Impact of {} is {}.", i, impact.getDouble(i));
            }
            File tempFile = new File("data", "tree.ser");
            try (FileOutputStream fileStream = new FileOutputStream(tempFile)) {
                ObjectOutputStream outStream = new ObjectOutputStream(fileStream);
                outStream.writeObject(tree);
            }
            log.info("Tree written to {}.", tempFile);
            DecisionTree tree0;
            try (FileInputStream fileStream = new FileInputStream(tempFile)) {
                ObjectInputStream inStream = new ObjectInputStream(fileStream);
                tree0 = (DecisionTree) inStream.readObject();
            }
            INDArray impact0 = tree0.computeImpact();
            assertThat(impact0, equalTo(impact));
            INDArray predictions0 = Nd4j.zeros(readSet.numExamples(), readSet.numOutcomes());
            tree0.vote(readSet.getFeatures(), predictions0);
            assertThat(predictions0, equalTo(predictions));
        }
    }

    @Test
    public void testLargeDataset() throws IOException {
        log.info("Running testLargeDataset().");
        File partFile = new File("data", "thr.tbl");
        List<String> outcomes = Arrays.asList("None", "Low", "High");
        List<String> meta = Arrays.asList("sample_id", "density", "production");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "prod_level",
                outcomes, meta)) {
            reader.setBatchSize(3000);
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            RandomForest.Parms parms = new RandomForest.Parms(readSet).setNumFeatures(readSet.numInputs());
            Iterator<TreeFeatureSelectorFactory> factoryIter = new NormalTreeFeatureSelectorFactory(142857,
                    reader.getWidth(), parms.getNumFeatures(), parms.getNumTrees());
            DecisionTree tree = new DecisionTree(readSet, parms, factoryIter.next());
            // Create a label array for output.
            INDArray predictions = Nd4j.zeros(readSet.numExamples(), readSet.numOutcomes());
            tree.vote(readSet.getFeatures(), predictions);
            // Get the actual labels.
            INDArray expectations = readSet.getLabels();
            // Compare output to actual.
            double good = 0.0;
            double total = 0.0;
            for (int i = 0; i < readSet.numExamples(); i++) {
                int actual = ClassPredictError.computeBest(expectations, i);
                int predicted = ClassPredictError.computeBest(predictions, i);
                total++;
                if (predicted == actual) good++;
            }
            log.info("BIGTREE: good = {}, total = {}.", good, total);
            assertThat(good / total, greaterThan(0.9));
        }
    }

    @Test
    public void testRandomForest() throws IOException, ClassNotFoundException {
        log.info("Running testRandomForest().");
        File partFile = new File("data", "thr.tbl");
        List<String> outcomes = Arrays.asList("None", "Low", "High");
        List<String> meta = Arrays.asList("sample_id", "density", "production");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "prod_level",
                outcomes, meta)) {
            reader.setBatchSize(3000);
            log.info("Reading dataset.");
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            RandomForest.Parms parms = new RandomForest.Parms(readSet);
            RandomForest.setSeed(142857);
            Iterator<TreeFeatureSelectorFactory> factoryIter = new NormalTreeFeatureSelectorFactory(142857,
                    reader.getWidth(), parms.getNumFeatures(), parms.getNumTrees());
            log.info("Creating random forest. Free memory = {} on {} processors.  Width = {}.",
                    Runtime.getRuntime().freeMemory(), Runtime.getRuntime().availableProcessors(),
                    reader.getWidth());
            RandomForest forest = new RandomForest(readSet, parms, factoryIter);
            // Create a label array for output.
            log.info("Creating predictions.");
            INDArray predictions = forest.predict(readSet.getFeatures());
            // Get the actual labels.
            log.info("Reading expectations.");
            INDArray expectations = readSet.getLabels();
            // Compare output to actual.
            double good = 0.0;
            double total = 0.0;
            for (int i = 0; i < readSet.numExamples(); i++) {
                int actual = ClassPredictError.computeBest(expectations, i);
                int predicted = ClassPredictError.computeBest(predictions, i);
                total++;
                if (predicted == actual) good++;
            }
            log.info("FOREST: good = {}, total = {}.", good, total);
            INDArray impact = forest.computeImpact();
            File tempFile = new File("data", "tree.ser");
            try (FileOutputStream fileStream = new FileOutputStream(tempFile)) {
                ObjectOutputStream outStream = new ObjectOutputStream(fileStream);
                outStream.writeObject(forest);
            }
            log.info("Forest written to {}.", tempFile);
            RandomForest forest0;
            try (FileInputStream fileStream = new FileInputStream(tempFile)) {
                ObjectInputStream inStream = new ObjectInputStream(fileStream);
                forest0 = (RandomForest) inStream.readObject();
            }
            INDArray impact0 = forest0.computeImpact();
            assertThat(impact0, equalTo(impact));
            INDArray predictions0 = forest.predict(readSet.getFeatures());
            assertThat(predictions0, equalTo(predictions));
        }
    }

    @Test
    public void testRandomizerMethods() throws IOException {
        log.info("Running testRandomizerMethods().");
        File partFile = new File("data", "thr.tbl");
        List<String> outcomes = Arrays.asList("None", "Low", "High");
        List<String> meta = Arrays.asList("sample_id", "density", "production");
        try (TabbedDataSetReader reader = new TabbedDataSetReader(partFile, "prod_level",
                outcomes, meta)) {
            reader.setBatchSize(3000);
            DataSet readSet = reader.next();
            INDArray features = readSet.getFeatures();
            readSet.setFeatures(features.reshape(features.size(0), features.size(3)));
            RandomForest.Parms parms = new RandomForest.Parms(readSet).setNumFeatures(14);
            File ratingFile = new File("data", "ratings.tbl");
            List<String> impactCols = TabbedLineReader.readColumn(ratingFile, "1");
            List<Iterator<TreeFeatureSelectorFactory>> finders = Arrays.asList(
                    new NormalTreeFeatureSelectorFactory(142857, reader.getWidth(), parms.getNumFeatures(), parms.getNumTrees()),
                    new RootedTreeFeatureSelectorFactory(142857, reader.getFeatureNames(),
                            impactCols, parms.getNumFeatures(), parms.getNumTrees()));
            for (Method method : Method.values()) {
                for (Iterator<TreeFeatureSelectorFactory> finder : finders) {
                    parms.setMethod(method);
                    log.info("Processing method {} with finder {}.", method, finder);
                    RandomForest forest = new RandomForest(readSet, parms, finder);
                    // Create a label array for output.
                    INDArray predictions = forest.predict(readSet.getFeatures());
                    // Get the actual labels.
                    INDArray expectations = readSet.getLabels();
                    // Compare output to actual.
                    double good = 0.0;
                    double total = 0.0;
                    for (int i = 0; i < readSet.numExamples(); i++) {
                        int actual = ClassPredictError.computeBest(expectations, i);
                        int predicted = ClassPredictError.computeBest(predictions, i);
                        total++;
                        if (predicted == actual) good++;
                    }
                    log.info("FOREST {}[{}]: good = {}, total = {}.", method, finder.toString(), good, total);
                }
            }
        }
    }
}
