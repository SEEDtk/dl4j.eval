/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.KFoldIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.dl4j.TabbedDataSetReader;
import org.theseed.dl4j.decision.RandomForest;
import org.theseed.io.TabbedLineReader;

/**
 * This class contains utilities for training and testing evaluators.
 *
 * @author Bruce Parrello
 *
 */
public class EvalUtilities {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(EvalUtilities.class);


    /**
     * Create a new dataset by removing the specified column and converting it to one-hot labels.
     *
     * @param inputSet		dataset to convert
     * @param colIdx		index of the column to make into labels
     * @param labelString	string of occurrence counts for the column being converted
     *
     * @return a new dataset ready for training or testing
     */
    public static DataSet prepareData(DataSet inputSet, int colIdx, String labelString) {
        String[] counts = StringUtils.split(labelString, ",");
        return prepareData(inputSet, colIdx, counts);
    }

    /**
     * Create a new dataset by removing the specified column and converting it to one-hot labels.
     *
     * @param inputSet		dataset to convert
     * @param colIdx		index of the column to make into labels
     * @param counts		array of occurrence counts for the column being converted
     *
     * @return a new dataset ready for training or testing
     */
    public static DataSet prepareData(DataSet inputSet, int colIdx, String[] counts) {
        // Get the features from the input dataset.
        INDArray features = inputSet.getFeatures();
        // Analyze the output labels.  Insure we allow the minimum number of possibilities.
        int maxCount = counts.length - 1;
        while (counts[maxCount].contentEquals("0")) maxCount--;
        // Build the label array.
        INDArray labels = Nd4j.zeros(features.rows(), maxCount + 1);
        for (int row = 0; row < features.rows(); row++) {
            int output = features.getInt(row, colIdx);
            if (output < maxCount) labels.put(row, output, 1.0);
        }
        // Get a features array without the output column.
        int[] columns = new int[features.columns() - 1];
        for (int i = 0; i < colIdx; i++) columns[i] = i;
        for (int i = colIdx; i < columns.length; i++) columns[i] = i + 1;
        INDArray featuresOut = features.getColumns(columns);
        // Return the features and labels as a dataset.
        DataSet retVal = new DataSet(featuresOut, labels);
        return retVal;
    }

    /**
     * Build a roles.to.use structure from a training/testing file.  The training/testing file
     * has genome IDs in the first column and role IDs as the headers for all the other
     * columns.  For each role ID it has an occurrence counts for the identifier genomes.
     * In some cases, instead of occurrence counts there are 0/1 indicators.
     *
     * A duplicate column ID will cause this method to fail, because it indicates an input error.
     *
     * The output structure is a list mapping each role to an array of occurrence counts.
     *
     * @param dataFile		name of the file containing the training/testing data
     * @param max			maximum permission occurrence-count value
     *
     * @return a map from role IDs to occurrence counts
     *
     * @throws IOException
     */
    public static Map<String, int[]> buildRolesToUse(File dataFile, int max) throws IOException {
        Map<String, int[]> retVal;
        try (TabbedLineReader dataStream = new TabbedLineReader(dataFile)) {
            // Get a hash of all the roles.  For each role we will map to an array of counts.
            // Note that the first label (and therefore the first data column) is an ID, not a role.
            String[] labels = dataStream.getLabels();
            retVal = IntStream.range(1, labels.length).mapToObj(i  -> labels[i])
                    .collect(Collectors.toMap(x -> x, x -> new int[max+1]));
            // Now we total up the columns in each record.
            for (TabbedLineReader.Line line :  dataStream) {
                for (int i = 1; i  <  labels.length; i++) {
                    String role = labels[i];
                    int count = line.getInt(i);
                    if (count > max) count = max;
                    int[] counters = retVal.get(role);
                    counters[count]++;
                }
            }
        }
        return retVal;
    }

    /**
     * Read all the data from a training/testing file.
     *
     * @param inFile		file to read
     * @param metaCols		list of metadata column names
     *
     * @return a dataset of the data in the file; it will not contain labels, only features
     *
     * @throws IOException
     */
    public static DataSet readDataSet(File inFile, List<String> metaCols) throws IOException {
        TabbedDataSetReader dataReader = new TabbedDataSetReader(inFile, metaCols);
        DataSet retVal = dataReader.readAll();
        RandomForest.flattenDataSet(retVal);
        dataReader.close();
        return retVal;
    }

    /**
     * Cross-validate a random forest training/testing set.
     *
     * @param hParms		hyper-parameters
     * @param trainingSet	training/testing set to validate
     * @param kFold			number of folds
     *
     * @return the interquartile range for this dataset
     */
    public static double crossValidate(RandomForest.Parms hParms, DataSet trainingSet, int kFold) {
        KFoldIterator validator = new KFoldIterator(kFold, trainingSet);
        DescriptiveStatistics stats = new DescriptiveStatistics();
        int foldI = 1;
        while (validator.hasNext()) {
            log.info("Testing fold {}.", foldI);
            DataSet trainFold = validator.next();
            RandomForest classFold = new RandomForest(trainFold, hParms);
            double accFold = classFold.getAccuracy(validator.testFold());
            stats.addValue(accFold);
            foldI++;
        }
        // Compute the IQR.
        double iqr = stats.getPercentile(75.0) - stats.getPercentile(25.0);
        return iqr;
    }

}
