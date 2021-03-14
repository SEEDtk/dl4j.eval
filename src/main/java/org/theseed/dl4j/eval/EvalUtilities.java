/**
 *
 */
package org.theseed.dl4j.eval;

import org.apache.commons.lang3.StringUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;

/**
 * This class contains utilities for training and testing evaluators.
 *
 * @author Bruce Parrello
 *
 */
public class EvalUtilities {

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
        // Get the features from the input dataset.
        INDArray features = inputSet.getFeatures();
        // Analyze the output labels.  Insure we allow the minimum number of possibilities.
        String[] counts = StringUtils.split(labelString, ",");
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

}
