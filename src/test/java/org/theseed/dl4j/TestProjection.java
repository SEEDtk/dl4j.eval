/**
 *
 */
package org.theseed.dl4j;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.junit.jupiter.api.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dimensionalityreduction.PCA;
import org.nd4j.linalg.ops.transforms.Transforms;
import org.theseed.dl4j.decision.RandomForest;

/**
 * @author Bruce Parrello
 *
 */
class TestProjection {

    @Test
    void testPCA() throws IOException {
        File inFile = new File("data", "continuous.tbl");
        try (TabbedDataSetReader inStream = new TabbedDataSetReader(inFile,
                Arrays.asList("sample_id", "density", "production"))) {
            DataSet inData = inStream.readAll();
            RandomForest.flattenDataSet(inData);
            INDArray inFeatures = inData.getFeatures();
            INDArray workMatrix = inFeatures.dup();
            INDArray projector = PCA.pca_factor(workMatrix, 2, false);
            INDArray projected = inFeatures.mmul(projector);
            assertThat(projected.rows(), equalTo(inFeatures.rows()));
            assertThat(projected.columns(), equalTo(2));
            File outFile = new File("data", "distances.ser");
            int rN = projected.rows();
            int size = rN * (rN + 1) / 2;
            int used = 0;
            double[] oldD = new double[size];
            double[] newD = new double[size];
            try (PrintWriter writer = new PrintWriter(outFile)) {
                writer.println("r1\tr2\told\tnew");
                for (int r1 = 0; r1 < rN; r1++) {
                    for (int r2 = r1 + 1; r2 < rN; r2++) {
                        double oldDist = Transforms.euclideanDistance(inFeatures.slice(r1, 0), inFeatures.slice(r2, 0));
                        double newDist = Transforms.euclideanDistance(projected.slice(r1, 0), projected.slice(r2, 0));
                        writer.format("%d\t%d\t%6.4f\t%6.4f%n", r1, r2, oldDist, newDist);
                        oldD[used] = oldDist;
                        newD[used] = newDist;
                        used++;
                    }
                }
            }
            PearsonsCorrelation pCorr = new PearsonsCorrelation();
            double coeff = pCorr.correlation(oldD, newD);
            assertThat(coeff, greaterThan(0.5));
        }
    }

}
