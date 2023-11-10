/**
 *
 */
package org.theseed.dl4j.eval.reports;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.theseed.basic.ParseFailureException;
import org.theseed.dl4j.eval.stats.GenomeStats;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;

/**
 * @author Bruce Parrello
 *
 */
class TestRefComputer {

    @Test
    void testP3RefComputer() throws IOException, ParseFailureException {
        // Test a ref computer with a hash.
        RefGenomeComputer computer = new PatricRefGenomeComputer(new File("data", "Eval.Fake"));
        // We need to read in the test genomes and build the genome reports.
        GenomeSource testGenomes = GenomeSource.Type.DIR.create(new File("data", "genomes"));
        GenomeStats[] reports = new GenomeStats[testGenomes.size()];
        int rIdx = 0;
        for (Genome genome : testGenomes) {
            final var report = new GenomeStats(genome);
            // We need to fill in the seed protein.
            genome.getPegs().stream().filter(x ->
                    x.getPegFunction().contentEquals("Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)"))
                    .forEach(x -> report.countSeed(x.getProteinTranslation()));
            reports[rIdx] = report;
            rIdx++;
        }
        // Set up the reference computer.
        computer.setupReferences(reports);
        testReferences(computer, reports);
        // Now try again without a cache.
        computer = new PatricRefGenomeComputer(new File("data", "Eval.P3"));
        computer.setupReferences(reports);
        testReferences(computer, reports);
    }

    /**
     * Test the reference genomes for the specified reference computer.
     *
     * @param computer	reference computer
     * @param reports	GenomeStats array
     */
    protected void testReferences(RefGenomeComputer computer, GenomeStats[] reports) {
        // Test the reference genomes.
        assertThat(reports[0].getId(), equalTo("1001240.4"));
        Genome genome = computer.ref(reports[0].getGenome());
        assertThat(genome.getId(), equalTo("1001240.6"));
        assertThat(genome.getName(), equalTo("Cryobacterium roopkundense strain DSM 21065"));
        assertThat(reports[0].getId(), equalTo("1001240.4"));
        assertThat(reports[1].getId(), equalTo("1001240.6"));
        genome = computer.ref(reports[1].getGenome());
        assertThat(genome.getId(), equalTo("1001240.4"));
        assertThat(genome.getName(), equalTo("Cryobacterium roopkundense RuG17"));
        assertThat(reports[2].getId(), equalTo("1123349.3"));
        genome = computer.ref(reports[2].getGenome());
        assertThat(genome.getId(), equalTo("1123350.4"));
        assertThat(genome.getName(), equalTo("Tepidibacter thalassicus DSM 15285"));
        assertThat(reports[3].getId(), equalTo("1144273.3"));
        genome = computer.ref(reports[3].getGenome());
        assertThat(genome, nullValue());
        assertThat(reports[4].getId(), equalTo("29391.22"));
        genome = computer.ref(reports[4].getGenome());
        assertThat(genome.getId(), equalTo("562982.3"));
        assertThat(genome.getName(), equalTo("Gemella morbillorum M424"));
        assertThat(reports[5].getId(), equalTo("511145.12"));
        genome = computer.ref(reports[5].getGenome());
        assertThat(genome, nullValue());
    }

}
