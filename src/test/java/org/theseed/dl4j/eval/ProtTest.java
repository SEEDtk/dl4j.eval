package org.theseed.dl4j.eval;

import org.junit.jupiter.api.Test;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.junit.jupiter.api.Assertions.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.lang3.StringUtils;
import org.theseed.proteins.RoleMatrix;

public class ProtTest  {


    /**
     * Test the role matrix
     */
    @Test
    public void testRoleMatrix() {
        RoleMatrix testMatrix = new RoleMatrix(2, 5);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(0));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(0));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiHigb"), equalTo(0));
        assertFalse(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        Collection<String> ab = Arrays.asList("AntiHiga", "ToxiHigb");
        Collection<String> ac = Arrays.asList("AntiHiga", "VapcToxiProt");
        Collection<String> abc = Arrays.asList("AntiHiga", "ToxiHigb", "VapcToxiProt");
        Collection<String> abcB = Arrays.asList("AntiHiga", "ToxiHigb", "VapbProt", "VapcToxiProt");
        testMatrix.register("83333.1", ab);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(1));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(0));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        testMatrix.register("83333.1", abc);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(1));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(1));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertTrue(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        testMatrix.register("100226.1", ac);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(2));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(2));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertTrue(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertTrue(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        testMatrix.register("11446.1", ac);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(3));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(3));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertTrue(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertTrue(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertTrue(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        assertThat(testMatrix.completeness(ab, "83333.1"), equalTo(1.0));
        assertThat(testMatrix.completeness(abc, "83333.1"), equalTo(1.0));
        assertThat(testMatrix.completeness(abcB, "83333.1"), closeTo(0.750, 0.001));
        assertThat(testMatrix.completeness(abcB, "100226.1"), closeTo(0.500, 0.001));
        assertThat(testMatrix.completeness(abcB, "11446.1"), closeTo(0.500, 0.001));
        Collection<String> roles = testMatrix.getMarkerRoles(0.75);
        assertThat(testMatrix.completeness(roles, "83333.1"), greaterThanOrEqualTo(0.75));
        assertThat(testMatrix.completeness(roles, "100226.1"), greaterThanOrEqualTo(0.75));
        assertThat(testMatrix.completeness(roles, "11446.1"), greaterThanOrEqualTo(0.75));
    }

    /**
     * Stress test for role matrix
     * @throws IOException
     */
    @Test
    public void testRoleMatrixStress() throws IOException {
        // rickettsia file (100 genomes)
        File inFile = new File("data", "rickettsia.roles.tbl");
        FileReader fileStream = new FileReader(inFile);
        RoleMatrix stressMatrix = new RoleMatrix(100, 100);
        Set<String> genomes = new HashSet<String>();
        try (BufferedReader reader = new BufferedReader(fileStream)) {
            // Skip the header.
            reader.readLine();
            // Loop through the file.
            for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                String[] fields = StringUtils.splitPreserveAllTokens(line, '\t');
                List<String> roles = Arrays.asList(Arrays.copyOfRange(fields, 2, fields.length));
                stressMatrix.register(fields[0], roles);
                genomes.add(fields[0]);
            }
        }
        Set<String> matrixGenomes = stressMatrix.genomes();
        for (String genome : genomes)
            assertTrue(matrixGenomes.contains(genome));
        assertThat(matrixGenomes.size(), equalTo(genomes.size()));
        Collection<String> universals = stressMatrix.getMarkerRoles(0.90);
        for (String genome : genomes) {
            double completeness = stressMatrix.completeness(universals, genome);
            assertThat(completeness, greaterThanOrEqualTo(0.90));
        }
        universals = stressMatrix.getCommonRoles(0.95);
        for (String role : universals) {
            int frequency = (int) (stressMatrix.roleCount(role) / 0.95);
            assertThat(frequency, greaterThanOrEqualTo(genomes.size()));
        }

    }

}
