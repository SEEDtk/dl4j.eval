package org.theseed.dl4j.eval;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.kmers.ProteinKmers;

/**
 * Unit test for evaluation structures
 */
public class AppTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * genome stats test
     * @throws IOException
     * @throws NumberFormatException
     */
    public void testGenomeStats() throws NumberFormatException, IOException {
        GenomeStats testGenome = new GenomeStats("100226.11", "Bacteria", "Streptomyces coelicolor clonal population");
        testGenome.setGroup("root");
        testGenome.completeRole("RoleG", 0);
        testGenome.completeRole("RoleH", 0);
        testGenome.completeRole("RoleA", 1);
        testGenome.completeRole("RoleB", 1);
        testGenome.completeRole("RoleC", 1);
        testGenome.completeRole("RoleD", 1);
        testGenome.completeRole("RoleE", 2);
        testGenome.completeRole("RoleF", 2);
        testGenome.consistentRole("RoleA", 1, 1);
        testGenome.consistentRole("RoleB", 1, 1);
        testGenome.consistentRole("RoleC", 1, 1);
        testGenome.consistentRole("RoleD", 1, 1);
        testGenome.consistentRole("RoleE", 1, 2);
        testGenome.consistentRole("RoleF", 1, 2);
        testGenome.consistentRole("RoleG", 1, 0);
        testGenome.consistentRole("RoleH", 2, 0);
        testGenome.consistentRole("RoleI", 0, 1);
        testGenome.consistentRole("RoleJ", 0, 1);
        assertThat(testGenome.getId(), equalTo("100226.11"));
        assertThat(testGenome.getName(), equalTo("Streptomyces coelicolor clonal population"));
        assertThat(testGenome.getGroup(), equalTo("root"));
        assertThat(testGenome.getCoarsePercent(), equalTo(60.0));
        assertThat(testGenome.getFinePercent(), equalTo(40.0));
        assertThat(testGenome.getCompletePercent(), equalTo(75.0));
        assertThat(testGenome.getContaminationPercent(), equalTo(20.0));
        // Make sure we can catch fractional percents.
        testGenome.completeRole("RoleK", 1);
        assertThat(testGenome.getCompletePercent(), closeTo(77.78, 0.005));
        assertThat(testGenome.getContaminationPercent(), closeTo(18.18, 0.005));
        testGenome.consistentRole("RoleL", 2, 2);
        assertThat(testGenome.getFinePercent(), closeTo(45.45, 0.005));
        assertThat(testGenome.getCoarsePercent(), closeTo(63.64, 0.005));
        // Get the list of problematic roles.
        Collection<String> pprs = testGenome.getProblematicRoles();
        assertThat(pprs, contains("RoleE", "RoleF", "RoleG", "RoleH", "RoleI", "RoleJ"));
        testGenome.completeRole("RoleM", 2);
        GenomeStats.ProblematicRole ppr = testGenome.getReport("RoleM");
        assertTrue(ppr.isUniversal());
        assertThat(ppr.getPredicted(), equalTo(1));
        assertThat(ppr.getActual(), equalTo(2));
        ppr = testGenome.getReport("RoleH");
        assertFalse(ppr.isUniversal());
        assertThat(ppr.getPredicted(), equalTo(2));
        assertThat(ppr.getActual(), equalTo(0));
        // Test seed proteins
        assertFalse(testGenome.isGoodSeed());
        String p1 = StringUtils.repeat('X', 250);
        testGenome.countSeed(p1);
        assertTrue(testGenome.isGoodSeed());
        assertThat(testGenome.getSeed(), equalTo(p1));
        testGenome.countSeed("MNNNNN");
        assertFalse(testGenome.isGoodSeed());
        assertThat(testGenome.getSeed(), equalTo(p1));
        String p2 = StringUtils.repeat('X', 350);
        testGenome.countSeed(p2);
        assertFalse(testGenome.isGoodSeed());
        assertThat(testGenome.getSeed(), equalTo(p2));
        testGenome = new GenomeStats("83333.1", "Bacteria", "Escherichia coli K12");
        testGenome.countSeed(p2);
        assertTrue(testGenome.isGoodSeed());
        testGenome = new GenomeStats("83333.1", "Bacteria", "Escherichia coli K12");
        testGenome.countSeed(StringUtils.repeat('X', 500));
        assertFalse(testGenome.isGoodSeed());
        testGenome = new GenomeStats("1100226.1", "Archaea", "unclassified archaeon");
        assertFalse(testGenome.isGoodSeed());
        testGenome.countSeed(p1);
        assertFalse(testGenome.isGoodSeed());
        assertThat(testGenome.getSeed(), equalTo(p1));
        testGenome.countSeed("MNNNNN");
        assertFalse(testGenome.isGoodSeed());
        assertThat(testGenome.getSeed(), equalTo(p1));
        testGenome.countSeed(p2);
        assertFalse(testGenome.isGoodSeed());
        assertThat(testGenome.getSeed(), equalTo(p2));
        testGenome = new GenomeStats("83333.1", "Archaea", "Archaeochia coli");
        testGenome.countSeed(p2);
        assertTrue(testGenome.isGoodSeed());
        testGenome = new GenomeStats("83333.1", "Archaea", "Archaeochia coli");
        testGenome.countSeed(StringUtils.repeat('X', 500));
        assertTrue(testGenome.isGoodSeed());
        // Test peg counting.
        Genome myGenome = new Genome(new File("src/test", "test.gto"));
        testGenome = new GenomeStats(myGenome.getId(), myGenome.getDomain(), myGenome.getName());
        for (Feature peg : myGenome.getPegs())
            testGenome.countPeg(peg);
        testGenome.computeMetrics(myGenome);
        assertThat(testGenome.getContigCount(), equalTo(157));
        assertThat(testGenome.getDnaSize(), equalTo(1111866));
        assertThat(testGenome.getCdsPercent(), closeTo(115.48, 0.005));
        assertThat(testGenome.getHypotheticalPercent(), closeTo(50.16, 0.005));
        assertThat(testGenome.getL50(), equalTo(44));
        assertThat(testGenome.getN50(), equalTo(8620));
        assertThat(testGenome.getPlfamPercent(), closeTo(78.12, 0.005));
        assertThat(testGenome.getPegCount(), equalTo(1284));
        assertThat(testGenome.getHypoCount(), equalTo(644));
        assertThat(testGenome.getPlfamCount(), equalTo(1003));
    }

    /**
     * test universal role object
     * @throws IOException
     */
    public void testUniversalRoles() throws IOException {
        UniversalRoles uniRoles = new UniversalRoles(100, "R100 (Acidaminococcus fermentans DSM 20731)",
                "MSDKLQELREKIQKDLSQVKSVEDLKNIRVQYLGKKGALTSILRSLGDVAAEERPKIGKMVNEVRAKMEQRINEQMKLLEAHQMEEKLASETLDFTLPGRKPALGHLHPVTQTLRDIKKVFMRMGFEVVEGPEIETDYFNFEALNLPKDHPARDMQDTFYITDDILLRTQTSGVQARTMQSREPNTPIRMICPGTVYRNDYDATHSPMFHQVEGLVVDKDISLADLKGTLELFCKEMFGDSVKIRLRPSFFPFTEPSCEVDISCVMCGGKGCRVCKNSGWLEILGAGMVHPNVLRMSGYDPDKMKGFAFGMGVERIAMLRYGIDDLRLFFENDLRFIRQF");
        uniRoles.addRole("RoleA");
        uniRoles.addRole("RoleB");
        uniRoles.addRole("RoleC");
        assertThat(uniRoles.size(), equalTo(3));
        assertTrue(uniRoles.contains("RoleA"));
        assertTrue(uniRoles.contains("RoleB"));
        assertTrue(uniRoles.contains("RoleC"));
        assertFalse(uniRoles.contains("RoleD"));
        Collection<String> markerRoles = uniRoles.getRoles();
        assertThat(markerRoles, containsInAnyOrder("RoleA", "RoleB", "RoleC"));
        String s1773_3521 = "MGDPPLESIVSMLSPEALTTAVDAAQQAIALADTLDVLARVKTEHLGDRSPLALARQALAVLPKEQRAEAGKRVNAARNAAQRSYDERLATLRAERDAAVLVAEGIDVTLPSTRVPAGARHPIIMLAEHVADTFIAMGWELAEGPEVETEQFNFDALNFPADHPARGEQDTFYIAPEDSRQLLRTHTSPVQIRTLLARELPVYIISIGRTFRTDELDATHTPIFHQVEGLAVDRGLSMAHLRGTLDAFARAEFGPSARTRIRPHFFPFTEPSAEVDVWFANKIGGADWVEWGGCGMVHPNVLRATGIDPDLYSGFAFGMGLERTLQFRNGIPDMRDMVEGDVRFSLPFGVGA";
        String s563191_3 = "MYDKLEELREKIRKDLGEVKCVDDLKNIRVQYLGKKGALTEILRGLGSVAAEERPKVGKMVNEVRSKVEARIADQMKLLEARQLEEKMASEKIDVTLPGRKAAEGHLHPVTLTLREIKKVFMRMGFEVAEGPEIENDYFNFEALNLPKDHPARDMQDTFYLTDEFLMRTQTSPVQARTMQSRTPNSPIRMICPGTVYRNDYDATHSPMFHQVEGLVIDKNISLADLKGTLELFCKEMFGSSVKIRLRPSFFPFTEPSAEVDISCVICGGKGCRVCKNSGWLEILGAGMVHPNVLRMSGYDPEKVSGFAFGMGVERIAMLRYGIDDLRLFFENDLRFIRQFK";
        String s1262687_3 = "MSDKLQELREKIQKDLSEVKSVEDLKNIRVQYLGKKGALTSILRSLGDVAAEERPKIGKMVNEVRAKVEQRINEQMKLLEAHQMEEKLASETLDFTLPGRKPALGHLHPVTQTLRDIKKVFMRMGFEVVEGPEIETDYFNFEALNLPKDHPARDMQDTFYITDDILLRTQTSGVQARTMQSREPNSPIRMICPGTVYRNDYDATHSPMFHQVEGLVVDKDISLGDLKGTLELFCKEMFGDSVKIRLRPSFFPFTEPSCEVDISCVMCGGKGCRVCKNSGWLEILGAGMVHPNVLRMSGYDPDKMKGFAFGMGVERIAMLRYGIDDLRLFFENDLRFIRQF";
        String s1104577_5 = "MIDKLIEKLHPLERKVLPFLKLKNVKEIIEKSGMQEIEVMRALQWLENKDVLKINIEIKENVYLGENGEKYLKEGLPERRFLECLDKELSLNEIKKNARLDKDEISACLGILKREKAIEFVYDKVKRTDKWKEVLSKIKKDEDFLKSLPSESKHLGSRLDEFRKRKELIEIKIEKIRDIKLNELGEKLIKTEIKDNFIDSVDSKILKNKEWKNKRFRAYDIKINVPKIYNGRKHFVNESIEYARQIWLEMGFKEMSGPLVQTSFWNFDALFTAQDHPVREMQDTFFIKNPAKGKLPKELVDKIRKVHENGWTTNSKGWGGKWNEEEAKKNVLRTHTTVLSARTIAALKKEDLPAKFFSVGRCFRNETLDWSHLFEFNQTEGIVVDPNANFRHLLGYLKEFFKKMGFENARFRPAYFPYCEPNVEIEVYHPVHKKWIELGGAGIFRPEVVKPLLGFECPVLAWGPGFDRMIMDYYRINDIRELYKNDLKQIREARVWLK";
        ProteinKmers p1773_3521 = new ProteinKmers(s1773_3521);
        ProteinKmers p563191_3 = new ProteinKmers(s563191_3);
        ProteinKmers p1262687_3 = new ProteinKmers(s1262687_3);
        ProteinKmers p1104577_5 = new ProteinKmers(s1104577_5);
        assertThat(uniRoles.inGroup(p1773_3521), equalTo(0));
        assertThat(uniRoles.inGroup(p563191_3), equalTo(137));
        assertThat(uniRoles.inGroup(p1262687_3), equalTo(301));
        assertThat(uniRoles.inGroup(p1104577_5), equalTo(0));
        uniRoles = new UniversalRoles(0, "root", "");
        uniRoles.addRole("RoleD");
        uniRoles.addRole("RoleE");
        assertThat(uniRoles.size(), equalTo(2));
        assertThat(uniRoles.inGroup(p1773_3521), equalTo(1));
        assertThat(uniRoles.inGroup(p563191_3), equalTo(1));
        assertThat(uniRoles.inGroup(p1262687_3), equalTo(1));
        assertThat(uniRoles.inGroup(p1104577_5), equalTo(1));
        File compFile = new File("src/test", "comp.tbl");
        List<UniversalRoles> compList = UniversalRoles.Load(compFile);
        assertThat(compList.size(), equalTo(4));
        uniRoles = compList.get(0);
        assertThat(uniRoles.inGroup(p1773_3521), equalTo(156));
        assertThat(uniRoles.inGroup(p563191_3), equalTo(0));
        assertThat(uniRoles.inGroup(p1262687_3), equalTo(0));
        assertThat(uniRoles.inGroup(p1104577_5), equalTo(0));
        assertThat(uniRoles.getName(), equalTo("R100 (Mycobacterium intracellulare MOTT-64)"));
        assertTrue(uniRoles.contains("PhenTrnaSyntAlph"));
        assertTrue(uniRoles.contains("PhosPhosSixa"));
        assertTrue(uniRoles.contains("ImidGlycPhosSynt4"));
        assertTrue(uniRoles.contains("InosUridPrefNucl"));
        assertFalse(uniRoles.contains("MbthLikeNrpsChap4"));
        assertThat(uniRoles.size(), equalTo(411));
        uniRoles = compList.get(1);
        assertThat(uniRoles.inGroup(p1773_3521), equalTo(ProteinKmers.INFINITY));
        assertThat(uniRoles.inGroup(p563191_3), equalTo(0));
        assertThat(uniRoles.inGroup(p1262687_3), equalTo(0));
        assertThat(uniRoles.inGroup(p1104577_5), equalTo(0));
        assertThat(uniRoles.getName(), equalTo("R200 (Mycobacterium bovis AF2122/97)"));
        assertTrue(uniRoles.contains("PhenTrnaSyntAlph"));
        assertTrue(uniRoles.contains("PhosPhosSixa"));
        assertTrue(uniRoles.contains("ImidGlycPhosSynt4"));
        assertTrue(uniRoles.contains("InosUridPrefNucl"));
        assertTrue(uniRoles.contains("MbthLikeNrpsChap4"));
        assertThat(uniRoles.size(), equalTo(1388));
        uniRoles = compList.get(2);
        assertThat(uniRoles.inGroup(p1773_3521), equalTo(0));
        assertThat(uniRoles.inGroup(p563191_3), equalTo(96));
        assertThat(uniRoles.inGroup(p1262687_3), equalTo(75));
        assertThat(uniRoles.inGroup(p1104577_5), equalTo(0));
        assertThat(uniRoles.getName(), equalTo("R50 (Phascolarctobacterium sp. CAG:266)"));
        assertTrue(uniRoles.contains("3DeoxDMannOctu8Phos"));
        uniRoles = compList.get(3);
        assertThat(uniRoles.inGroup(p1773_3521), equalTo(1));
        assertThat(uniRoles.inGroup(p563191_3), equalTo(1));
        assertThat(uniRoles.inGroup(p1262687_3), equalTo(1));
        assertThat(uniRoles.inGroup(p1104577_5), equalTo(1));
        assertThat(uniRoles.getName(), equalTo("root"));
        markerRoles = uniRoles.getRoles();
        assertThat(markerRoles, containsInAnyOrder("NLThreSynt", "LsuRiboProtL10e", "LsuRiboProtL18p", "LsuRiboProtL15p",
                "LsuRiboProtL14p", "LsuRiboProtL13a", "LsuRiboProtL11e", "SsuRiboProtS16e",
                "SsuRiboProtS15e", "SsuRiboProtS15a", "SsuRiboProtS13p", "SsuRiboProtS13e",
                "PhenTrnaSyntAlph", "SsuRiboProtS11p", "SsuRiboProtS11e", "ProtTranSubuSecy",
                "SsuRrnaNAdenNDime", "TrnaNMeth3", "LsuRiboProtL6p", "LsuRiboProtL3e",
                "LsuRiboProtL2p", "LsuRiboProtL1e", "SsuRiboProtS2p", "SsuRiboProtS2e"));
        // Test the searching.
        uniRoles = UniversalRoles.findGroup(s1104577_5, compList);
        assertThat(uniRoles.getName(), equalTo("root"));
        uniRoles = UniversalRoles.findGroup(s1773_3521, compList);
        assertThat(uniRoles.getName(), equalTo("R200 (Mycobacterium bovis AF2122/97)"));
        uniRoles = UniversalRoles.findGroup(s1262687_3, compList);
        assertThat(uniRoles.getName(), equalTo("R50 (Phascolarctobacterium sp. CAG:266)"));

    }

}
