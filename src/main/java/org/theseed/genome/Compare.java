/**
 *
 */
package org.theseed.genome;

import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.theseed.locations.Location;

/**
 * This class contains utilities for genome comparison.
 *
 * @author Bruce Parrello
 *
 */
public class Compare {

    // FIELDS

    /** comparator for sorting */
    private Comparator<Feature> orfSorter;

    /** number of ORFS annotated in the old genome but not the new one */
    private int oldOnly;
    /** number of ORFs annotated in the new genome but not the old one */
    private int newOnly;
    /** number of ORFs identical in both genomes */
    private int identical;
    /** number of ORFs that are commonly annotated but have different functions */
    private int differentFunctions;
    /** number of ORFs that have the same function, but are longer in the new genome */
    private int longer;
    /** number of ORFs that have the same function, but are shorted in the new genome */
    private int shorter;


    /** Initialize a comparison processor */
    public Compare() {
        this.orfSorter = new OrfSorter();
    }

    /**
     * Comparator for sorting Features by contig, strand, end point.
     */
    public static class OrfSorter implements Comparator<Feature> {

        @Override
        public int compare(Feature arg0, Feature arg1) {
            Location loc0 = arg0.getLocation();
            Location loc1 = arg1.getLocation();
            int retVal = loc0.getContigId().compareTo(loc1.getContigId());
            if (retVal == 0) {
                retVal = loc0.getDir() - loc1.getDir();
                if (retVal == 0)
                    retVal = loc0.getEnd() - loc1.getEnd();
            }
            return retVal;
        }
    }

    /**
     * Compare the two genomes and store the statistics in this object.
     *
     * @param newGenome		new genome
     * @param compareGenome	original genome with different annotations
     *
     * @return TRUE if successful, FALSE if the genomes cannot be compared.
     */
    public boolean compare(Genome newGenome, Genome oldGenome) {
        boolean retVal = true;
        // Determine if the genomes can be compared.
        for (Contig oldContig : oldGenome.getContigs()) {
            if (newGenome.getContig(oldContig.getId()) == null)
                retVal = false;
        }
        if (retVal) {
            // The genomes are comparable. Load all their features into a sorted set.
            SortedSet<Feature> newOrfs = this.sortFeatures(newGenome);
            SortedSet<Feature> oldOrfs = this.sortFeatures(oldGenome);
            // Clear the counters.
            this.oldOnly = 0;
            this.newOnly = 0;
            this.identical = 0;
            this.differentFunctions = 0;
            this.longer = 0;
            this.shorter = 0;
            // Get iterators through the sets.
            Iterator<Feature> newIter = newOrfs.iterator();
            Iterator<Feature> oldIter = oldOrfs.iterator();
            if (! newIter.hasNext()) {
                this.oldOnly = oldOrfs.size();
            } else if (! oldIter.hasNext()) {
                this.newOnly = oldOrfs.size();
            } else {
                // We have ORFs to compare. Prime the main loop.
                Feature oldFeature = newIter.next();
                Feature newFeature = newIter.next();
                while (newIter.hasNext() && oldIter.hasNext()) {
                    int comp = this.orfSorter.compare(oldFeature, newFeature);
                    if (comp < 0) {
                        // Old feature is an orphan.
                        this.oldOnly++;
                        oldFeature = oldIter.next();
                    } else if (comp > 0) {
                        // New feature is an orphan.
                        this.newOnly++;
                        newFeature = newIter.next();
                    } else {
                        // Both features match.  Check the annotations.
                        if (! newFeature.getFunction().contentEquals(oldFeature.getFunction())) {
                            differentFunctions++;
                        } else {
                            // Annotations match.  Check the lengths.
                            comp = newFeature.getLocation().getLength() - oldFeature.getLocation().getLength();
                            if (comp < 0) {
                                this.shorter++;
                            } else if (comp > 0) {
                                this.longer++;
                            } else {
                                this.identical++;
                            }
                        }
                        // Advance both features.
                        oldFeature = oldIter.next();
                        newFeature = newIter.next();
                    }
                }
                // Run out both iterators.
                while (newIter.hasNext()) {
                    newIter.next();
                    newOnly++;
                }
                while (oldIter.hasNext()) {
                    oldIter.next();
                    oldOnly++;
                }
            }
        }
        return retVal;
    }

    /**
     * Get all the protein features of the genome sorted by ORF.
     *
     * @param genome	genome whose features are to be processed
     *
     * @return a sorted set of all the pegs in the genome
     */
    private SortedSet<Feature> sortFeatures(Genome genome) {
        SortedSet<Feature> retVal = new TreeSet<Feature>(this.orfSorter);
        retVal.addAll(genome.getPegs());
        return retVal;
    }

    /**
     * @return the number of ORFs only annotated in the old genome
     */
    public int getOldOnly() {
        return oldOnly;
    }

    /**
     * @return the number of ORFs only annotated in the new genome
     */
    public int getNewOnly() {
        return newOnly;
    }

    /**
     * @return the number of ORFs identically annotated in both genomes
     */
    public int getIdentical() {
        return identical;
    }

    /**
     * @return the number of ORFs annotated in both genomes, but with different functions
     */
    public int getDifferentFunctions() {
        return differentFunctions;
    }

    /**
     * @return the number of ORFs annotated identically, but are longer in the new genome
     */
    public int getLonger() {
        return longer;
    }

    /**
     * @return the number of ORFs annotated identically, but are shorter in the new genome
     */
    public int getShorter() {
        return shorter;
    }

    /**
     * @return the number of ORFs annotated in both genomes
     */
    public int getCommon() {
        return (identical + differentFunctions + longer + shorter);
    }

}
