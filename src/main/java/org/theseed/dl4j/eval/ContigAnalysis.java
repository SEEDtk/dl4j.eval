/**
 *
 */
package org.theseed.dl4j.eval;

import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.locations.Location;

import j2html.tags.DomContent;
import static j2html.TagCreator.*;

import java.util.HashSet;
import java.util.Set;

/**
 * This object contains analytical information about a contig, and is used to generate comments about it.
 *
 * @author Bruce Parrello
 *
 */
public class ContigAnalysis {

    // FIELDS
    /** the base contig being analyzed */
    private Contig contig;
    /** TRUE if this contig is unusually short */
    private boolean tooShort;
    /** number of good features */
    private int goodFeatures;
    /** list of bad features */
    private Set<Feature> badFeatures;

    private static final int CONTIG_EDGE = 6;

    /**
     * Create the analysis object for a contig.
     *
     * @param contig	contig to analyze
     * @param gReport	quality report for the relevant genome
     */
    public ContigAnalysis(Contig contig, GenomeStats gReport) {
        this.contig = contig;
        this.tooShort = (contig.length() < gReport.getN70());
        this.goodFeatures = 0;
        this.badFeatures = new HashSet<Feature>();
    }

    /**
     * Count a feature on this contig.
     *
     * @param feat		feature of interest
     * @param isGood	status of the feature
     */
    public void countFeature(Feature feat, GenomeStats.FeatureStatus isGood) {
        switch (isGood) {
        case GOOD:
            this.goodFeatures++;
            break;
        case BAD:
            this.badFeatures.add(feat);
            break;
        default:
            break;
        }
    }

    /** @return a comment about this feature's relationship to this contig
     *
     * @param loc	location of the feature
     */
    public DomContent locationComment(Location loc) {
        // Determine how many edges we are close to.
        boolean leftEdge = (loc.getLeft() <= CONTIG_EDGE);
        boolean rightEdge = (loc.getRight() + CONTIG_EDGE >= this.contig.length());
        boolean forward = (loc.getDir() == '+');
        String position;
        if (leftEdge && rightEdge) {
            position = "fills contig";
        } else if (leftEdge && forward || rightEdge && ! forward) {
            position = "starts near the edge of contig";
        } else if (leftEdge || rightEdge) {
            position = "ends near the edge of contig";
        } else {
            position = "is in contig";
        }
        // Determine what we know about the contig.
        PhraseBuilder knowledge = new PhraseBuilder(", which", ".");
        if (this.tooShort) {
            knowledge.addPhrase("is short");
        }
        if (this.goodFeatures <= 0) {
            knowledge.addPhrase("has no good features");
        }
        // String it all together.
        DomContent retVal = join(text(position), this.link(), text(knowledge.toString()));
        return retVal;
    }

    /**
     * @return the list of bad features in this contig
     */
    public Set<Feature> getBadFeatures() {
        return this.badFeatures;
    }

    /**
     * @return the name of this contig, or optionally a hyperlink to it
     */
    public DomContent link() {
        return text(this.contig.getId());
    }

}
