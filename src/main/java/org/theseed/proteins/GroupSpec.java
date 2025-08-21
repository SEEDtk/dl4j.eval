/**
 *
 */
package org.theseed.proteins;

import java.util.Collection;

import org.theseed.proteins.kmers.reps.RepGenome;
import org.theseed.sequence.ProteinKmers;

/**
 * This object specifies a completeness group.  It contains the representative-genome object (null for the
 * root group), the score limit, and the role matrix.  The objects sort by score (descending) and then
 * representative genome ID.  This means that if we output the groups in order, we can stop at the first
 * match.
 *
 * @author Bruce Parrello
 *
 */
public class GroupSpec implements Comparable<GroupSpec> {

    // FIELDS
    /** representative-genome object (NULL for the root group) */
    private final RepGenome repGenome;
    /** minimum acceptable score */
    private final int score;
    /** role matrix */
    private final RoleMatrix roleMtx;
    /** special repGenome object for root group */
    public static final RepGenome ROOT_GENOME = new RepGenome("root", "prokaryotic organism group",
            "MQHLNELIEKAKLAIESIQDKSLTALDEIRVEYFGKKGHFTQLMQELRNVSAEERPAMGA"
            + "KINEAKQAALEFLNTKKAEWEQAELNSKLEKERVDVSLPGRKVETGGLHPVTMTINRVTK"
            + "FFSELGFSVENGPEIESDYYNFDALNIPKHHPARADHDTFWFNPELLLRTQTSGVQIRTM"
            + "EKMQPPIRIMAPGRVYRNDYDQTHTPMFHQIELLYVDKKANFTELKGLLHDFLRAFFEED"
            + "LQVRFRPSYFPFTEPSAEVDVMGKNGKWLEVLGCGMVHLNVLRNVGIDPNEYSGFAVGMG"
            + "VERLTMLRYNVTDLRSFFENDLRFLKQFK");
    /** default role-matrix dimensions */
    private static final int MATRIX_SIZE = 100;

    /**
     * Create a group specification from the specified representative genome.
     *
     * @param repGenome		representative-genome object
     * @param score			minimum acceptable score
     */
    public GroupSpec(RepGenome repGenome, int score) {
        this.repGenome = repGenome;
        this.score = score;
        this.roleMtx = new RoleMatrix(MATRIX_SIZE, MATRIX_SIZE);
    }

    @Override
    public int compareTo(GroupSpec o) {
        int retVal = o.score - this.score;
        if (retVal == 0)
            retVal = this.repGenome.compareTo(o.repGenome);
        return retVal;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.repGenome == null) ? 0 : this.repGenome.hashCode());
        result = prime * result + this.score;
        return result;
    }


    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof GroupSpec)) {
            return false;
        }
        GroupSpec other = (GroupSpec) obj;
        if (this.repGenome == null) {
            if (other.repGenome != null) {
                return false;
            }
        } else if (!this.repGenome.equals(other.repGenome)) {
            return false;
        }
        return (this.score == other.score);
    }

    /**
     * @return TRUE if the specified genome is in this group, else FALSE
     *
     * @param seedProtein	seed protein of the specified genome
     */
    public boolean isClose(ProteinKmers seedProtein) {
        return (this.repGenome.similarity(seedProtein) >= this.score);
    }

    /**
     * Add the specified genome roles to this representative's role matrix.
     *
     * @param genomeId	ID of the source genome
     * @param roleIds	list of roles to add
     */
    public void recordGenome(String genomeId, Collection<String> roleIds) {
        this.roleMtx.register(genomeId, roleIds);
    }

    /**
     * @return the number of genomes in this group.
     */
    public int gCount() {
        return this.roleMtx.gCount();
    }

    /**
     * @return the role matrix for this group
     */
    public RoleMatrix getMatrix() {
        return this.roleMtx;
    }

    /**
     * @return the ID of this group's representative genome
     */
    public String getId() {
        return this.repGenome.getGenomeId();
    }

    /**
     * @return the display header for this group
     */
    public String getHeader() {
        return String.format("%s\t%d\tR%d (%s)\t%s", this.repGenome.getGenomeId(), this.score,
                this.score, this.repGenome.getName(), this.repGenome.getProtein());
    }

    @Override
    public String toString() {
        return String.format("R%d %s (%s)", this.score, this.repGenome.getGenomeId(), this.repGenome.getName());
    }

}
