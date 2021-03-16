/**
 *
 */
package org.theseed.proteins;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import org.theseed.counters.CountMap;

/**
 * This class maintains a matrix of roles to genomes.  A CountMap tracks the total number of occurrences
 * for each role. For each genome, there is a HashMap identifying the roles in it.
 *
 * @author Bruce Parrello
 *
 */
public class RoleMatrix {

    // FIELDS
    /** map of role IDs to counts */
    private CountMap<String> roleCounts;
    /** map of genome IDs to role hash sets */
    private Map<String, Set<String>> genomeRoles;
    /** expected number of roles per genome */
    private int roleCapacity;

    /**
     * Create a new role matrix with the specified capacities for genomes and roles per genome.
     *
     * @param genomes	estimated number of genomes
     * @param roles		estimated number of roles per genome
     */
    public RoleMatrix(int genomes, int roles) {
        this.roleCapacity = roles;
        this.roleCounts = new CountMap<String>();
        this.genomeRoles = new HashMap<String, Set<String>>(genomes);
    }

    /**
     * Register a genome and its roles.
     *
     * @param genome_id		ID of the genome to register
     * @param roles			collection of roles in the genome
     */
    public void register(String genome_id, Collection<String> roles) {
        Set<String> genomeSet = this.genomeRoles.get(genome_id);
        if (genomeSet == null) {
            genomeSet = new HashSet<String>(this.roleCapacity);
            this.genomeRoles.put(genome_id, genomeSet);
        }
        for (String role : roles)
            if (genomeSet.add(role)) this.roleCounts.count(role);
    }

    /**
     * @return the number of occurrences for a role
     *
     * @param role	role to check
     */
    public int roleCount(String role) {
        return this.roleCounts.getCount(role);
    }

    /**
     * @return TRUE if the specified genome contains the specified role
     *
     * @param genome_id		ID of the genome of interest
     * @param role			role of interest
     */
    public boolean rolePresent(String genome_id, String role) {
        boolean retVal = false;
        Set<String> genomeSet = this.genomeRoles.get(genome_id);
        if (genomeSet != null)
            retVal = genomeSet.contains(role);
        return retVal;
    }

    /**
     * @return the fraction of roles in the specified list found in the specified genome
     *
     * @param roles			collection of roles to check
     * @param genome_id		genome to check
     */
    public double completeness(Collection<String> roles, String genome_id) {
        double retVal = 0.0;
        Set<String> genomeSet = this.genomeRoles.get(genome_id);
        if (genomeSet != null && roles.size() > 0) {
            int count = 0;
            for (String role : roles)
                if (genomeSet.contains(role)) count++;
            retVal = ((double) count) / roles.size();
        }
        return retVal;
    }
    /**
     * Compute a marker role set.  A marker role set is a set of roles that occurs at a threshold fraction
     * or better in every genome.
     *
     * @return a marker role set for the genomes in this matrix.
     *
     * @param threshold		the minimum completeness fraction for the role set
     */
    public Collection<String> getMarkerRoles(double threshold) {
        // Get all the roles in order from most frequent to least frequent.
        List<String> retVal = this.roleCounts.sortedCounts().stream().map(CountMap<String>.Count::getKey).collect(Collectors.toList());
        // Loop until we run out of roles or meet the threshold.  If we fail, pop a role off the end.
        while (retVal.size() > 0 && ! this.threshold(retVal, threshold)) {
            retVal.remove(retVal.size() - 1);
        }
        // Return the result.
        return retVal;
    }

    /**
     * Compute a common role set.  A common role set is a set of roles each of which occurs in a given fraction
     * of the genomes.
     *
     * @return a common role set for the genomes in this matrix.
     *
     * @param threshold		the minimum frequency fraction for the role set
     */
    public Collection<String> getCommonRoles(double threshold) {
        // Create an empty list of roles.
        ArrayList<String> retVal = new ArrayList<String>(1000);
        // Fill it with roles that occur often enough.
        int minCount = (int) Math.ceil(this.genomeRoles.size() * threshold);
        for (CountMap<String>.Count counter : this.roleCounts.counts()) {
            if (counter.getCount() >= minCount) retVal.add(counter.getKey());
        }
        return retVal;
    }

    /**
     * @return TRUE if every genome meets the specified completeness threshold
     *
     * @param roles			list of roles
     * @param threshold		minimum fraction of roles that should be in every genome
     */
    private boolean threshold(Collection<String> roles, double threshold) {
        boolean retVal = true;
        // Iterate through all the genome IDs.
        Iterator<String> iter = this.genomeRoles.keySet().iterator();
        // Loop until we fail or run out of genomes.
        while (retVal && iter.hasNext()) {
            String genome = iter.next();
            double actual = this.completeness(roles, genome);
            retVal = (actual >= threshold);
        }
        return retVal;
    }

    /**
     * @return a list of all the genomes in this matrix
     */
    public Set<String> genomes() {
        return this.genomeRoles.keySet();
    }

    /**
     * @return the number of genomes in this matrix
     */
    public int gCount() {
        return this.genomeRoles.size();
    }

}
