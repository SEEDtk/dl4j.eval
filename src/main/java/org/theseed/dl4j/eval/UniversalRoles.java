/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.theseed.sequence.ProteinKmers;

/**
 * This object describes a completeness class and lists its universal roles.
 *
 * @author Bruce Parrello
 *
 */
public class UniversalRoles {

    // FIELDS
    /** minimum seed similarity for a genome to qualify (0 for the root group) */
    int score;
    /** seed protein kmers */
    ProteinKmers seedKmers;
    /** name of the group */
    String name;
    /** set of marker roles for this group */
    Set<String> markerRoles;

    /**
     * Create the universal role set for a specified group.
     *
     * @param score		minimum similarity score for group membership
     * @param name		name of the group
     * @param protein	seed protein sequence (use empty string for the root group)
     */
    public UniversalRoles(int score, String name, String protein) {
        this.score = score;
        this.name = name;
        this.seedKmers = new ProteinKmers(protein);
        this.markerRoles = new HashSet<String>(600);
    }

    /**
     * Add a new marker role to the list of this group's universal roles.
     *
     * @param role		ID of the role to add
     */
    public void addRole(String role) {
        this.markerRoles.add(role);
    }

    /**
     * Given a genome's seed protein, get its score for this group, or 0 if the genome does not belong.
     *
     * @param seed	seed protein kmers
     * @return the score for the protein relative to this group, or 0 if it does not belong
     */
    public int inGroup(ProteinKmers seed) {
        // Return 1 for the root group.
        int retVal = 1;
        if (this.score > 0) {
            // Here we are not the root group.
            retVal = this.seedKmers.similarity(seed);
            // Convert the score to 0 if it is too low.
            if (retVal < this.score)
                retVal = 0;
        }
        return retVal;
    }

    /**
     * @return TRUE if the specified role is a marker for this group, else FALSE
     *
     * @param role	ID of the role to check
     */
    public boolean contains(String role) {
        return this.markerRoles.contains(role);
    }

    /**
     * @return the number of marker roles for this group
     */
    public int size() {
        return this.markerRoles.size();
    }

    /**
     * @return the name of this group
     */
    public String getName() {
        return this.name;
    }

    /**
     * @return the list of roles in this group
     */
    public Set<String> getRoles() {
        return this.markerRoles;
    }

    /**
     * Load a list of universal role sets from a file.
     *
     * @param dataFile	input file to read (produced by kmers.reps.RoleProcessor)
     *
     * @return an ordered list of universal role sets from the file
     *
     * @throws IOException
     */
    public static List<UniversalRoles> Load(File dataFile) throws IOException {
        List<UniversalRoles> retVal = new ArrayList<UniversalRoles>(500);
        try (Scanner inStream = new Scanner(dataFile)) {
            inStream.useDelimiter("\\t|[\\r\\n]+\\s*");
            // Loop until we run out of file.
            while (inStream.hasNext()) {
                // We are positioned on a header.  Skip the genome ID and read the score.
                inStream.next();
                int score = inStream.nextInt();
                String name = inStream.next();
                String sequence = inStream.next();
                // Create the universal role object.
                UniversalRoles uniRoles = new UniversalRoles(score, name, sequence);
                // Loop through the member roles, stopping on the trailer.
                String line = inStream.next();
                while (! line.contentEquals("//")) {
                    uniRoles.addRole(line);
                    line = inStream.next();
                }
                // Push us onto the list.
                retVal.add(uniRoles);
            }
        }
        return retVal;
    }

    /**
     * @return the universal role group appropriate for the genome with the specified seed protein
     *
     * @param seedProtein	the seed protein for the relevant genome
     * @param compList		the comprehensive list of universal role groups
     */
    public static UniversalRoles findGroup(String seedProtein, List<UniversalRoles> compList) {
        // Parse the seed protein.
        ProteinKmers seedKmers = new ProteinKmers(seedProtein);
        // We want the highest similarity with the tightest group (based on minScore).  We are
        // guaranteed to find the root if nothing else matches.
        UniversalRoles retVal = null;
        int bestScore = 0;
        int bestSim = 0;
        for (UniversalRoles uniRoles : compList) {
            if (uniRoles.score >= bestScore) {
                // This group is the same or tighter than the current group.  Check the similarity.
                int sim = uniRoles.inGroup(seedKmers);
                if (sim > 0 && uniRoles.score > bestScore || sim > bestSim) {
                    // Here the new group is better.
                    retVal = uniRoles;
                    bestScore = uniRoles.score;
                    bestSim = sim;
                }
            }
        }
        return retVal;
    }

}
