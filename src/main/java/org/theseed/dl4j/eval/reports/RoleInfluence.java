/**
 *
 */
package org.theseed.dl4j.eval.reports;

/**
 * This is a simple class that contains a role ID, a role name, and a floating-point number.  It is used
 * to accumulate the significane of a role.
 *
 * @author Bruce Parrello
 *
 */
public class RoleInfluence implements Comparable<RoleInfluence> {

    // FIELDS
    /** ID of the role */
    private String id;
    /** name of the role */
    private String name;
    /** influence total */
    private double rating;

    /**
     * Construct a role influence object from an ID and name.
     *
     * @param id		ID of this role
     * @param name		name of this role
     */
    public RoleInfluence(String id, String name) {
        this.id = id;
        this.name = name;
        this.rating = 0.0;
    }

    /**
     * We sort from highest rating to lowest, then by role ID.
     */
    @Override
    public int compareTo(RoleInfluence o) {
        int retVal = Double.compare(o.rating, this.rating);
        if (retVal == 0)
            retVal = this.id.compareTo(o.id);
        return retVal;
    }

    /**
     * Increment the role rating.
     *
     * @param incr	value to add to the rating
     */
    public void increment(double incr) {
        this.rating += incr;
    }

    /**
     * @return the role id
     */
    public String getId() {
        return this.id;
    }

    /**
     * @return the role name
     */
    public String getName() {
        return this.name;
    }

    /**
     * @return the role rating
     */
    public double getRating() {
        return this.rating;
    }

}
