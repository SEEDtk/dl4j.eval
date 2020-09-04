/**
 *
 */
package org.theseed.dl4j.eval;

import java.io.File;

/**
 * @author Bruce Parrello
 *
 */
public interface IConsistencyChecker {

    /**
     * @return the directory containing the role models
     */
    File getRoleDir();

    /**
     * Store the actual count of a role
     *
     * @param iRole		index of the role in the role list
     * @param role		ID of the role
     * @param iGenome	index of the genome in the matrix
     * @param count		number of occurrences found
     */
    void storeActual(int iRole, String role, int iGenome, int count);



}
