/**
 *
 */
package org.theseed.dl4j.eval.reports;

/**
 * This interface is used for reporters that support reference genomes.
 *
 * @author Bruce Parrello
 *
 */
public interface IRefReporter {

    /**
     * Specify the engine for computing reference genomes.
     *
     * @param refEngine		object for computing reference genome IDs
     */
    public void setEngine(RefGenomeComputer refEngine);


}
