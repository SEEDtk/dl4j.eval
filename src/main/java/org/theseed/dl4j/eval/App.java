package org.theseed.dl4j.eval;


/**
 * Hello world!
 *
 */
public class App
{
    public static void main( String[] args )
    {
        TrainProcessor runObject = new TrainProcessor();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }
    }
}
