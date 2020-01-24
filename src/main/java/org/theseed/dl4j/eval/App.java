package org.theseed.dl4j.eval;

import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * Hello world!
 *
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        // Parse the parameters.
        switch (command) {
        case "train" :
            processor = new TrainProcessor();
            break;
        case "predict" :
            processor = new PredictProcessor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ": must be \"train\" or \"predict\".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}

