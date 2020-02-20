package org.theseed.dl4j.eval;

import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * DL4J-based genome evaluation
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
        case "eval" :
            processor = new EvalProcessor();
            break;
        case "p3eval" :
            processor = new P3EvalProcessor();
            break;
        case "gto" :
            processor = new GtoEvalProcessor();
            break;
        case "compare" :
            processor = new CompareProcessor();
            break;
        case "bins" :
            processor = new BinProcessor();
            break;
        case "build" :
        	processor = new BuildProcessor();
        	break;
        default :
            throw new RuntimeException("Invalid command " + command + ".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}

