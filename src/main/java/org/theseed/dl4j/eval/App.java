package org.theseed.dl4j.eval;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * DL4J-based genome evaluation
 *
 * train	train models to predict role frequencies
 * eval		evaluate a directory of GTO files
 * p3eval	evaluate a list of PATRIC genomes
 * gto		evaluate a single GTO
 * compare	produce an ORF-by-ORF comparison of genomes with their reference genomes
 * bins		evaluate binning results
 * build	build a genome consistency evaluator from training data
 * analyze	analyze significant contributions among input features
 *
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Parse the parameters.
        switch (command) {
        case "analyze" :
            processor = new AnalyzeProcessor();
            break;
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

