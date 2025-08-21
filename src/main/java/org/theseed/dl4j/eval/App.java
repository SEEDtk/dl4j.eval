package org.theseed.dl4j.eval;

import java.util.Arrays;

import org.theseed.basic.BaseProcessor;

/**
 * DL4J-based genome evaluation
 *
 * train		train models to predict role frequencies
 * eval			evaluate a directory of GTO files
 * p3eval		evaluate a list of PATRIC genomes
 * gto			evaluate a single GTO
 * compare		produce an ORF-by-ORF comparison of genomes with their reference genomes
 * build		build a genome consistency evaluator from the PATRIC CoreSEED dump
 * analyze		analyze significant contributions among input features
 * p3All		create a master directory of evaluation GTOs for PATRIC prokaryotes
 * comp			produce completeness engine from master directory
 * mass			produce a summary evaluation report on a group of genomes
 * updateMass	update a summary evaluation report on a genome master directory
 * sort			sort the summary evaluation report
 * rRoles		create a report on the roles present or absent in a group of genomes
 * rTrain		process a role-training file to generate presence/absence classifiers
 * rBuild		generate role-training files from a representative-genome list file
 * binEval		evaluate binning results and optionally attempt to improve the bins
 * binReport	produce a standalone binning summary report
 *
 */
public class App
{

    /** static array containing command names and comments */
    protected static final String[] COMMANDS = new String[] {
             "train", "train models to predict role frequencies",
             "eval", "evaluate a directory of GTO files",
             "p3eval", "evaluate a list of PATRIC genomes",
             "gto", "evaluate a single GTO",
             "compare", "produce an ORF-by-ORF comparison of genomes with their reference genomes",
             "build", "build a genome consistency evaluator from the PATRIC CoreSEED dump",
             "analyze", "analyze significant contributions among input features",
             "p3All", "create a master directory of evaluation GTOs for PATRIC prokaryotes",
             "comp", "produce completeness engine from master directory",
             "mass", "produce a summary evaluation report on a group of genomes",
             "updateMass", "update a summary evaluation report on a genome master directory",
             "sort", "sort the summary evaluation report",
             "rRoles", "create a report on the roles present or absent in a group of genomes",
             "rTrain", "process a role-training file to generate presence/absence classifiers",
             "rBuild", "generate role-training files from a representative-genome list file",
             "binEval", "evaluate binning results and optionally attempt to improve the bins",
             "binReport", "produce a standalone binning summary report",
    };

    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Parse the parameters.
        switch (command) {
        case "analyze" -> processor = new AnalyzeProcessor();
        case "train" -> processor = new TrainProcessor();
        case "eval" -> processor = new EvalProcessor();
        case "p3eval" -> processor = new P3EvalProcessor();
        case "gto" -> processor = new GtoEvalProcessor();
        case "compare" -> processor = new CompareProcessor();
        case "build" -> processor = new BuildProcessor();
        case "p3All" -> processor = new P3AllProcessor();
        case "comp" -> processor = new CompletenessRolesProcessor();
        case "mass" -> processor = new MassEvalProcessor();
        case "sort" -> processor = new EvalSortProcessor();
        case "rTrain" -> processor = new RoleTrainProcessor();
        case "rRoles" -> processor = new RoleReportProcessor();
        case "rBuild" -> processor = new RoleBuildProcessor();
        case "binEval" -> processor = new BinEvalProcessor();
        case "binReport" -> processor = new BinReportProcessor();
        case "updateMass" -> processor = new UpdateMasterProcessor();
        case "-h", "--help" -> processor = null;
        default -> throw new RuntimeException("Invalid command " + command + ".");
        }
        if (processor == null)
            BaseProcessor.showCommands(COMMANDS);
        else {
            processor.parseCommand(newArgs);
            processor.run();
        }
    }
}

