/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.bamparsermaven;

import java.io.File;
import java.time.LocalDate;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * parses 10x bam file and serializes n requested top cells (most Umis) into a
 * java object file. Format is HashMap<GeneName,
 * Hashmap<CellBarcode,HashMap<UMIseq,UmiReadCount>>> @author rainer
 *
 */
public class Main {
    private static final String VERSION = "1.001";
    private static final String BUILD_DATE = "  Build Dec 4, 2019";
    //private static 
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
       // System.out.println(Main.class.getResource("/versionbuild.properties"));
       System.out.println("+++++++   ILLUMINAPARSER " + VERSION + BUILD_DATE + "++++++++++++");
       LocalDate exp = LocalDate.of(2020, 4, 4);
        if(LocalDate.now().isAfter(exp)){
            System.out.println("****************************************************************************************************************************************");
            System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  EXPIRED " + exp.toString() +  " - Download a new version from Github !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            System.out.println("****************************************************************************************************************************************");
            System.exit(0);
        }
        File inFile;
        File outFile;
        Options options = cli_otions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException ex) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("10x bam reader\n parses 10x bam files and creates Hashmaps that are used to identify corresponding Nanopore reads ", options);
            System.exit(1);
        }
        inFile = new File(cmd.getOptionValue("i"));
        if (inFile.exists() == false) {
            System.out.println(inFile.getName() + "input File/Directory does not exist");
            System.exit(1);
        }
        outFile = new File(cmd.getOptionValue("o"));
        if (outFile.exists()) {
            System.out.println(outFile.getName() + " !!!!!!!!!!!!!!!!!!!!!!!!!!!!  output File exists !!!!!!!!!!!!!!!!!!");
            System.exit(1);
        }
        Integer nCells = null;
        if (cmd.hasOption("n")) {
            nCells = Integer.parseInt(cmd.getOptionValue("n"));
        }
        String tsv = cmd.getOptionValue("t");
        File tsvFile = null;
        if (tsv != null) {
            tsvFile = new File(tsv);
            if(tsvFile.exists() == false || tsvFile.canRead() == false){
                System.out.println("tsv file with barcodes " + tsvFile.toString() + " does not exist or not readable");
                System.exit(1);
            }
        }
        if (nCells == null && tsvFile == null || nCells != null && tsvFile != null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Must provide cellranger tsv file(-t)", options);
            System.exit(1);
        }
        int windowsize;
        if (cmd.hasOption("w")) {
            windowsize = Integer.parseInt(cmd.getOptionValue("w"));
        } else {
            windowsize = 500;
        }
        new Parser(inFile, outFile, tsvFile, nCells, cmd.getOptionValue("g"),
                cmd.getOptionValue("b"), cmd.getOptionValue("u"))
                .parse(windowsize);
    }

    /**
     *
     * @return
     */
    private static Options cli_otions() {
        Options options = new Options();

        options.addOption(Option.builder("i").
                longOpt("inFileIllumina").
                required(true).
                desc("full path of 10x bam file").
                numberOfArgs(1)
                .build());
        options.addOption(Option.builder("o").
                longOpt("outFile").
                required(true).
                desc("full path of object output file").
                numberOfArgs(1)
                .build());
        options.addOption(Option.builder("w").
                longOpt("windowSize").
                required(false).
                desc("Window size for genomic regions, defaults to 500. Don't set too small -> risk that nanopore does not have alignmentblock in region").
                numberOfArgs(1)
                .build());
//        options.addOption(Option.builder("n").
//                longOpt("nCells").
//                required(false).
//                desc("use n cells with the most umis\nsupply this or the 10x tsv file with the list of cell BCs"
//                        + "bam file contains all barcodes also from drops without cells").
//                numberOfArgs(1)
//                .build());
        options.addOption(Option.builder("t").
                longOpt("tsv").
                required(true).
                desc("use this 10x tsv FILE to define the cell barcodes to use ").
                numberOfArgs(1)
                .build());
        options.addOption(Option.builder("b").
                longOpt("cellBCflag").
                required(true).
                desc("SAM tag for cell BC").
                numberOfArgs(1)
                .build());
        options.addOption(Option.builder("u").
                longOpt("umiFlag").
                required(true).
                desc("SAM tag for umi").
                numberOfArgs(1)
                .build());
        options.addOption(Option.builder("g").
                longOpt("geneFlag").
                required(true).
                desc("SAM tag for Gene name").
                numberOfArgs(1)
                .build());
        return options;
    }

}
