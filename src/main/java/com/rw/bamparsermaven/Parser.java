/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.bamparsermaven;

import com.rw.nuc.encoding.TwoBit.NucleicAcidTwoBitPerBase;
import com.rw.nuc.reads.Illumina.All10xselectedCells;
import com.rw.nuc.reads.Illumina.Hashing;
import com.rw.nuc.reads.Illumina.ParsedIlluminaData;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import it.unimi.dsi.fastutil.longs.LongHash;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author rainer
 */
public class Parser {

    private final static String cellsUMIsLogFileSuffix = ".cellsUmicounts.txt";
    private final File inFile;
    private final File outFile;
    private final File tsvFile;
    private final Integer nCells;
    private final String illuminaGeneFlag;
    private final String cellBCFlag;
    private final String umiFlag;

    public Parser(File inFile, File outFile, File tsvFile, Integer nCells,
            String illuminaGeneFlag, String cellBCFlag, String umiFlag) {
        this.inFile = inFile;
        this.outFile = outFile;
        this.tsvFile = tsvFile;
        this.nCells = nCells;
        this.illuminaGeneFlag = illuminaGeneFlag;
        this.cellBCFlag = cellBCFlag;
        this.umiFlag = umiFlag;
    }

    /**
     *
     * @param windowSize
     */
    public void parse(int windowSize) {
        //int nUMIsFound = 0;
        //File logFilePath = outFile.getParentFile();
        //key is cell BC , value is list of all UMIs

        //all cells in 10x tsv file
        All10xselectedCells all10xselCells;
//        /**
//         * hashmap with all UMIs for each 10x selected cell
//         */
//        TLongLongHashMap allUmisEachCell = new TLongLongHashMap();

        //for each cell all UMIs that were found - for checking reads where cell BC found but no gene ... key is ell bc , value is set of UMIs
        //Long2ObjectOpenHashMap<LongOpenHashSet> allUMIsforEachCell = new Long2ObjectOpenHashMap<>(all10xselCells.size());
        long starttime = System.currentTimeMillis();
        //******************************    Generate data for serialized object
        //HashSet<NucleicAcidTwoBitPerBase> dummyCellCounting = new HashSet<>();
        ParsedIlluminaData illuminaGeneDat = null;
        int debugreads = 0;
        //int debugNcells =0;
        final SamReaderFactory factory
                = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
        SamReader sr = null;
        try {
            sr = factory.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(inFile), 10000000)));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        }
        SAMRecordIterator samReadIterator = sr.iterator(); //define maxscore for adapter

        while (samReadIterator.hasNext()) {
            debugreads++;
            SAMRecord sam = samReadIterator.next();
            //gene
            String geneAttribute = (String) sam.getAttribute(illuminaGeneFlag);
            String cellString = sam.getStringAttribute(cellBCFlag);
            String umiString = sam.getStringAttribute(umiFlag);
            if (illuminaGeneDat == null)// need cellBC length to construct it
            {
                all10xselCells = tsvFile != null ? getCellsToUseFromCellRangerTSV(tsvFile, Hashing.get_LONG_HASH_STRATEGY(cellString.length() - 2))
                        : getCellsToUseDefineNcells(nCells, Hashing.get_LONG_HASH_STRATEGY(cellString.length() - 2));
                illuminaGeneDat = new ParsedIlluminaData(all10xselCells, windowSize);
            }
            illuminaGeneDat.addSamRecord(sam, geneAttribute, cellString, umiString);
        }

        OutputStream streamOut;
        ObjectOutputStream oos = null;
        try {
            streamOut = new BufferedOutputStream(new FileOutputStream(outFile), 1000000);
            oos = new ObjectOutputStream(streamOut);
            oos.writeObject(illuminaGeneDat);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        }
        if (oos != null) {
            try {
                oos.close();
            } catch (IOException ex) {
                Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        System.out.println("Scan took: " + (new SimpleDateFormat("mm:ss:SSS")).format(new Date(System.currentTimeMillis() - starttime)));
        System.out.println("Getting some stats");
        int umiCountGenes = illuminaGeneDat.illuminaGenesData.entrySet().stream().map(x->x.getValue()).
                flatMap(v -> v.stream()).mapToInt(f-> f.size()).sum();
        System.out.println("Found " + umiCountGenes + " UMIs associated with genes.");
        if(umiCountGenes == 0)
             System.out.println("!!!!!!!! WARNING - NO UMI ASSOCIATED WITH GENES FOUND -- CHECK PARAMETERS !!!!!!");
        //int umiCountRegions = illuminaGeneDat.illuminaChromosomesData.values().stream().flatMap(x -> x.stream()).
        //       System.out.println("took " + (System.currentTimeMillis() - starttime) / 1000 + " secs\n type enter to exit");
//        System.out.println("TES reading file");
//        //illuminaData = null;
//        IlluminaData dummy = null;
//        InputStream streamIn;
//        ObjectInputStream ois = null;
//
//        //int debugNcells =0;
//        try {
//            streamIn = new BufferedInputStream(new FileInputStream(outFile), 1000000);
//            ois = new ObjectInputStream(streamIn);
//        } catch (FileNotFoundException ex) {
//            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
//        } catch (IOException ex) {
//            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        try {
//            dummy = (IlluminaData) ois.readObject();
//            ois.close();
//        } catch (IOException ex) {
//            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
//        } catch (ClassNotFoundException ex) {
//            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        System.out.println("Ois reading took " + (System.currentTimeMillis() - starttime) / 1000 + " secs\n type enter to exit");
//        Scanner sc = new Scanner(System.in);
//        String s = sc.nextLine();
        //write XML
//        if (this.createXML) {
//            XMLEncoder encoder = null;
//            try {
//                encoder = new XMLEncoder(new BufferedOutputStream(new FileOutputStream("D:\\Nanopore_13.6\\chr11.xml")));
//            } catch (FileNotFoundException fileNotFound) {
//                System.out.println("ERROR: While Creating or Opening the File");
//            }
//            encoder.writeObject(illuminaData);
//            encoder.close();
//        }
    }

    /**
     *
     * @param tsvFile cellranger file with list of cell BC use
     * @return
     */
    private All10xselectedCells getCellsToUseFromCellRangerTSV(File tsvFile, LongHash.Strategy strategy) {
        All10xselectedCells retval;
        retval = new All10xselectedCells();
        BufferedReader reader;
        try {
//            InputStream is;
//            if(tsvFile.getName().endsWith(".gz")){
//             is = new GZIPInputStream(new FileInputStream(tsvFile));   
//            } else
//                is = new FileInputStream(tsvFile);
//            reader = new BufferedReader(new InputStreamReader(is));
           reader = new BufferedReader(new FileReader(tsvFile));
            String s;
            while ((s = reader.readLine()) != null) {
                int dashindex = s.indexOf('-');//CB contains -1 in the end
                NucleicAcidTwoBitPerBase cellbc = new NucleicAcidTwoBitPerBase(dashindex == -1 ? s : s.substring(0, dashindex));
                retval.add(cellbc.getSequence());
            }
            reader.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        }
        return retval;
    }

    /**
     * count Umis per cell and define n cells with most umis
     *
     * @param n
     * @return
     */
    private All10xselectedCells getCellsToUseDefineNcells(int n, LongHash.Strategy strategy) {
        long starttime = System.currentTimeMillis();
        //don't treat UMIs per gene, simply count UMIs for cells, corrected 10 nt UMIs should not collide for transcripts of cells
        //key is cell BC , value is Hashset of UMIs - entered just the long value not the object to use efficient Trove TlongHashset instead of HashSet<NucleicAcidTwoBitPerBase>
        HashMap<NucleicAcidTwoBitPerBase, LongOpenHashSet> umiCountingMap = new HashMap<>();
        final SamReaderFactory factory
                = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader sr = null;
        try {
            sr = factory.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(inFile), 10000000)));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        }
        //only use entries with just one gene (uniquely matching)
        //htsjdk samfile reading is by far the most rate limiting step (>95% of CPUtime )
        sr.iterator().stream().filter((r) -> {
            String s = (String) r.getAttribute("GN");
            return s != null && s.split(";").length == 1;
        }).forEach(
                //sr.iterator().forEachRemaining(
                (r) -> {
                    String cellString = r.getStringAttribute("CB");
                    String u = r.getStringAttribute("UB");
                    if (cellString != null && u != null) {
                        int dashindex = cellString.indexOf('-');//CB contains -1 in the end
                        NucleicAcidTwoBitPerBase cellbc = new NucleicAcidTwoBitPerBase(dashindex == -1 ? cellString : cellString.substring(0, dashindex));
                        LongOpenHashSet oneCellData;
                        if ((oneCellData = umiCountingMap.get(cellbc)) == null) {
                            umiCountingMap.put(cellbc, (oneCellData = new LongOpenHashSet()));
                        }
                        oneCellData.add((new NucleicAcidTwoBitPerBase(u)).getSequence());//just add the long val for hashing to save ram
                    }
                });
        try {
            sr.close();
        } catch (IOException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("parsing BAM took " + (new SimpleDateFormat("mm:ss:SSS")).format(new Date(System.currentTimeMillis() - starttime)));
        starttime = System.currentTimeMillis();
        List<Integer> umiCountsList = new ArrayList<>();//a sorted list for median umi calculation
        List<NucleicAcidTwoBitPerBase> bcList;
        bcList = umiCountingMap.entrySet()
                .stream()
                .sorted((e1, e2) -> new Integer(e2.getValue().size()).compareTo(e1.getValue().size())) // custom Comparator                    
                .limit(nCells)
                .peek(e -> umiCountsList.add(e.getValue().size())) // need this one for stats
                .map(e -> e.getKey())
                .collect(Collectors.toCollection(ArrayList::new)); //https://stackoverflow.com/questions/30425836/java-8-stream-map-to-list-of-keys-sorted-by-values 
        saveSelectedCellsInfo(bcList, umiCountsList); // save umi counts for selected cells
       All10xselectedCells cellsToUse = new All10xselectedCells();
        bcList.forEach(e -> cellsToUse.add(e.getSequence()));
//            all10xselectedCells = umiCountsCellBCmap.entrySet().stream()
//                    .peek(e -> umiCountsList.add(e.getKey()))
//                    .map(e -> e.getValue())
//                    .collect(Collectors.toCollection(HashSet::new)); //https://stackoverflow.com/questions/30425836/java-8-stream-map-to-list-of-keys-sorted-by-values           

//            List<Integer> listForUMIcountstats = new ArrayList<>();     
//   .collect(Collectors.toMap(p -> p.getValue().size(), p -> p.getKey()));    
        System.out.println("nCells : \t" + cellsToUse.size());
        System.out.println("total Umis : \t" + umiCountsList.stream().mapToInt(Integer::intValue).sum());
        System.out.println("median Umis/cell : \t" + umiCountsList.get(Math.round(umiCountsList.size() / 2)));
        System.out.println("mean Umis/cell : \t" + umiCountsList.stream().reduce(0, Integer::sum) / umiCountsList.size());
        //System.out.println("took " + (new SimpleDateFormat("mm:ss:SSS")).format(new Date(System.currentTimeMillis() - starttime)));
        return cellsToUse;
    }

    /**
     * writes file with cell BCs, Umi counts
     *
     *
     */
    private void saveSelectedCellsInfo(List<NucleicAcidTwoBitPerBase> bcs, List<Integer> counts) {
        File out = new File(this.outFile.getAbsolutePath() + cellsUMIsLogFileSuffix);
        BufferedWriter writer;
        try {
            writer = new BufferedWriter(new FileWriter(out));
            writer.write("Cell BC \t nUMIs");
            writer.newLine();
            for (int i = 0; i < bcs.size() && i < counts.size(); i++) {
                writer.write(bcs.get(i).toString() + "\t" + counts.get(i));
                writer.newLine();
            }
            writer.close();
        } catch (IOException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
