// -----------------------------------------------
// ---
// ---       Kean GSSRP 2016 Bioinformatics
// ---             Victor Jacobson
// ---               
// ---     Tests seq. search and index suffix tree
// ---     with reverse and complimentary searching
// ---
// -----------------------------------------------
import java.util.*;
import java.util.Date;
import java.io.*;

public class GnomeRCompSearchTester
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        printl("Welcome to Search Tester 0.01 !");

        //String genome_file = "genomes/Zika.NC_012532.1.fasta";
        String genome_file = "genomes/EcoliK12DH10B.fasta";

        String log_file = ("searchresults/GenomeSearchLog" + new Date().toString() + ".dat").replaceAll("\\s","").replaceAll(":","");

        File file = new File(log_file);
        PrintWriter writer = null;

        printl("Save Log File? y/n");
        
        if( in.next().charAt(0) == 'y' )
        {
            try
            {
                file.createNewFile();

                writer = new PrintWriter( file );  

                printl("Createed log file: " + log_file);
            }
            catch(Exception e)
            {
                printl("Log file creation failed!");
            } 

            if(writer != null)
            {
                writer.println("Kean GSSRP 2016 Bioinformatics");
                writer.println(" Victor Jacobson ");
                writer.println("Searching " + genome_file + " using:");
                writer.println("Sequential Search, Naive Suffx Tree, and Optimized Suffix Tree.");
            }
        }

        DNA_NSS_RComp zikaSeqSearch = new DNA_NSS_RComp();
        zikaSeqSearch.newGenomeFromFASTA(genome_file);

        DNA_OST_RComp zikaOptSuffixTree = new DNA_OST_RComp(genome_file);

        KmerGenerator kmers = new KmerGenerator();
        NanoTimer timer = new NanoTimer();

        while( true )
        {
            ArrayList<String> printQueue = new ArrayList<String>();
            printl("Enter a kmer length to generate and search:");
            int length;
            if(in.hasNextInt())
                length = in.nextInt();
            else
                break;

            String kmer = kmers.getKmer(length);
            printQueue.add("Search length: " + kmer.length() + " - " + kmer);

            // Print results
            timer.reset();
            timer.start();
            zikaSeqSearch.searchGenome(kmer);
            timer.stop();
            printQueue.add("  NSS: \t" + //results.size() + " matches in " + 
                timer.elapsedFormated()// + results.toString()
            );

            timer.reset();
            timer.start();
            zikaOptSuffixTree.searchGenomeRComp(kmer);
            timer.stop();
            printQueue.add("  OST: \t" + //results.size() + " matches in " + 
                timer.elapsedFormated()// + results.toString()
            );
            
            zikaSeqSearch.printResults();
            zikaOptSuffixTree.printResults();
            
            for(String s : printQueue)
                printl(s);

            if(writer != null)
            {
                for(String s : printQueue)
                    writer.println(s);
            }

            printQueue.clear();
        }
        if(writer != null)
            writer.close();

        printl("Finished");
    }

    public static void line(){System.out.println(); }

    public static String printl(Object o){System.out.println(o); return (String)o; }

    public static void print(Object o){System.out.print(o); }

    public static int randomInt(int max){return (int)(Math.random() * max);}
}