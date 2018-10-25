// -----------------------------------------------
// ---
// ---             Victor Jacobson
// ---            
// ---        Naive Sequential Search
// ---     
// ---    searhces kmer for exact matches 
// ---         in one direction only
// ---     
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class DNA_NaiveSequentialSearch
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        print("DNA_NaiveSequentialSearch!");
        
        DNA_NaiveSequentialSearch data = new DNA_NaiveSequentialSearch();
        
        data.newGenomeFromFASTA("genomes/Zika.NC_012532.1.fasta");
        
        data.print();
        
        String stringToSearch = "AAAAA";
        
        ArrayList<Integer> matches = data.searchGenome(stringToSearch);
        
        
        for(int i = 0; i < matches.size(); i++)
        {
            printl(matches.get(i));
        }
        printl("There were " + matches.size() + " matches for " + stringToSearch);
    }
    public static void line(){System.out.println(); }
    public static void printl(Object o){System.out.println(o); }
    public static void print(Object o){System.out.print(o); }
    public static int randomInt(int max){return (int)(Math.random() * max);}
    

    public String genome;
    
    public DNA_NaiveSequentialSearch() { }
    
    public void newGenomeFromFASTA(String file_name)
    {
        StringBuilder stringB = new StringBuilder();
        
        File data = new File(file_name);
        try
        {
            Scanner in = new Scanner(data);
      
            line();
            printl("Reading: " + in.nextLine());
            
            while(in.hasNext())
            {
                stringB.append(in.nextLine());
            }
            
            printl(this.toString() + " contains: " + file_name);
        }
        catch(FileNotFoundException fnfe)
        {
            fnfe.printStackTrace();
        }
        
        genome = stringB.toString();
    }
    
    public ArrayList<Integer> searchGenome(String kmer)
    {
        if(genome.length() == 0)
        {
            printl("The genome is empty!");
            return null;
        }
        
        ArrayList<Integer> matchIndexes = new ArrayList<Integer>();
        
        for(int i = 0; i < genome.length() - kmer.length(); i++)
        {
            // Search for the first character match
            boolean match = true;
            for(int j = 0; j < kmer.length(); j++)
            {
                // if it matches, check each subsequent character for a match
                if(genome.charAt(i + j) != kmer.charAt(j))
                {
                    match = false;
                    continue;
                }
            }
            // if we get to the end of the kmer and have a match, record the index of the match
            if(match)
            {
                matchIndexes.add(i);
            }
        }
        return matchIndexes;
    }

    public void print()
    {
        printl(genome);
    }
}