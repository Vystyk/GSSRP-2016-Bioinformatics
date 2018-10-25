// -----------------------------------------------
// ---
// ---             Victor Jacobson
// ---            
// ---        Naive Sequential Search
// ---     
// ---    searhces kmer for exact matches in 
// ---    both directions and complementary matches
// ---     
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class DNA_NSS_WildCodes
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        print("DNA_NSS_WildCodes! - Reverse and complimentary searches");

        DNA_NSS_WildCodes data = new DNA_NSS_WildCodes();

        //String genome_file = "genomes/Zika.NC_012532.1.fasta";
        String genome_file = "genomes/EcoliK12DH10B.fasta";
        
        data.newGenomeFromFASTA(genome_file);

        //data.print();

        printl("F = forard, RC = Reverse Compliment");

        String kmer = "AGCCTGA";

        printl("Searching for: " + kmer);

        data.searchGenome(kmer);

        data.printResults();
    }

    public static void line(){System.out.println(); }

    public static void printl(Object o){System.out.println(o); }

    public static void print(Object o){System.out.print(o); }

    public static int randomInt(int max){return (int)(Math.random() * max);}

    public String genome;

    public String lastSearch = "";
    public ArrayList<Match> matches;

    public DNA_NSS_WildCodes() { }

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

    public void printResults()
    {
        if( lastSearch.length() == 0)
        {
            printl("Nothing has been searched yet.");
            return;
        }

        if( matches.size() == 0)
        {
            printl("No matches for: " + lastSearch);
            return;
        }

        for(int i = 0; i < matches.size(); i++)
        {
            printl(matches.get(i).toString());// + ", ");
        }
        printl("There were " + matches.size() + " matches for " + lastSearch);
    }

    public void printDetailedResults()
    {
        if( lastSearch.length() == 0)
        {
            printl("Nothing has been searched yet.");
            return;
        }

        if( matches.size() == 0)
        {
            printl("No matches for: " + lastSearch);
            return;
        }

        for(int i = 0; i < matches.size(); i++)
        {
            printl(matches.get(i).toStringDetailed());
        }
        printl("There were " + matches.size() + " matches for " + lastSearch);
    }

    public void searchGenome(String kmer)
    {
        lastSearch = kmer;

        if(genome.length() == 0)
        {
            printl("The genome is empty!");
            return;
        }

        matches = new ArrayList<Match>();

        for(int i = 0; i < genome.length() - kmer.length(); i++)
        {
            // Search for the first character match
            boolean match_f = true;
            boolean match_rc = true;

            int kmer_rev_len = kmer.length() - 1;
            
            //Keep track of match exactness. Zero = perfect, higher numbers = less perfect 
            int matchQuality_F = 0;
            int matchQuality_RC = 0;
            
            byte current_check = -1;

            // if it matches, check each subsequent character for a match
            for(int j = 0; j < kmer.length(); j++)
            {
                // Find forward match
                current_check = DNA.compare( kmer.charAt(j), genome.charAt(i + j));
                if( current_check == -1 )
                    match_f = false;
                else
                    matchQuality_F += current_check;

                // Find complimentary reverse match
                current_check = DNA.compare( DNA.comp(kmer.charAt(kmer_rev_len - j)), genome.charAt(i + j));
                if( current_check == -1 )
                    match_rc = false;
                else
                    matchQuality_RC += current_check;
                    
                // If both directions are false, stop searching this index
                if( !match_f && !match_rc )
                    continue;

            }
            // if we get to the end of the kmer and have a match, record the index of the match
            if(match_f)
            {
                matches.add(new Match("F", i, kmer.length(), matchQuality_F));
            }

            if(match_rc)
            {
                matches.add(new Match("RC", i, kmer.length(), matchQuality_RC));//  + kmer_rev_len));
            }
        }
        Collections.sort(matches);
    }

    public void print()
    {
        printl(genome);
    }

    class Match implements Comparable 
    {
        public String type;
        public int index;
        public int length;
        public int exactness;

        public Match(String type, int index, int length, int exactness)
        {
            this.type = type; 
            this.index = index;
            this.length = length;
            this.exactness = exactness;
        }

        public String toString()
        {
            return type + " @ " + index + " Q: " + exactness;
        }

        public String toStringDetailed()
        {
            StringBuilder str = new StringBuilder();
            int start = Math.max( 0, index - 3);
            int end = Math.min( index + length + 3 , genome.length() - 1);
            str.append(type + " @ " + index);
            str.append(" Q: " + exactness );
            str.append(" ..." + genome.substring( start, index));
            str.append("<" + genome.substring(index, index + length));
            str.append(">" + genome.substring(index + length, end));
            str.append("...");
            return str.toString();
        }

        public int compareTo(Object o)
        {
            Match match = (Match)o;

            if(matchDirSorting(type) < matchDirSorting(match.type))
            {
                return -1;
            }
            else if(matchDirSorting(type) > matchDirSorting(match.type))
            {
                return 1;
            }
            else
            {
                if( index < match.index ) 
                    return -1;
                else
                    return 1;
            }
        }

        int matchDirSorting(String dir)
        {
            switch(dir)
            {
                case "F": return 0;
                
                case "RC": return 1;
                
                default: return -1;
            }
        }
    }
}