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

public class DNA_NSS_RComp
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        print("DNA_NSS_RComp! - Reverse and complimentary searches");

        DNA_NSS_RComp data = new DNA_NSS_RComp();

        data.newGenomeFromFASTA("genomes/Zika.NC_012532.1.fasta");

        data.print();

        printl("F = forard, R = reverse, C = Compliment");

        String kmer = "GAAAAA";

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

    public DNA_NSS_RComp() { }

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
            print(matches.get(i).toString() + ", ");
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
            //             boolean match_r = true;
            //             boolean match_cf = true;
            boolean match_cr = true;

            int kmer_rev_len = kmer.length() - 1;

            // if it matches, check each subsequent character for a match
            for(int j = 0; j < kmer.length(); j++)
            {
                // Find forward match
                if(genome.charAt(i + j) != kmer.charAt(j))
                {
                    match_f = false;
                    //continue;
                }

                // Find complimentary reverse match
                if(genome.charAt(i + j) != comp(kmer.charAt(kmer_rev_len - j)))
                {
                    match_cr = false;
                    //continue;
                }

                if( !match_f && !match_cr )
                    continue;

            }
            // if we get to the end of the kmer and have a match, record the index of the match
            if(match_f)
            {
                matches.add(new Match("F", i, kmer.length()));
            }

            if(match_cr)
            {
                matches.add(new Match("CR", i, kmer.length()));//  + kmer_rev_len));
            }
        }
        Collections.sort(matches);
    }

    public static char comp(char n)
    {
        switch(n)
        {
            case 'A':
            return 'T';
            case 'T':
            return 'A';
            case 'C':
            return 'G';
            case 'G':
            return 'C';
            case 'R':
            return 'Y';
            case 'Y':
            return 'R';
            
        }
        return '$';
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

        public Match(String type, int index, int length)
        {
            this.type = type; 
            this.index = index;
            this.length = length;
        }

        public String toString()
        {
            return type + " @ " + index;
        }

        public String toStringDetailed()
        {
            StringBuilder str = new StringBuilder();
            int start = Math.max( 0, index - 3);
            int end = Math.min( index + length + 3 , genome.length() - 1);
            str.append(type + " @ " + index);
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
                case "F":
                return 0;
                case "R":
                return 1;
                case "CF":
                return 2;
                case "CR":
                return 3;
                default:
                return -1;
            }
        }
    }
}