// -----------------------------------------------
// ---
// ---             Victor Jacobson
// ---            
// ---       Generates random kmers for 
// ---          testing algorithms
// ---               
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class KmerGenerator
{
    // -------- Main method for testing ----------
    public static void main(String [] args) throws Exception
    {
        Scanner in = new Scanner(System.in);
        print("Hello!");
        
        KmerGenerator kmers = new KmerGenerator(100);
        
        kmers.printAll();
        
    }
    
    
    // ------------- Data fields -----------------
    
    public ArrayList<String> kmers = new ArrayList<String>(); 
    
    
    public KmerGenerator(){ }
    
    public KmerGenerator(int init_size)
    {
        for(int i = 0; i < init_size; i++)
        {
            kmers.add(generateKmer());
        }
    }
    
    public String generateKmer()
    {
        StringBuilder kmer = new StringBuilder();
        
        for(int i = 0; i < 30; i++)
        {
            int nucleotide = randomInt(4);
            switch(nucleotide)
            {
                case 0: kmer.append('A');break;
                case 1: kmer.append('C');break;
                case 2: kmer.append('G');break;
                case 3: kmer.append('T');break;
            }
        }
        
        return kmer.toString();
    }
    
    public String getKmer(int length)
    {
        StringBuilder kmer = new StringBuilder();
        
        for(int i = 0; i < length; i++)
        {
            int nucleotide = randomInt(4);
            switch(nucleotide)
            {
                case 0: kmer.append('A');break;
                case 1: kmer.append('C');break;
                case 2: kmer.append('G');break;
                case 3: kmer.append('T');break;
            }
        }
        
        return kmer.toString();
    }
    
    
    public String get(int i)
    {
        return kmers.get(i);
    }
    
    public void printAll()
    {
        for(String kmer : kmers)
        {
            print(kmer);
        }
        
    }
    
    
    public static void line(){System.out.println(); }
    public static void print(Object o){System.out.println(o); }
    public static int randomInt(int max){return (int)(Math.random() * max);}
}