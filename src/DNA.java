// -----------------------------------------------
// ---
// ---             Victor Jacobson
// ---            
// ---             
// ---               
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class DNA
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        printl("Hello!");

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
            case 'S':
            return 'W';
            case 'W':
            return 'S';
            case 'K':
            return 'M';
            case 'M':
            return 'K';

            case 'B':
            return 'V';
            case 'V':
            return 'B';
            case 'D':
            return 'H';
            case 'H':
            return 'D';
        }
        return 'N';
    }

    public static byte compare(char kmer, char genome)
    {
        switch(kmer)
        {
            case 'A': switch(genome)
            {
                case 'A': return 0;// Exact match

                case 'T':
                case 'C':
                case 'G': return -1;// No match

                case 'R':
                case 'W':
                case 'M': return 1;// Wildcard match

                case 'D':
                case 'H':
                case 'V': return 2;

                case 'N': return 3;
            }
            case 'T': switch(genome)
            {
                case 'T': return 0;// Exact match

                case 'A':
                case 'C':
                case 'G': return -1;// No match

                case 'Y':
                case 'W':
                case 'K': return 1;// Wildcard match

                case 'B':
                case 'D':
                case 'H': return 2;

                case 'N': return 3;
            }
            case 'C': switch(genome)
            {
                case 'C': return 0;// Exact match

                case 'A':
                case 'T':
                case 'G': return -1;// No match

                case 'Y':
                case 'S':
                case 'M': return 1;// Wildcard match

                case 'B':
                case 'H':
                case 'V': return 2;

                case 'N': return 3;
            }
            case 'G': switch(genome)
            {
                case 'G': return 0;// Exact match

                case 'A':
                case 'T':
                case 'C': return -1;// No match

                case 'R':
                case 'S':
                case 'K': return 1;// Wildcard match

                case 'B':
                case 'D':
                case 'V': return 2;

                case 'N': return 3;
            }
            default:
            return -1;
        }
    }

    public static void line(){System.out.println(); }

    public static void printl(Object o){System.out.println(o); }

    public static void print(Object o){System.out.print(o); }

    public static int randomInt(int max){return (int)(Math.random() * max);}
}