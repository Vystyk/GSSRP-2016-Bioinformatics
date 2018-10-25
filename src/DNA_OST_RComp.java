// -----------------------------------------------
// ---
// ---         GSSRP 2016 Bioinformatics
// ---             Victor Jacobson
// ---                2016.07.14
// ---
// ---              Using Strings
// ---            Major Optimization:
// ---          Instead of adding complete substrings 
// ---          of each suffix in nodes, save only the
// ---          start and end indexes and reference 
// ---          the genome string directly, reducing 
// ---          suffix tree size from (n^2)/2 to 
// ---          approx. n*10 to n*20
// ---          
// ---          Added: Search for reverse and compliment matches
// ---              
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class DNA_OST_RComp
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        printl("Hello!");

        // Initialize the suffix tree with a FASTA file
        //String file_name = "genomes/Zika.NC_012532.1.fasta"; // 11k
        String file_name = "genomes/EcoliK12DH10B.fasta"; // 4.7m

        DNA_OST_RComp ecoliSuffixTree = new DNA_OST_RComp(file_name);

        // Print to file
        //zikaSuffixTree.print();
        ecoliSuffixTree.saveSuffixTreeToFile("genomes/EcoliOptimizedSuffixTree.dat");

        // TODO Test search
        String kmer = "ATCAACAGGTTTAATTTGGATTT";//"GGGGTTTTT";

        ecoliSuffixTree.searchGenomeRComp(kmer);

        ecoliSuffixTree.printResults();

    }

    // Data fields
    public String genome;
    public Node root = new Node();

    public String lastSearch = "";
    public ArrayList<Match> matches;

    public DNA_OST_RComp() { }

    public DNA_OST_RComp(String file_name)
    {

        printl("Scanning file: " + file_name + " ...");
        genome = parseSequenceFromFile(file_name);

        printl("Generating Suffix Tree...");

        generateSuffixTree();

    }

    private void generateSuffixTree()
    {
        System.gc();
        
        root = new Node();
        Node currentNode;

        // Start adding suffixes from smallest to largest
        for(int i = genome.length() - 1; i >= 0; i--)
        {
            // Get the current suffix
            //String suffix = genome.substring(i);
            int suffix_start = i;

            // Start at the Root
            currentNode = root;
            int nodeSearchIndex = 0;
            boolean newNodeCreated = false;
            int last_index = genome.length() - 1;

            while( !newNodeCreated )
            {
                // Does the current node string contain a char to compare to the suffix
                if( nodeSearchIndex >= currentNode.length() )// No more chars to compare (should never be greater than)
                {
                    // Remove the portion of the suffix that is already contained in the current node
                    suffix_start += nodeSearchIndex;
                    nodeSearchIndex = 0;

                    // TODO account for aditional codes R, Y, etc...
                    //                     if( charToIndex(genome.charAt(suffix_start)) == -1)
                    //                     {
                    //                         printl("Error, File contains: " + genome.charAt(suffix_start));
                    //                         return;
                    //                     }
                    //-----------------------------------------------------------------

                    // Is there a child node for the first char of the remaining suffix
                    if( currentNode.children[charToIndex(genome.charAt(suffix_start))] == null )// No corresponding child node
                    {
                        // Create new node containing the suffix of the suffix
                        currentNode.children[charToIndex(genome.charAt(suffix_start))] = new Node(suffix_start, last_index, i);
                        newNodeCreated = true;
                    }
                    else// There is a node
                    {
                        // Go to the next node
                        currentNode = currentNode.children[charToIndex(genome.charAt(suffix_start))];
                    }
                }
                else // The node has more characters to compare
                {
                    // Compare the corresponding characters of the suffix and the node string
                    if ( genome.charAt(suffix_start + nodeSearchIndex ) == genome.charAt(currentNode.start + nodeSearchIndex) )// Characters match
                    {
                        // Itterate to next position
                        nodeSearchIndex++;
                    }
                    else // Characters do not match
                    {
                        // Split the nodes
                        // Create a new node with the second half and transfer the children.
                        Node newNode = new Node(currentNode.start + nodeSearchIndex, currentNode.end, currentNode.index);
                        newNode.children = currentNode.children;

                        // Set the current node's string to the first half
                        currentNode.end = newNode.start;
                        currentNode.index = -1;

                        // Reset the current node's children and set the new node to the corresponding child node
                        currentNode.children = new Node[4];
                        currentNode.children[charToIndex(genome.charAt(newNode.start))] = newNode;

                        //Add remainder of s to the corresponding child node
                        currentNode.children[charToIndex(genome.charAt(suffix_start + nodeSearchIndex))] = new Node(suffix_start + nodeSearchIndex , last_index, i);
                        newNodeCreated = true;
                    }
                }
            }// Node checking

        }// Suffix
        System.gc();
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
        // If the kmer is longer than the genome, truncate the kmer
        int kmer_length = kmer.length();
        if( kmer_length > genome.length() )
        {
            kmer = kmer.substring( 0, genome.length());
        }

        // Search through suffix tree
        int searchIndex = 0;
        Node currentNode = root;

        while( true )
        {
            // Is there a char to compare in the kmer
            if( searchIndex < kmer.length() )
            {
                // Does the current node contain a char to compare to the kmer
                if( searchIndex >= currentNode.length())// Reached end of node
                {
                    kmer = kmer.substring(searchIndex);
                    searchIndex = 0;

                    // Are there more nodes to search
                    if( currentNode.children[charToIndex(kmer.charAt(0))] == null )// No more nodes to search, match found
                    {
                        // Collect all the indexes of all child nodes

                        // TODO Partial matchs found
                        return;
                    }
                    else// More nodes to search
                    {
                        currentNode = currentNode.children[charToIndex(kmer.charAt(0))];
                    }
                }
                else // The suffix tree has more characters to compare
                {
                    // Compare the corresponding characters of the suffix and the node string
                    if ( kmer.charAt(searchIndex) == genome.charAt(currentNode.start + searchIndex) )// Characters match
                    {
                        // Itterate to next position
                        searchIndex++;
                    }
                    else// No matches
                    {
                        return;
                    }
                }
            }
            else // Reached end of kmer with no conflicts
            {
                // Collect all the indexes of all child nodes
                matches = getLeafIndexes(currentNode, kmer_length, "F");
            }
        }// Searching

        //return matchIndexes;
    }

    public void searchGenomeRComp(String kmer)
    {
        lastSearch = kmer;

        matches = new ArrayList<Match>();

        matches.addAll( searchGenome(kmer, "F"));
        //matches.addAll( searchGenome(((new StringBuilder(kmer)).reverse()).toString(), "R"));
        kmer = compliment( (new StringBuilder(kmer)).reverse().toString());
        //matches.addAll( searchGenome(kmer, "CF"));
        matches.addAll( searchGenome(kmer, "CR"));

        Collections.sort(matches);
    }

    private ArrayList<Match> searchGenome(String kmer, String searchDir)
    {
        // If the kmer is longer than the genome, truncate the kmer
        int kmer_length = kmer.length();
        if( kmer_length > genome.length() )
        {
            kmer = kmer.substring( 0, genome.length());
        }

        // Search through suffix tree
        int searchIndex = 0;
        Node currentNode = root;

        while( true )
        {
            // Is there a char to compare in the kmer
            if( searchIndex < kmer.length() )
            {
                // Does the current node contain a char to compare to the kmer
                if( searchIndex >= currentNode.length())// Reached end of node
                {
                    kmer = kmer.substring(searchIndex);
                    searchIndex = 0;

                    // Are there more nodes to search
                    if( currentNode.children[charToIndex(kmer.charAt(0))] == null )// No more nodes to search, match found
                    {
                        // Collect all the indexes of all child nodes

                        // TODO Partial matchs found
                        return new ArrayList<Match>();
                    }
                    else// More nodes to search
                    {
                        currentNode = currentNode.children[charToIndex(kmer.charAt(0))];
                    }
                }
                else // The suffix tree has more characters to compare
                {
                    // Compare the corresponding characters of the suffix and the node string
                    if ( kmer.charAt(searchIndex) == genome.charAt(currentNode.start + searchIndex) )// Characters match
                    {
                        // Itterate to next position
                        searchIndex++;
                    }
                    else// No matches
                    {
                        return new ArrayList<Match>();
                    }
                }
            }
            else // Reached end of kmer with no conflicts
            {
                // Collect all the indexes of all child nodes
                return getLeafIndexes(currentNode, kmer_length, searchDir);
            }
        }// Searching

        //return matchIndexes;
    }

    public String compliment(String kmer)
    {
        StringBuilder comp = new StringBuilder();
        for(int i = 0; i < kmer.length(); i++)
        {
            comp.append(compChar(kmer.charAt(i)));
        }
        return comp.toString();
    }

    public char compChar(char n)
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
        }
        return '$';
    }

    private ArrayList<Match> getLeafIndexes(Node node, int kmer_length, String searchDir)
    {
        ArrayList<Match> matches = new ArrayList<Match>();

        LinkedList<Node> nodeQueue = new LinkedList<Node>();
        nodeQueue.add(node);

        // Itterate through all leaf nodes of the node where the match was confirmed
        while(nodeQueue.peek() != null)
        {
            node = nodeQueue.poll();
            for(Node n: node.children)
            {
                if(n != null)
                {
                    nodeQueue.add(n);
                }
            }
            if(node.index >= 0)
            {
                matches.add(new Match(searchDir, node.index, kmer_length));
            }

        }
        Collections.sort(matches);
        return matches;
    }

    private String parseSequenceFromFile(String file_name)
    {
        StringBuilder genome = new StringBuilder();

        File data = new File(file_name);

        try
        {
            Scanner in = new Scanner(data);

            line();
            printl(in.nextLine());

            while(in.hasNext())
            {
                genome.append(in.nextLine());
            }

        }
        catch(FileNotFoundException fnfe)
        {
            fnfe.printStackTrace();
        }

        return genome.toString();
    }

    public void print()
    {
        LinkedList<Object> printQueue = new LinkedList<Object>();

        printQueue.add(root);

        Node nodeToPrint = null;

        while(printQueue.peek() != null)
        {
            if(printQueue.peek() instanceof Node)
            {
                nodeToPrint = (Node)printQueue.poll();
                print(nodeToPrint.toString() + " ");
                for(Node n: nodeToPrint.children)
                {
                    if(n != null)
                    {
                        printQueue.add(n);
                    }
                }
                printQueue.add("/n");
            }
            else
            {
                String s = (String)printQueue.poll();
                line();

            }
        }
    }

    public void saveSuffixTreeToFile(String file_name)
    {
        //FileOutputStream file = new FileOutputStream(file_name);
        // Create file if it does not exist yet
        //PrintWriter writer;

        try(
            OutputStream file = new FileOutputStream(file_name);
            OutputStream buffer = new BufferedOutputStream(file);
            ObjectOutput writer = new ObjectOutputStream(buffer);
            )
        {

            printl("Writing file...");
            
            LinkedList<Object> writeQueue = new LinkedList<Object>();
            writeQueue.add(root);
            Node nodeToPrint = null;
            while(writeQueue.peek() != null)
            {
                nodeToPrint = (Node)writeQueue.poll();
                writer.writeObject(nodeToPrint);
                for(Node n: nodeToPrint.children)
                {
                    if(n != null)
                    {
                        writeQueue.add(n);
                    }
                }
            }

            printl("File Created: " + file_name);
        }
        catch(IOException e)
        {
            printl("Output failed. " + e.toString());
        }
        catch(Exception e)
        {
            printl("file creation failed");
        }

    }

    public static int charToIndex(char c)
    {
        switch(c)
        {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
        }
        // TODO add more cases for aditional codes 
        return 0;
    }

    // Helpers
    public static void line(){System.out.println(); }

    public static void printl(Object o){System.out.println(o); }

    public static void print(Object o){System.out.print(o); }

    class Node// implements Serializable
    {
        //public String string = "";
        public int start, end;
        //public boolean segmentIndex;

        public int index = -1;

        public Node[] children = new Node[4];

        public Node(){}// Root node has empty string

        public Node(int start_index, int end_index, int index)
        {
            this.start = start_index;
            this.end = end_index;
            this.index = index;
        }

        public int length()
        {
            return end - start;
        }

        public String toString()
        {
            //return string + "$" + index;
            return start + "," + end + "$" + index;
        }

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
    }
}