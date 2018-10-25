// -----------------------------------------------
// ---
// ---         GSSRP 2016 Bioinformatics
// ---             Victor Jacobson
// ---            
// ---              Not optimized
// ---           (using chars/strings)
// ---          
// ---      FixedBuggs:
// ---      2016.07.14  -Fixed incorrectly returning partial matches in some cases
// ---                  -Sort match indexes
// ---  
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class DNA_NaiveSuffixTreeSearch
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        printl("Hello!");

        // Initialize the suffix tree with a FASTA file
        //String file_name = "genomes/EcoliK12DH10B.fasta";
        String file_name = "genomes/Zika.NC_012532.1.fasta";
        DNA_NaiveSuffixTreeSearch zikaSuffixTree = new DNA_NaiveSuffixTreeSearch(file_name);

        // Print to file
        //zikaSuffixTree.print();
        //zikaSuffixTree.saveSuffixTreeToFile("genomes/EcoliSuffixTree.dat");

        // TODO Test search
        String stringToSearch = "AAAAA";
        ArrayList<Integer> matches = zikaSuffixTree.searchGenome(stringToSearch);

        for(int i = 0; i < matches.size(); i++)
        {
            printl(matches.get(i));
        }
        printl("There were " + matches.size() + " matches for " + stringToSearch);

    }

    // Data fields
    public String genome;
    public Node root = new Node();

    public DNA_NaiveSuffixTreeSearch() { }

    public DNA_NaiveSuffixTreeSearch(String file_name)
    {
        genome = parseSequenceFromFile(file_name);

        //printl( genome.substring(0));

        generateSuffixTree();

    }

    public void generateSuffixTree()
    {

        root = new Node();
        Node currentNode;

        // Start adding suffixes from smallest to largest
        for(int i = genome.length() - 1; i >= 0; i--)
        {
            // Ger the current suffix
            String suffix = genome.substring(i);

            // Start at the Root
            currentNode = root;
            int nodeStringIndex = 0;
            boolean newNodeCreated = false;

            while( !newNodeCreated )
            {
                // Does the current node string contain a char to compare to the suffix
                if( nodeStringIndex >= currentNode.string.length())// No more chars to compare (should never be greater than)
                {
                    // Remove the portion of the suffix that is already contained in the current node
                    suffix = suffix.substring(nodeStringIndex);
                    nodeStringIndex = 0;
                    // Is there a child node for the first char of the remaining suffix 
                    if( currentNode.children[charToIndex(suffix.charAt(0))] == null )// No corresponding child node
                    {
                        // Create new node containing the suffix of the suffix
                        currentNode.children[charToIndex(suffix.charAt(0))] = new Node(suffix, i);
                        newNodeCreated = true; 
                    }
                    else// There is a node
                    {
                        // Go to the next node
                        currentNode = currentNode.children[charToIndex(suffix.charAt(0))];
                    }
                }
                else // The node has more characters to compare
                {
                    // Compare the corresponding characters of the suffix and the node string
                    if ( suffix.charAt(nodeStringIndex) == currentNode.string.charAt(nodeStringIndex) )// Characters match
                    {
                        // Itterate to next position
                        nodeStringIndex++;
                    }
                    else // Characters do not match
                    {
                        // Split the nodes
                        String stringToSplit = currentNode.string;
                        String firstHalf = stringToSplit.substring(0, nodeStringIndex);
                        String secondHalf = stringToSplit.substring(nodeStringIndex);

                        // Create a new node with the second half and transfer the children.
                        Node newNode = new Node(secondHalf, currentNode.index);
                        newNode.children = currentNode.children;

                        // Set the current node's string to the first half
                        currentNode.string = firstHalf;
                        currentNode.index = -1;

                        // Reset the current node's children and set the new node to the corresponding child node
                        currentNode.children = new Node[4];
                        currentNode.children[charToIndex(secondHalf.charAt(0))] = newNode;

                        //Add remainder of s to the corresponding child node
                        suffix = suffix.substring(nodeStringIndex);
                        currentNode.children[charToIndex(suffix.charAt(0))] = new Node(suffix, i);
                        newNodeCreated = true;
                    }
                }
            }// Node checking

        }// Suffix

    }

    public ArrayList<Integer> searchGenome(String kmer)
    {
        // If the kmer is longer than the genome, truncate the kmer
        if( kmer.length() > genome.length() )
        {
            kmer = kmer.substring( 0, genome.length());
        }

        // TODO search through suffix tree
        int searchIndex = 0;
        Node currentNode = root;

        while( true )
        {
            // Is there a char to compare in the kmer
            if( searchIndex < kmer.length() )
            {
                // Does the current node contain a char to compare to the kmer
                if( searchIndex >= currentNode.string.length())// Reached end of node
                {
                    kmer = kmer.substring(searchIndex);
                    searchIndex = 0;

                    // Are there more nodes to search
                    if( currentNode.children[charToIndex(kmer.charAt(0))] == null )// No more nodes to search, match found
                    {
                        // Collect all the indexes of all child nodes

                        //TODO Partial matches found
                        return new ArrayList<Integer>();
                    }
                    else// More nodes to search
                    {
                        currentNode = currentNode.children[charToIndex(kmer.charAt(0))];
                    }
                }
                else // The suffix tree has more characters to compare
                {
                    // Compare the corresponding characters of the suffix and the node string
                    if ( kmer.charAt(searchIndex) == currentNode.string.charAt(searchIndex) )// Characters match
                    {
                        // Itterate to next position
                        searchIndex++;
                    }
                    else// No matches
                    {
                        return new ArrayList<Integer>();
                    }
                }
            }
            else // Reached end of kmer with no conflicts
            {
                // Collect all the indexes of all child nodes
                return getLeafIndexes(currentNode);
            }
        }// Searching

        //return matchIndexes;
    }

    private ArrayList<Integer> getLeafIndexes(Node node)
    {
        ArrayList<Integer> matchIndexes = new ArrayList<Integer>();

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
                matchIndexes.add(node.index);
            }

        }
        Collections.sort(matchIndexes);
        return matchIndexes;
    }

    public static String parseSequenceFromFile(String file_name)
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

            printl("Processed: " + file_name);
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
        File file = new File(file_name);
        // Create file if it does not exist yet
        PrintWriter writer;

        try
        {
            file.createNewFile();

            writer = new PrintWriter( file);   
        }
        catch(Exception e)
        {
            printl("file creation failed");
            return;
        }

        LinkedList<Object> printQueue = new LinkedList<Object>();
        printQueue.add(root);
        Node nodeToPrint = null;
        while(printQueue.peek() != null)
        {
            nodeToPrint = (Node)printQueue.poll();
            writer.println(nodeToPrint.toString());
            for(Node n: nodeToPrint.children)
            {
                if(n != null)
                {
                    printQueue.add(n);
                }
            }

            printl("File Created: " + file_name);
            writer.close();
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
            //case '$': return 4;
        }
        return 0;
    }

    // Helpers
    public static void line(){System.out.println(); }

    public static void printl(Object o){System.out.println(o); }

    public static void print(Object o){System.out.print(o); }

    class Node
    {
        public String string = "";

        public int index = -1;

        public Node[] children = new Node[4];

        public Node(){}// Root node has empty string

        public Node(String string, int index)
        {
            this.string = string;
            this.index = index;
        }

        public String toString()
        {
            return string + "$" + index;
        }

    }
}
