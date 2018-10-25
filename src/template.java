// -----------------------------------------------
// ---
// ---             Victor Jacobson
// ---            
// ---             
// ---               
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class template
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        printl("Hello!");
        
        
    }
    public static void line(){System.out.println(); }
    public static void printl(Object o){System.out.println(o); }
    public static void print(Object o){System.out.print(o); }
    public static int randomInt(int max){return (int)(Math.random() * max);}
}