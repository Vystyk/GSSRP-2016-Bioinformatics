// -----------------------------------------------
// ---
// ---             Victor Jacobson
// ---            
// ---   For System.nanoTime(), use double, not float
// ---               
// -----------------------------------------------
import java.util.*;
import java.io.*;

public class NanoTimer
{
    public static void main(String [] args)
    {
        Scanner in = new Scanner(System.in);
        printl("Hello!");

        NanoTimer timer = new NanoTimer();

        timer.start();
        double sum = 0;
        for(int i = 0; i < 10000; i++)
        {
            for(int j = 0; j < 10000; j++)
                sum += Math.sqrt(i);
        }

        timer.stop();

        printl("Time elapsed: " + timer.elapsed());

    }
    private double timeElapsed;
    private double startTime;

    public NanoTimer()
    {

    }

    public void reset()
    {
        timeElapsed = 0;
    }

    public void start()
    {
        startTime = System.nanoTime() + timeElapsed;
        //printl(startTime);
    }

    public void stop()
    {
        //printl(System.nanoTime());
        timeElapsed = System.nanoTime() - startTime;
        //printl(timeElapsed);
    }

    public double elapsed()
    {
        return timeElapsed / 1000000d;
    }
    
    public String elapsedFormated()
    {
        if( timeElapsed > 1000000000d )
            return timeElapsed / 1000000000d + "sec.";
        else if( timeElapsed > 1000000d )
            return timeElapsed / 1000000d + "ms.";
        else if( timeElapsed > 1000d )
            return timeElapsed / 1000d + "μs.";
        else  
            return timeElapsed + "ns.";
            
    }

    public static void line(){System.out.println(); }
    public static void printl(Object o){System.out.println(o); }
    public static void print(Object o){System.out.print(o); }
    public static int randomInt(int max){return (int)(Math.random() * max);}
}