package com.cs.sjsu.edu.HMMPayl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Random;

public class HMM {

	private static final int MAX_CHARS = 0;
	public static double pi[];
	public static double piBar[];
	public static double A[][];
	public static double Abar[][];
	public static double B[][];
	public static double Bbar[][];
	
	public static StepStruct step[];
	
	public static void main(String[] args) 
	{
	
		int N = 2 ;   // Number of states 
		int M = 27;   // Number of observations
		int seed;
		
		
		int startPos,startChar,maxChars,maxIters,i,j,T,iter;
		
		startPos = Integer.parseInt(args[1]);
		startChar = Integer.parseInt(args[2]);			// The character from which to start saving observations for training .
		maxChars =Integer.parseInt(args[3]);			// maximum number of characters to be saved .
		maxIters = Integer.parseInt(args[4]);
		seed     = Integer.parseInt(args[5]);
		

	    double logProb,
	           newLogProb;

	    // The matrices of HMM
	     pi = new double [N];
	    piBar = new double[N];
	    A= new double [N][N] ;
	    Abar = new double [N][N];
	    B = new double[N][M];
	    Bbar = new double[N][M];
	    
	  
	    
	    // input file
	    String fname = new String();
	    fname = args[0];
	    
	    // Test for command line arguments 
	    //System.out.println("The input parameters are " + fname + " " + startPos + " " + startChar + " " + maxChars + " " + maxIters + " "+ seed );
	    
	    
	    ////////////////////////
	    // read the data file //
	     ////////////////////////

	    // determine number of observations
	    System.out.println("GetT... ");
	    T = GetT(fname,startPos,startChar,maxChars);
	    System.out.println("Number of Observations T = "+T);
	    
	    
	    /////////////////////////////////////////////////////////////////
	    //Determine the number of Steps ie total number of observations//
	    /////////////////////////////////////////////////////////////////
	    
	    step =new StepStruct[T];
	    
	    
	    //////////////////////////
	    //Read the observations //
	    //////////////////////////
	    
	    System.out.println("GetObservations...");
	    int obs = GetObservations(fname,T,startPos,startChar,maxChars);
	    
	    System.out.println("Number of Observations "+obs);
	    
	    /* Print the contents of step observation on console . Verify whether step is populated correctly or not 
	    for(i=0;i<step.length;i++)
	    {
	    	System.out.println("The observation "+step[i].obs+" and the number "+i);
	    }
	    */
	    
	    ///////////////////////
	    //Hidden Markov Model//
	    //////////////////////
	    
	    initMatrices(seed);
	    
	    System.out.println("The Pi Matrix");
	    printPi();
	    System.out.println("\nThe A Matrix");
	    printA();
	    System.out.println("\nThe B Matrix");
	    printB();
	    
	    // initialization
	    iter = 0;
	    logProb = -1.0;
	    newLogProb = 0.0;
	    
	    
	    // main loop
	    int count =1;
	    
	    // Do 100 iterations during the training phase of HMM .
	    while (iter <maxIters)
	    {
	    	System.out.println("\nBegin iteration"+iter);
	    	
	    	logProb = newLogProb;
	    	
	    	System.out.print("\nAlpha pass..");
	    	alphaPass(T);
	    	System.out.println("Done");
	    	
	    	System.out.print("\nBeta pass..");
	    	betaPass(T);
	    	System.out.println("Done");
	    	
	    	System.out.println("\nCompute Gammas and DiGammas");
	    	computeGammas(T);
	    	System.out.println("Done");
	    	
	    	System.out.println("\nReestimate Pi");
	    	reEstimatePi();
	    	System.out.println("Done");
	    	
	    	System.out.println("\nRestimate A");
	    	reEstimateA(T);
	    	System.out.println("Done");
	    	
	    	System.out.println("\nRestimate B");
	    	reEstimateB(T);
	    	System.out.println("Done");
	    	
	    	// assign pi, A and B corresponding "new" values
	        for(i = 0; i < N; ++i)
	        {
	            pi[i] = piBar[i];
	        
	            for(j = 0; j < N; ++j)
	            {
	                A[i][j] = Abar[i][j];
	            }

	            for(j = 0; j < M; ++j)
	            {
	                B[i][j] = Bbar[i][j];
	            }
	            
	        }// next i
	        
	        // Print out the new values of the matrices 
	        System.out.println("New Pi");
	        printPi();
	        System.out.println("New A");
	        printA();
	        System.out.println("New B");
	        printB();
	        
	        System.out.println();
	        
	        // compute log [P(observations | lambda)], where lambda = (A,B,pi)
	        newLogProb = 0.0;
	        for(i = 0; i < T; ++i)
	        {
	            newLogProb += Math.log(step[i].c);
	        }
	        newLogProb = -newLogProb;

	        // a little trick so that no initial logProb is required
	        if(iter == 0)
	        {
	            logProb = newLogProb - 1.0;
	        }
	    	
	        System.out.println("completed iteration = " + iter +" log [P(observation | lambda)]"+newLogProb);
	    	
	        // Increment the iterations
	        iter ++; 
	    	
	    }
	    
	    
	    
	    
	}

	//
	// read (but don't save) observations get T
	//
	public static int GetT(String fname,
	         int startPos,
	         int startChar,
	         int maxChars)
	{
	    // FILE *in;
	  
		
		try
		{
		File file = new File(fname);
	
	    BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));

	    int i,j,len,thisStartPos,totalNum,num;
	        
	    char [] temp = new char [MAX_CHARS + 1];
	    
	    char space[] = {' '};
	    
	    // List of observations
	    char alphabetLower[] = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',' '};
	    
	    char alphabetUpper[] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',' '};
	    
	    // Number of Observations
	    int M =27 ; 
	    int N=2 ;

	    // count 'em
	    totalNum = num = 0;
	    
	    String line = new String();
	    
	    // Count the lines to add 'space' characters at end of each line 
	    int countlines =0 ;
	    
	    while((line=reader.readLine())!=null)
	    {
	        // Convert String 'line' to a character array 'temp'
	        temp = line.toCharArray();
	        
	        //Increment line count
	        countlines ++;
	        
	        
	        len = temp.length;

	        // each line should end with a single space
	        
	        while(temp[len -1]==space[0] && len>0)
	        {
	        	len --;
	        }
	      
	       /* while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
	        {
	            --len;
	        }*/

	        //System.out.println("The size of temp is "+temp.length);
	        //System.out.println("The length of the line is "+len);
	        
	        // Printing the line character by character to the console
	        
	        //for(i=0;i<temp.length;i++)
	        	//System.out.print(temp[i]);
	        //System.out.println();
	        
	        
	        //Append a single space character at the end of the line . This solution does not work here .
	        
	        //temp[len]=space[0];
	        //strncpy(&temp[len], space, 1);
	        
	        thisStartPos = startPos;
	        
	        // ignore leading spaces
	        while(temp[thisStartPos]==space[0])
	        {
	        	thisStartPos++;
	        }
	        
	        /*while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
	        {
	            ++thisStartPos;
	        }*/
	        
	        for(i = thisStartPos; i < len; ++i)
	        {
	            // find alphabetic characters, ignoring case
	            // also drop all non-alphabet characters other than space
	        	
	            for(j = 0; j < M; ++j)
	            {
	                if(temp[i]==alphabetLower[j] || temp[i]==alphabetUpper[j])
	                {
	                    ++totalNum;
	                    if(totalNum >= startChar)
	                    {
	
	                        //printf("%c %d\n", alphabet[j], num);
	                       
	                    	++num;
	                    	
	                        if((maxChars > 0) && (num >= maxChars))
	                        {
	                        	// return the counted characters + number of lines = to add exact number of spaces to the observation sequence .
	                        	// Temperory fix - No more !!
	                            return(num);
	                        }
	                        
	                    }// end if
	                    
	                    break;
	                    
	                }// end if
	                
	            }// next j

	        }// next i
	        
	        // Add the spaces at end of the lines
	        num++ ;
	    }// end while

	    // Print the total number of lines 
	    System.out.println("Total number of lines is "+countlines);
	    
	    // return the counted characters + number of lines = to add exact number of spaces to the observation sequence .
    	// Temperory fix -No more !!
	    return(num);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		
		return 0 ;
	}// end GetT

	
	
	//
	// read and save in memory observations get T
	//
	public static int GetObservations(String fname,int T, 
	         int startPos,
	         int startChar,
	         int maxChars)
	{
	    // FILE *in;
	  
		
		try
		{
		File file = new File(fname);
	
	    BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));

	    int i,j,len,thisStartPos,totalNum,num;
	        
	    char [] temp = new char [MAX_CHARS + 1];
	    
	    char space[] = {' '};
	    
	    // List of observations
	    char alphabetLower[] = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',' '};
	    
	    char alphabetUpper[] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',' '};
	    
	    // Number of Observations
	    int M =27 ; 
	    int N=2;

	    // count 'em
	    totalNum = num = 0;
	    
	    String line = new String();
	    
	    // Count the lines to add 'space' characters at end of each line 
	    int countlines =0 ;
	    
	    while((line=reader.readLine())!=null)
	    {
	        // Convert String 'line' to a character array 'temp'
	        temp = line.toCharArray();
	        
	        //Increment line count
	        countlines ++;
	        
	        
	        len = temp.length;

	        // each line should end with a single space
	        
	        while(temp[len -1]==space[0] && len>0)
	        {
	        	len --;
	        }
	      
	       /* while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
	        {
	            --len;
	        }*/

	        //System.out.println("The size of temp is "+temp.length);
	        //System.out.println("The length of the line is "+len);
	        
	        // Printing the line character by character to the console
	        
	        //for(i=0;i<temp.length;i++)
	        	//System.out.print(temp[i]);
	        //System.out.println();
	        
	        
	        //Append a single space character at the end of the line . This solution does not work here .
	        
	        //temp[len]=space[0];
	        //strncpy(&temp[len], space, 1);
	        
	        thisStartPos = startPos;
	        
	        // ignore leading spaces
	        while(temp[thisStartPos]==space[0])
	        {
	        	thisStartPos++;
	        }
	        
	        /*while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
	        {
	            ++thisStartPos;
	        }*/
	        
	        for(i = thisStartPos; i < len; ++i)
	        {
	            // find alphabetic characters, ignoring case
	            // also drop all non-alphabet characters other than space
	        	
	            for(j = 0; j < M; ++j)
	            {
	                if(temp[i]==alphabetLower[j] || temp[i]==alphabetUpper[j])
	                {
	                    ++totalNum;
	                    if(totalNum >= startChar)
	                    {
	
	                        //printf("%c %d\n", alphabet[j], num);
	                    	
	                    	//Print the alphabets on the console for verifying whether they are read correctly or not 
	                    	// System.out.println(temp[i]);
	                    	
	                    	// Record the character in the step 
	                    	step[num] = new StepStruct();
	                    	step[num].obs = j;
	                        
	                    	++num;
	                    	
	                    	if(num > T) 
	                    	{
	                    		System.out.println("T exceeded !! No more memory");
	                    		return 0;
	                    	}
	                    	
	                        if((maxChars > 0) && (num >= maxChars))
	                        {
	                        	// return the counted characters + number of lines = to add exact number of spaces to the observation sequence .
	                        	// Temperory fix - No more !!
	                            return(num);
	                        }
	                        
	                    }// end if
	                    
	                    break;
	                    
	                }// end if
	                
	            }// next j

	        }// next i
	        
	        step[num] = new StepStruct();
	        step[num].obs = 26;
	        num ++; // Add the space at the end of the line
	    }// end while

	    // Print the total number of lines 
	    System.out.println("Total number of lines is "+countlines);
	    
	    // return the counted characters + number of lines = to add exact number of spaces to the observation sequence .
    	// Temperory fix - No more !!
	    return(num);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		
		return 0 ;
	}// end GetObservations
	
	
	public static void initMatrices(int seed)
	{
		int N = 2 ;
		int M=27;
		
		 int i,j;
	        
	    double prob,ftemp,ftemp2;
	    
	    Random rand = new Random((long)seed);
	    
	    
	    // initialize pi
	    prob = 1.0 / (double)N;
	    ftemp = prob / 10.0;
	    ftemp2 = 0.0;
	    for(i = 0; i < N; ++i)
	    {
	        if((rand.nextInt() & 0x1) == 0)
	        {
	            pi[i] = prob + (double)(rand.nextInt() & 0x7) / 8.0 * ftemp;
	        }
	        else
	        {
	            pi[i] = prob - (double)(rand.nextInt() & 0x7) / 8.0 * ftemp;
	        }
	        ftemp2 += pi[i];
	        
	    }// next i
	    
	    for(i = 0; i < N; ++i)
	    {
	        pi[i] /= ftemp2;
	    }
		
		
	 // initialize A[][]
	    prob = 1.0 / (double)N;
	    ftemp = prob / 10.0;
	    for(i = 0; i < N; ++i)
	    {
	        ftemp2 = 0.0;
	        for(j = 0; j < N; ++j)
	        {
	            if((rand.nextInt() & 0x1) == 0)
	            {
	                A[i][j] = prob + (double)(rand.nextInt() & 0x7) / 8.0 * ftemp;
	            }
	            else
	            {
	                A[i][j] = prob - (double)(rand.nextInt() & 0x7) / 8.0 * ftemp;
	            }
	            ftemp2 += A[i][j];
	            
	        }// next j
	        
	        for(j = 0; j < N; ++j)
	        {
	            A[i][j] /= ftemp2;
	        }
	        
	    }// next i
	    
	    
	 // initialize B[][]
	    prob = 1.0 / (double)M;
	    ftemp = prob / 10.0;
	    for(i = 0; i < N; ++i)
	    {
	        ftemp2 = 0.0;
	        for(j = 0; j < M; ++j)
	        {
	            if((rand.nextInt() & 0x1) == 0)
	            {
	                B[i][j] = prob + (double)(rand.nextInt() & 0x7) / 8.0 * ftemp;
	            }
	            else
	            {
	                B[i][j] = prob - (double)(rand.nextInt() & 0x7) / 8.0 * ftemp;
	            }
	            ftemp2 += B[i][j];
	            
	        }// next j
	        
	        for(j = 0; j < M; ++j)
	        {
	            B[i][j] /= ftemp2;
	        }
	        
	    }// next i
	    
	} // end initMatrices()

	
	public static void printPi()
	{
		int i;
		int N=2;
		int M=27;
        
	    double ftemp;

	    ftemp = 0.0;
	    for(i = 0; i < N; ++i)
	    {
	        System.out.print(pi[i] + " ");
	        ftemp += pi[i];
	    }
	    System.out.println(" sum = "+ftemp);
	}//end printPi()
	
	public static void printA()
	{
		int i,j ;
		int N=2;
		int M=27;
		
		double ftemp;

	    for(i = 0; i < N; ++i)
	    {
	        ftemp = 0.0;
	        for(j = 0; j < N; ++j)
	        {
	            System.out.print(A[i][j]+ " ");
	            ftemp += A[i][j];
	        }
	        System.out.println("sum = "+ ftemp);
	        
	    }// next i
	}// end printA()
	
	public static void printB()
	{
		int i,j ;
		int N=2;
		int M=27;
		
		double ftemp;

	    for(i = 0; i < N; ++i)
	    {
	        ftemp = 0.0;
	        for(j = 0; j < M; ++j)
	        {
	            System.out.print(B[i][j]+ " ");
	            ftemp += B[i][j];
	        }
	        System.out.println("sum = "+ ftemp);
	        
	    }// next i
	}// end printB
	
	public static void alphaPass(int T)
	{
		int i,j,t;
        double ftemp;
        
        int N =2 ; //Number of states
        int M=27;  //Number of observations
        
        // compute alpha[0]'s
	    ftemp = 0.0;
	    for(i = 0; i < N; ++i)
	    {
	        step[0].alpha[i] = pi[i] * B[i][step[0].obs];
	        ftemp += step[0].alpha[i];
	    }
	    step[0].c = 1.0 / ftemp;
	
	    // scale alpha[0]'s
	    for(i = 0; i < N; ++i)
	    {
	        step[0].alpha[i] /= ftemp;
	    }
	
	    // alpha pass
	    for(t = 1; t < T; ++t)
	    {
	        ftemp = 0.0;
	        for(i = 0; i < N; ++i)
	        {
	            step[t].alpha[i] = 0.0;
	            for(j = 0; j < N; ++j)
	            {
	                step[t].alpha[i] += step[t - 1].alpha[j] * A[j][i];
	            }
	            step[t].alpha[i] *= B[i][step[t].obs];
	            ftemp += step[t].alpha[i];
	        }
	        step[t].c = 1.0 / ftemp;
	        
	        // scale alpha's
	        for(i = 0; i < N; ++i)
	        {
	            step[t].alpha[i] /= ftemp;
	        }
	    
	    }// next t
	}// end alphaPass
	
	//
	// beta pass (or backwards pass) including scaling
	//
	public static void betaPass(int T)
	{
		int i,j,t;
		int N =2; //Number of states
		int M=27; // Number of observations

	    // compute scaled beta[T - 1]'s
	    for(i = 0; i < N; ++i)
	    {
	        step[T - 1].beta[i] = 1.0 * step[T - 1].c;
	    }
	
	    // beta pass
	    for(t = T - 2; t >= 0; --t)
	    {
	        for(i = 0; i < N; ++i)
	        {
	            step[t].beta[i] = 0.0;
	            for(j = 0; j < N; ++j)
	            {
	                step[t].beta[i] += A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j];
	            }
	            
	            // scale beta's (same scale factor as alpha's)
	            step[t].beta[i] *= step[t].c;
	        }
	
	    }// next t
			
	}// end betaPass
	
	public static void computeGammas(int T)
	{
		int i,j,t;
		int N=2 ; //Number of states
		int M=27 ;//Number of Observations
        
    double denom;


    double ftemp,
           ftemp2;


	    // compute gamma's and diGamma's
	    for(t = 0; t < T - 1; ++t)
	    {
	        denom = 0.0;
	        for(i = 0; i < N; ++i)
	        {
	            for(j = 0; j < N; ++j)
	            {
	                denom += step[t].alpha[i] * A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j];
	            }
	        }
	            
	
	        ftemp2 = 0.0;
	
	        for(i = 0; i < N; ++i)
	        {
	            step[t].gamma[i] = 0.0;
	            for(j = 0; j < N; ++j)
	            {
	                step[t].diGamma[i][j] = (step[t].alpha[i] * A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j])
	                                        / denom;
	                step[t].gamma[i] += step[t].diGamma[i][j];
	            }
	
	
	            // verify that gamma[i] == alpha[i]*beta[i] / sum(alpha[j]*beta[j])
	            ftemp2 += step[t].gamma[i];
	            ftemp = 0.0;
	            for(j = 0; j < N; ++j)
	            {
	                ftemp += step[t].alpha[j] * step[t].beta[j];
	            }
	            ftemp = (step[t].alpha[i] * step[t].beta[i]) / ftemp;
	            
	            /*if(DABS(ftemp - step[t].gamma[i]) > EPSILON)
	            {
	                System.out.println("gamma[%d] = %f (%f) ", i, step[t].gamma[i], ftemp);
	                System.out.println("********** Error !!!\n");
	            }*/
	
	
	        }// next i
	            
	
	        /*if(DABS(1.0 - ftemp2) > EPSILON)
	        {
	            printf("sum of gamma's = %f (should sum to 1.0)\n", ftemp2);
	        }*/
	
	            
	    }// next t

	}//end computeGammas
	
	public static void reEstimatePi()
	{
		int i;
		int N =2 ; // Number of states
		int M=27 ; // number of observations
	    
	    // reestimate pi[]        
	    for(i = 0; i < N; ++i)
	    {
	        piBar[i] = step[0].gamma[i];
	    }
	} //end reEstimatePi()
	
	public static void reEstimateA(int T)
	{
		int i,j,t;
		int N= 2; //Number of states
		int M= 27 ; //Number of observations
    
		double numer,denom;
           
    // reestimate A[][]
	    for(i = 0; i < N; ++i)
	    {
	        for(j = 0; j < N; ++j)
	        {
	            numer = denom = 0.0;
	
	            // t = 0,1,2,...,T-1
	            for(t = 0; t < T - 1; ++t)
	            {
	                numer += step[t].diGamma[i][j];
	                denom += step[t].gamma[i];
	                
	            }// next t
	
	            Abar[i][j] = numer / denom;
	        
	        }// next j
	        
	    }// next i
	}// end reEstimateA()
	
	public static void reEstimateB(int T)
	{
		int i,j,t;
		
		int N= 2; //Number of states
		int M =27 ; // Number of observations
    
    double numer,denom;
           
    	// reestimate B[][]
    	for(i = 0; i < N; ++i)
    	{
    		for(j = 0; j < M; ++j)
    		{
    			numer = denom = 0.0;
    			
    			// t = 0,1,2,...,T-1
    			for(t = 0; t < T - 1; ++t)
    			{
    				if(step[t].obs == j)
    				{
    					numer += step[t].gamma[i];
    				}
    				denom += step[t].gamma[i];

    			}// next t

    			Bbar[i][j] = numer / denom;
        
    		}// next j
        
    	}// next i
        
	}
}
