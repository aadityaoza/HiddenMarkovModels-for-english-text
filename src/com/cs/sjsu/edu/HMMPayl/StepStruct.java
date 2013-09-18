package com.cs.sjsu.edu.HMMPayl;

public class StepStruct 
{
	
	public int obs;
    public double c;
    public double alpha[];
    public double beta[];
    public double gamma[];
    public double diGamma[][];
    
    StepStruct()
    {
    	int N = 2 ;
    	alpha = new double[N];
    	beta  = new double[N];
    	gamma = new double[N];
    	diGamma= new double[N][N];
    }
    
    

}
