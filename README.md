# About

This repo distributes the code for the Differential Network Enrichment Analysis (DNEA) project. 

# How to Use

1. Input file format : CSV format.   Must have a sample name column followed by a column 
   labelled "GROUP"(all uppercase), followed by metabolite columns e.g.
	
		SAMPLE, GROUP, M1, M2, M3, M4, M5
		S1, 1, 0.514711469, 0.39476782, 0.50289543, 0.651998092, 0.312626386
		S2, 1, 0.616713343, -0.38504144, -0.53754845, -0.303423873, 0.274196675
		S3, 2, -1.119975243, 0.97394155, 0.48339630, -0.170513684, 0.266857172
		S4, 2, 0.246490267, -0.49576216, -0.39343332, 0.255589595, -0.379552788
    
   * Only two values for "GROUP" will be used for the analysis.  More than two produces a warning.
   * Metabolite names should be unique -- duplicate column names are flagged as an error


2. Specify a results directory in the parameters section. For example, if you'd like result files to be written to "/home/wiggie/diffnet_release/results" you'd set 
	
		results_dir = "/home/wiggie/diffnet_release/results"
	
   * All output files will be written to the results directory and script will look for binary files from previous step here.


3. Specify the location of your input file, using the full path e.g. 
	
		input_data = "/home/wiggie/r_data/Main-Standardized.csv"


4. Currently, script is run in three steps by changing flags at the top of `diffnet.R`.  Only one option should be set to `TRUE`; more than one is flagged as an error.

   Analysis steps :   
		a. Tuning parameter selection:  set BIC = TRUE, SS = FALSE, runNetGSA = FALSE
		
		b. Stability selection : set BIC = FALSE, SS = TRUE, runNetGSA = FALSE
		
		c. Netgsa + final output : set BIC = FALSE, SS = FALSE, runNetGSA = TRUE
    

5. Parallel processing : To run 10 iterations each on 10 cluster nodes, set `nCores = 10`  and `nreps = 10`.  Recommended/default is 100 iterations.  


6. Current script screens for the following errors : 

	a.	Missing "GROUP" column
	b.	Results directory not specified or is missing
	c.  Corrects trailing "/" in directory name
	d.	No input file specified, missing input file
	e.	Duplicate columns in input file
	f.	At each step, missing binaries from previous step
	g.	More or less than 2 groups -- warning only
	h.	More than 1 analysis step selected
	i.	No analysis specified
		

