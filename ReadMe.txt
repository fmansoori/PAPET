PAPET v0.1
-------------------------------

These are the instructions to run the PAPET.

You need to install R release 2.0 and PRISM release 4.5 or later.

Before running PAPET do the following:

- Install PRISM from http://www.prismmodelchecker.org/download.php
- Add the path of the PRISM in your system (e.g. C:\Program Files\prism-4.5) to environment variable PRISM_DIR
- Add the "%PRISM_DIR%\lib;%PRISM_DIR%\bin" to the environment variable Path
- Add "%PRISM_DIR%\lib\prism.jar;%PRISM_DIR%\classes;%PRISM_DIR;%PRISM_DIR\lib\*; to the environment variable CP
- Copy the 'myprism' file to the [PRISM directory]\bin ('myprism' file is in the [PAPET directory]\PRISM folder)

- Change the PAPET.program.location in the 'PAPET_v01.R' file to the location of the PAPET in your local machine

---------------------------------
How To use "PAPET":
  
  PAPET ( esetm, group, regulated = "UP", fast_run = TRUE, perms = 100, adjusted = 'NO', threshold = 0.05, max_number_of_genes= 3000 )
  
  - esetm : A matrix containing log transformed and normalized gene expression data. Rows correspond to genes and columns to samples. 
            The first column of the matrix indicates the ENTREZ Id of the corresponding gene.
  - group : A character vector with the class labels of the samples. 
  - regulated : A character value to indicate considering up regulated or down regulated genes.
  - fast_run : A logical value to indicate if using the precomputed random distributions.
  - perms : Number of iterations to determine the gene set score significance p-values
  - adjusted : A logical value to indicate if adjusted p-values of differentially expressed genes for multiple testing
  - threshold : The threshold for computing differentially expressed genes
  - max_number_of_genes : The maximum number of genes considered as differentially expressed. 
  
 See the Example.R as a sample for running the PAPET.


-------------------------------------------------------------------
