# An Interpretable Framework for scRNA-seq analysis(scAIG)
#### OVERVIEW
In this repository we provide implementations in both Matlab and Python of scAIG. The main branch of the repository, scAIG, contains the code for the method (available in both Matlab and Python) and sc-seq data called Yan. It is important to note that the provided data are solely intended as examples and should not be substituted for the data presented in the relevant publications. The main content of the code is divided into two parts
  1. demo_scAIG_C: The code for clustering cells into different types. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene, output a structure called out。
  2. demo_scAIG_V: The code for sc-seq data visulization. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene. Output a cell array called Res_Visual_AIG.
   
