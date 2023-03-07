# An Interpretable Framework for scRNA-seq analysis(scAIG)
#### OVERVIEW
In this repository we provide implementations in both Matlab and Python of scAIG. The main branch of the repository, scAIG, contains the code for the method (available in both Matlab and Python) and sc-seq data called Yan. It is important to note that the provided data are solely intended as examples and should not be substituted for the data presented in the relevant publications. The main content of the code is divided into two parts
  1. demo_scAIG_C: The code for clustering cells into different types. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene. Output a structure called out, including

     num_cluster: Estimation of cell class numbers.
     
     idx_cell:  Index of remaining cells after outlier removal.
     
     idx_cell: Index of remaining gene after gene selection.
     nc: The number of connected branches of the similarity matrix for extimating number of clusterings.
     grp: Cell clustering result.
     grps: Cell clustering results under different estimated number of classes.
     Q: The spectral projection.
     C: The center coodinate according to the grp.
     D_: The distance matrix constructed in gene selection.
     A_: The similarity matrix constructed in gene selection.
     D: The distance matrix constructed for clustering.
     A: The similarity matrix constructed for clustering.
     In python, we adopt the same notation for the same meaning.
  2. demo_scAIG_V: The code for sc-seq data visulization. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene. Output a cell array called Res_Visual_AIG, including four variable:
     The first variable is the original cell data, the second variable is the result of 2D visualization, the third variable is the actual classification result of the cell (removing abnormal cells), and the fourth variable is the average contour index of the projection.
