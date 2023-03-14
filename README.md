# An Interpretable Framework for scRNA-seq analysis(scAIG)
#### OVERVIEW
In this repository we provide implementations in both Matlab and Python of scAIG. The main branch of the repository, scAIG, contains the code for the method (available in both Matlab and Python) and sc-seq data called Yan. It is important to note that the provided data are solely intended as examples and should not be substituted for the data presented in the relevant publications. The main content of the code is divided into two parts
  1. demo_scAIG_C: The code for clustering cells into different types. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene, output a structure called out.
  2. demo_scAIG_V: The code for sc-seq data visulization. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene. Output a cell array called Res_Visual_AIG.
#### CITATION
The latest version of the manuscript related to scAIG is published on Arxiv at...

When using ScAIG, please cite 

#### RUNNING THE MATLAB IMPLEMENTATION
We provide the MATLAB code to run scAIG on eight example data sets in demo_scAIG_C.m and demo_scAIG_V.m. We give the mat format of the Yan data set, and the sources of the rest of the data are as follows:
  1. Yan: Yan, L., Yang, M., Guo, H. et al. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nat Struct Mol Biol 20, 1131–1139 (2013). https://doi.org/10.1038/nsmb.2660.
  2. Goolam: Goolam M, Scialdone A, Graham S J L, et al. Heterogeneity in Oct4 and Sox2 targets biases cell fate in 4-cell mouse embryos. Cell, 2016, 165(1): 61-74. https://doi.org/10.1016/j.cell.2016.01.047.
  3. Deng: Deng, Qiaolin, et al. "Single-cell RNA-seq reveals dynamic, random monoallelic gene expression in mammalian cells. Science 343.6167 (2014): 193-196. http://doi.org/10.1126/science.1245316.
  4. Darmanis: Darmanis, Spyros, et al. "A survey of human brain transcriptome diversity at the single cell level." Proceedings of the National Academy of Sciences 112.23 (2015): 7285-7290. https://doi.org/10.1073/pnas.1507125112.
  5. Usoskin: Usoskin, D., Furlan, A., Islam, S. et al. Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. Nat Neurosci 18, 145–153 (2015). https://doi.org/10.1038/nn.3881.
  6. Xin: Xin, Yurong, et al. "RNA sequencing of single human islet cells reveals type 2 diabetes genes." Cell metabolism 24.4 (2016): 608-615.https://doi.org/10.1016/j.cmet.2016.08.018.
  7. Muraro: Muraro, Mauro J., et al. "A single-cell transcriptome atlas of the human pancreas." Cell systems 3.4 (2016): 385-394. https://doi.org/10.1016/j.cels.2016.09.002.
  8. Lake: Lake, Blue B., et al. "Neuronal subtypes and diversity revealed by single-nucleus RNA sequencing of the human brain." Science 352.6293 (2016): 1586-1590. http://doi.org/10.1126/science.aaf1204.
