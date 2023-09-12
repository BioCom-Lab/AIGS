# Single cell Analyzer with Intelligent Gene selection(AIGS)
#### OVERVIEW
In this repository we provide implementations in both Matlab and Python of AIGS. The main branch of the repository, AIGS, contains the code for the method (available in both Matlab and Python) and sc-seq data called Yan. It is important to note that the provided data are solely intended as examples and should not be substituted for the data presented in the relevant publications. The main content of the code is divided into two parts
  1. demo_AIGS_C: The code for clustering cells into different types. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene, output a structure called out.
  2. demo_AIGS_V: The code for sc-seq data visulization. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene. Output a cell array called Res_Visual_AIG.

#### AIGS
The recent emergence of single-cell RNA sequencing (sc-RNA seq) technology has made it possible to detect cell types at the molecular level. However, clustering cells based on sc-RNA seq data poses two major challenges. The first challenge is related to the high dimensionality of the data, which makes it difficult to use conventional distance measures like Euclidean and cosine distance to accurately capture the spatial distribution of the data. The second challenge arises from the limitations of the sequencing technology, which can result in missing gene expression readings (dropouts) and outliers. While traditional single-cell analysis frameworks like SEURAT and SC3 are widely used, they may not achieve high clustering accuracy. Deep network clustering approaches have shown significant improvements, but they are often considered black-box models. 

To address these challenges, we introduce AIGS, an interpretable framework for accurate and efficient sc-RNA seq analysis. AIGS employs an intelligent gene selection algorithm that targets genes indicative of cell types, which effectively reduces data dimensionality and mitigates the negative impact of dropouts. Additionally, AIGS utilizes a novel scaled-invariant distance metric that emphasizes the similarity of homogeneous cells while increasing the distance between heterogeneous cells. These features make AIGS a promising framework for the accurate and interpretable clustering of sc-RNA seq data.

#### CITATION
The latest version of the manuscript related to AIGS is published on Research Square at https://doi.org/10.21203/rs.3.rs-2738257/v1.

When using AIGS, please cite https://doi.org/10.21203/rs.3.rs-2738257/v1.

#### RUNNING THE MATLAB IMPLEMENTATION
We provide the MATLAB code to run AIGS on eight example data sets in demo_AIGS_C.m, demo_AIGS_V.m and demo_AIGS_M.m for clustering cells, embedding the cells into two dimension and selecting the marker genes through our clustering result. We give the mat format of the Yan data set, and the sources of the rest of the data are as follows:
  1. Yan: Yan, L., Yang, M., Guo, H. et al. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nat Struct Mol Biol 20, 1131–1139 (2013). https://doi.org/10.1038/nsmb.2660.
  2. Goolam: Goolam M, Scialdone A, Graham S J L, et al. Heterogeneity in Oct4 and Sox2 targets biases cell fate in 4-cell mouse embryos. Cell, 2016, 165(1): 61-74. https://doi.org/10.1016/j.cell.2016.01.047.
  3. Deng: Deng, Qiaolin, et al. "Single-cell RNA-seq reveals dynamic, random monoallelic gene expression in mammalian cells. Science 343.6167 (2014): 193-196. http://doi.org/10.1126/science.1245316.
  4. Darmanis: Darmanis, Spyros, et al. "A survey of human brain transcriptome diversity at the single cell level." Proceedings of the National Academy of Sciences 112.23 (2015): 7285-7290. https://doi.org/10.1073/pnas.1507125112.
  5. Usoskin: Usoskin, D., Furlan, A., Islam, S. et al. Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. Nat Neurosci 18, 145–153 (2015). https://doi.org/10.1038/nn.3881.
  6. Xin: Xin, Yurong, et al. "RNA sequencing of single human islet cells reveals type 2 diabetes genes." Cell metabolism 24.4 (2016): 608-615.https://doi.org/10.1016/j.cmet.2016.08.018.
  7. Muraro: Muraro, Mauro J., et al. "A single-cell transcriptome atlas of the human pancreas." Cell systems 3.4 (2016): 385-394. https://doi.org/10.1016/j.cels.2016.09.002.
  8. Lake: Lake, Blue B., et al. "Neuronal subtypes and diversity revealed by single-nucleus RNA sequencing of the human brain." Science 352.6293 (2016): 1586-1590. http://doi.org/10.1126/science.aaf1204.

#### RUNNING THE PYTHON IMPLEMENTATION
We also provide the PYTHON code to run AIGS on the Yan data set in the demo_AIGS_C.py for clustering and demo_AIGS_V.py for visulization. 

Please refer to the directory MATLAB and the file README.txt within for further detail.

#### COMPARISON OF THESE TWO IMPLEMENTATION
Since the program involves a large number of matrix operations, we recommend using Matlab version for faster clustering and visualization.

#### DEBUG
If you encounter any issues while running our code, please do not hesitate to contact us at 12035009@zju.edu.cn. We will be more than happy to assist you.
