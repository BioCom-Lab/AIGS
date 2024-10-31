# Single cell Analyzer with Intelligent Gene selection(AIGS)
#### OVERVIEW
In this repository we provide implementations in both Matlab and Python of AIGS. The main branch of the repository, AIGS, contains the code for the method (available in both Matlab and Python) and scRNA-seq data called Yan. It is important to note that the provided data are solely intended as examples and should not be substituted for the data presented in the relevant publications. The main content of the code is divided into two parts
  1. demo_AIGS_C: The code for clustering cells into different types. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene, output a structure called out.
  2. demo_AIGS_V: The code for scRNA-seq data visulization. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene. Output a cell array called Res_Visual_AIG.
  3. demo_AIGS_M: The code for scRNA-seq data visulization. For Matlab, the input is a single-cell RNA-seq data matrix, where each column is a cell, and each row is the log10 expression of a gene. The other input is the number of marker genes in each cell clusters. Output two gene index array called idx_marker_gene and idx_marker.
#### AIGS
The recent rise of single-cell RNA sequencing (scRNA-seq) has enabled the detection of cell types at the molecular level.   Cell annotation, also known as cell clustering, is a crucial step in scRNA-seq analysis. It enables the identification of cell populations 
\cite{segerstolpe2016single}, determination of cell population topology, and characterization of cellular heterogeneity in complex diseases . The main challenges of clustering come from two aspects. The first challenge arises from the extremely high dimensionality of the RNA-seq data, which presents obstacles in accurately characterizing the spatial distribution of data using conventional distances such as Euclidean and cosine distance.  The second challenge arises from the technical limitations of current sequencing methods, leading to missing gene expression readings (dropouts) and outliers. While classical single-cell analysis frameworks like SEURAT and SIMLR are widely used, they may fall short of achieving high clustering accuracy.  Deep network clustering approaches like scDHA have shown significant improvements in performance due to their extreme representational power. Nonetheless, deep networks are commonly perceived as black-box models, and the mechanism of extracting clustering-relevant information from high-dimensional data remains unclear, which hinders the further development of deep networks.

To overcome these challenges, we present AIGS (single cell Analyzer with Intelligent Gene selection), an interpretable framework designed for accurate and efficient scRNA-seq analysis. The AIGS pipeline comprises modules for gene selection, dimensionality reduction, clustering, visualization, and marker gene identification. AIGS distinguishes itself from other frameworks by utilizing an intelligent gene selection algorithm that targets genes which indicate cell types, a minority of all genes that provide the most informative data on cell types. This gene selector systematically identifies class-indicating genes based on the normalized mutual information (NMI) between the learned pseudo-labels and quantified genes, effectively reducing data dimensionality and mitigating the negative impact of dropouts. Furthermore, AIGS incorporates a novel scale-invariant distance metric that highlights the similarity among homogeneous cells while differentiating heterogeneous cells. Notably, AIGS emphasizes functional clustering, wherein cells with similar functions or states are grouped together, yielding deeper insights into the underlying biological processes. This combination of features establishes AIGS as a promising framework for the analysis of scRNA-seq data. And furthermore, AIGS's outstanding performance does not rely on unpredictable outcomes from deep learning frameworks, making it applicable to other analytical tasks, including spatial transcriptomics and multimodal single-cell genomics

#### CITATION
The latest version of the manuscript related to AIGS is published on Biorxiv at [https://doi.org/10.21203/rs.3.rs-2738257/v1](https://doi.org/10.1101/2024.09.01.610665).


If you use this code, please cite our paper:

@article {Ni2024.09.01.610665,
	author = {Ni, Tianhao and Zhang, Xinyu and Jin, Kaixiu and Pei, Guanxiong and Xue, Nan and Yan, Guanao and Li, Taihao and Li, Bingjie},
	title = {Interpretable scRNA-seq Analysis with Intelligent Gene Selection},
	elocation-id = {2024.09.01.610665},
	year = {2024},
	doi = {10.1101/2024.09.01.610665},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/09/03/2024.09.01.610665},
	journal = {bioRxiv}
}



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
If you encounter any issues while running our code, please do not hesitate to contact us at thni@zju.edu.cn. We will be more than happy to assist you.
