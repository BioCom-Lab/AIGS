We here provide the Python implementation of AIGS. 

**DEMOS**
We provide two demos for the usage of AIGS in both clustering and visulization implementations.

In demo_AIGS_C.py, we provide the script for clustering on eight data sets with the output as running time, ARI. The other output are saved in the sturct called 'out' as follows.

1. num_cluster: Number of cell classes estimated by AIGS.

2. idx_cell: After removing abnormal cells by AIGS, the index of remaining cell.

3. grp: The final clustering result.

4 Q: The spectrol projection for final clustering.

5 C: The category center obtained when Q is clustered with K-means.

6 A: The similarity matrix generated for clustering.

7 idx_gene: After AIGS gene selection, the index of remaining genes.

In demo_AIGS_V.m, we provide the script for clustering on eight data sets with the output as the figure of 2D embedding, ARI, running time and the mean silhouette coefficient between 2D embedding and golden label of each data set. 

#### ATTENTION
Due to the use of stochastic gradient descent in the visualization, the results of python may be different from those of Matlab.
