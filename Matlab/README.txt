We here provide the Matlab implementation of scAIG. 

**DEMOS**
We provide two demos for the usage of scAIG in both clustering and visulization implementations.

In demo_scAIG_C.m, we provide the script for clustering on eight data sets with the output as running time, ACC, NMI and ARI. The other output are saved in the sturct called 'out' as follows.

1. num_cluster: Number of cell classes estimated by scAIG.

2. idx_cell: After removing abnormal cells by scAIG, the index of remaining cell.

3. idx_gene: After scAIG gene selection, the index of remaining genes.

4. nc: Number of connected branches for estimating number of clustering

5. grp: The final clustering result.

6. grps: The clustering result under different number of clustering.

7. Q: The spectrol projection for final clustering.

8. C: The category center obtained when Q is clustered with K-means.

9. D_: The distance matrix generated in the gene selection.

10. A_: The similarity matrix generated in the gene selection.

11. D: The distance matrix generated for clustering.

12. A: The similarity matrix generated for clustering.
