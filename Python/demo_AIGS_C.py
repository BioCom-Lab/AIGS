import numpy as np
import time
from AIGS_C import AIGS_C
from sklearn.metrics.cluster import adjusted_rand_score as ari
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
Datasets = ['Yan']
n_datasets = len(Datasets)
for id in range(1):
    data = Datasets[id]
    print(data)
    fea = np.load(data+'_cell.npy')
    gnd = np.load(data+'_label.npy')
    for rep in range(1):
        t1 = time.time()
        [num_cluster,idx_cell,grp,Q,C,A,idx_gene] = AIGS_C(fea)
        t2 = time.time()
        labels = gnd-np.min(gnd)+1
        labels = labels[idx_cell]
        labels = labels-1
        labels = labels.reshape(grp.shape)
        ARI = ari(labels,grp)
        print('ARI=%2,4f.\n',ARI)
        print('Clustering CPU Time = %2,4f.\n',t2-t1)
