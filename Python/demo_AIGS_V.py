import numpy as np
import time
from AIGS_C import AIGS_C
import matplotlib.pyplot as plt
import warnings
import random
from AIGS_V import AIGS_V
import warnings
from sklearn.metrics.cluster import adjusted_rand_score as ari
from sklearn import metrics
warnings.filterwarnings("ignore")
# seed = 42
# random.seed(seed)
Datasets = ['Yan']
n_datasets = len(Datasets)
Res_Visual_AIG = []
for id in range(1):
    data = Datasets[id]
    print(data)
    fea = np.load(data+'_cell.npy')
    gnd = np.load(data+'_label.npy')
    t1 = time.time()
    [Y,grp,idx_cell] = AIGS_V(fea)
    t2 = time.time()
    labels = gnd-np.min(gnd)+1
    labels = labels[idx_cell]
    labels = labels-1
    labels = labels.reshape(grp.shape)
    ARI = ari(labels,grp)
    for ii in range(np.max(labels)):
        plt.scatter(Y[labels==ii,0],Y[labels==ii,1])
    plt.show()
    si = metrics.silhouette_score(Y,labels)
    print('ARI=%s.\n' %(ARI))
    print('Embedding CPU Time = %s.\n' %(t2-t1))
    print('Silhouette Coefficient = %s.\n' %(si))