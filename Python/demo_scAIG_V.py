import numpy as np
import time
from scAIG_C import scAIG_C
import matplotlib.pyplot as plt
import warnings
import random
from scAIG_V import scAIG_V
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
    [Y,grp,idx_cell] = scAIG_V(fea)
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
    print('ARI=%2,4f.\n',ARI)
    print('Embedding CPU Time = %2,4f.\n',t2-t1)
    print('Silhouette Coefficient = %2,4f.\n',si)