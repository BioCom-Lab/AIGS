import numpy as np
import networkx as nx               
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
import scipy.stats
import scipy
def mink(A,k):
    order = np.argsort(A,axis = 0)   
    n = np.size(A,axis = 0)
    res = np.zeros([k,n])
    for ii in range(n):
        res[:,ii] = A[order[range(k),ii],ii]
    return [res,order]
def maxk(A,k):
    order = np.argsort(A,axis = 0)
    n = np.size(A,axis = 0)
    res = np.zeros([k,n])
    for ii in range(n):
        res[:,ii] = A[order[-k-1:-1,ii],ii]
    return [res,order]
def labels2graph(grp):
    n = len(grp)
    A = np.zeros([n,n])
    for ii in range(n):
        for jj in range(ii):
            if grp[ii]==grp[jj]:
                A[ii,jj] = 1
    A = A+A.T
    return A
def calMI(idx,label):
    idx = idx-np.min(idx)
    label = label-np.min(idx)
    idx = idx.reshape(len(idx),1)
    label = label.reshape(len(label),1)
    label = label.astype(int)
    n = max(np.max(idx),np.max(label))+1
    Perf = np.zeros([n,n])
    P1 = np.zeros([n,1])
    P2 = np.zeros([n,1])
    for ii in range(len(idx)):
        Perf[idx[ii],label[ii]] = Perf[idx[ii],label[ii]]+1
        P1[idx[ii]] = P1[idx[ii]]+1
        P2[label[ii]] = P2[label[ii]]+1
    Perf = Perf/len(idx)
    P1 = P1/len(idx)
    P2 = P2/len(idx)
    mi = np.sum(Perf*np.log2((Perf+np.spacing(1))/(np.dot(P1,P2.T)+np.spacing(1))))
    H1 = -np.sum(P1*np.log2(P1+np.spacing(1)))
    H2 = -np.sum(P2*np.log2(P2+np.spacing(1)))
    mi = mi/(max(H1,H2)+np.spacing(1))
    return mi
def OBPD(D,knn):
    n = np.size(D,0)
    Do = n*np.ones_like(D)
    n_knns = np.linspace(1,knn,knn)
    for ii in range(n):
        temp = np.argsort(D[ii,:],axis = 0)
        a = temp[0:knn]
        Do[ii,a] = n_knns-1
    Do = np.minimum(Do,Do.T)
    return Do

def conncomp(A): 
    graph = nx.from_numpy_matrix(A)
    return nx.number_connected_components(graph)


def Dissm(fea):
    D = 1-pdist(fea.T,lambda u,v:scipy.stats.spearmanr(u,v)[0])
    D = squareform(D)
    Do = OBPD(D,10)
    D_min = mink(Do,3)[0]
    D_min = D_min[-1,:]
    D_min = D_min.reshape(len(D_min),1)
    sigma_ = D_min*D_min
    Ao = (Do<=np.minimum(sigma_,sigma_.T))
    m = conncomp(Ao)
    return [Do,D,m]

def GraphEmb(Do,n_class):
    sigma = mink(Do,7)[0]
    sigma = sigma[-1,:]
    sigma = np.reshape(sigma,[len(sigma),1])
    sigma = sigma*sigma
    temp = np.maximum(sigma,sigma.T)+np.spacing(1)
    A = np.exp(-Do*Do/temp)
    A[np.where(Do>2*temp)]=0
    d = np.sqrt(np.sum(A,axis = 0))
    d = np.diag(d**(-1)) 
    L = d@A@d
    L = (L+L.T)/2
    n = np.size(L,1)
    L[np.where(L<=np.spacing(1))]=0
    L = scipy.sparse.csr_matrix(L)
    [value,vector] = eigsh(L,n_class,ncv = min(1000,n),which = 'LM')
    temp = np.argsort(value.real)[::-1]
    Q = vector[:,temp]
    Q = Q.real
    return [A,Q]
def GeneSele(fea,ng,idg):
    ##Determining which genes to select by self-supervision
    ##Pseudolabeling from clustering
    n = np.size(fea,1)
    X = np.zeros((3*ng,n))
    D = Dissm(fea)[0]
    [_,Q] = GraphEmb(D,5)
    Q = Q+np.spacing(1)
    idx_gene = np.zeros((3*ng))
    for ii in range(1,4):
        Q1 = Q[:,1:ii+2]
        Q1 = Q1/np.linalg.norm(Q,ord = 2,axis = 1,keepdims = True)
        grp = KMeans(n_clusters=ii+2,max_iter = 1000,init = 'k-means++', n_init = 50).fit(Q1).labels_
        fea_ = (np.max(grp))*(fea-np.min(fea,axis = 0))/(np.max(fea,axis = 0)-np.min(fea,axis = 0))+1
        fea_ = np.round(fea_)
        d = np.zeros(fea_.shape[0])
        for jj in range(fea_.shape[0]):
            d[jj] = calMI(grp,fea_[jj,:])
        order = np.argsort(-d)
        X[range((ii-1)*ng,ii*ng),:] = fea[order[0:ng],:]
        idx_gene[range((ii-1)*ng,ii*ng)] = order[0:ng].reshape(-1)
    fea = np.unique(X,axis = 0)
    idx_gene = np.unique(idx_gene)
    idg = idg.astype(int)
    idx_gene = idg[idx_gene.astype(int)]
    return [fea,D,idx_gene]
def OutRem(fea,num_knbr,outrate):
    D = 1-pdist(fea.T,lambda u,v:scipy.stats.spearmanr(u,v)[0])
    D = squareform(D)
    S = np.exp(-D*D)
    n = np.size(D,0)
    Snk = maxk(S,num_knbr+1)[0]
    rhos = np.sum(Snk,axis=0)
    outrate = int(np.ceil(n*outrate/100))
    id_rem = np.argsort(rhos)
    id_rem = id_rem[0:outrate]
    idx_cell = list(set(range(n))-set(id_rem))
    fea = fea[:,idx_cell]
    return [fea,idx_cell]
def scAIG_C(fea,
    filter_thr=[np.log2(3),3/2],
    num_gene_min = 100,
    num_knbr = 10,
    outrate = 5):
    ## step 1: initial gene filtering
    idg = np.linspace(0,fea.shape[0]-1,fea.shape[0])
    thr_1 = filter_thr[0]
    thr_2 = filter_thr[1]
    absmax = np.max(fea,axis = 1)
    fea = fea[absmax>thr_1,:]
    idg = idg[absmax>thr_1]
    absvar = np.var(fea,axis = 1,ddof = 1)
    temp = np.argsort(absvar,axis = 0)
    i_pod = temp[-2*num_gene_min:][::-1]
    pod = absvar[i_pod]
    if pod[-1]>thr_2:
        fea = fea[absvar>thr_2,:]
        idg = idg[absvar>thr_2]
    else:
        fea = fea[i_pod,:]
        idg = idg[i_pod]
    ## gene selection 
    [fea,_,idx_gene] = GeneSele(fea,num_gene_min,idg)
    ## Outlier Removing
    [fea,idx_cell] = OutRem(fea,num_knbr,outrate)
    ## Clustering
    [D,_,nc] = Dissm(fea)
    [A,Q] = GraphEmb(D,nc+3)
    A_s = np.sign(A)
    Q = Q+np.spacing(1)
    ratio_max = 0
    grps = np.zeros((fea.shape[1],4))
    for rep in range(1,5):
        num_cluster = nc+rep-1
        Q1 = Q[:,1:num_cluster]
        Q1= Q1/np.linalg.norm(Q1,ord = 2,axis = 1,keepdims = True)
        KMeans_res = KMeans(n_clusters=num_cluster,max_iter = 1000,init = 'k-means++', n_init = 50).fit(Q1)
        C = KMeans_res.cluster_centers_
        grp = KMeans_res.labels_
        grps[:,rep-1] = grp
        B = labels2graph(grp)
        ratio = np.sqrt(num_cluster)*np.sum(A_s*B)*(np.sum(1-B))/np.sum(A_s*(1-B))/np.sum(B)
        if ratio > ratio_max:
            grp_final = grp
            num_cluster_final = num_cluster
            Q_final = Q[:,range(1,num_cluster)]
            ratio_max = ratio
            C_final = C
    return [num_cluster_final,idx_cell,grp_final,Q_final,C_final,A,idx_gene]
