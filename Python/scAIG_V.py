from scipy.optimize import curve_fit
import scipy
import networkx as nx     
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy import sparse
from sklearn.utils import check_random_state
from scAIG_C import scAIG_C
def mink(A,k):
    order = np.argsort(A,axis = 0)   
    n = np.size(A,axis = 0)
    res = np.zeros([k,n])
    for ii in range(n):
        res[:,ii] = A[order[range(k),ii],ii]
    return [res,order]
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
def find_ab_params(spread, min_dist):
    def curve(x, a, b):
        return 1.0 / (1.0 + a * x ** (2 * b))
    xv = np.linspace(0, spread * 3, 300)
    yv = np.zeros(xv.shape)
    yv[xv < min_dist] = 1.0
    yv[xv >= min_dist] = np.exp(-(xv[xv >= min_dist] - min_dist))
    params, covar = curve_fit(curve, xv, yv)
    return params[0], params[1]
def vis_rand(C,grp,Q):
    n = C.shape[0]
    m = len(grp)
    Qe = np.zeros([m,2])
    for kk in range(n):
        C[kk,:] = np.mean(Q[grp==kk,:],axis = 0)
    theta = np.linspace(0,2*np.pi-2*np.pi/n,n)
    T = np.array([np.cos(theta),np.sin(theta)])
    T = T.T
    dist = 0.01*np.linalg.norm(T[0,:]-T[1,:])
    for ii in range(n):
        Dist = Q[grp==ii,:]-C[ii,:]
        weight = np.sum(Dist*Dist,axis = 1)**(0.1)
        weight = dist*weight/np.max(weight)
        for jj in range(2):
            Qe[grp==ii,jj]=T[ii,jj]+weight*np.random.randn(np.sum(grp==ii))
    return Qe
def clip(val):
    if val>4:
        val = 4
    elif val<-4:
        val = -4
    return val
def optimize_each_point(embedding,A,a,b,alpha,loop1=20,loop2=20):
    n = embedding.shape[0]
    for ii in range(n):
        temp = np.where(A[ii,:]!=0)[0]
        current = embedding[ii,:]
        for loop in range(loop1):
            j = np.random.randint(len(temp))
            jj = temp[j]
            if np.random.rand()<A[ii,jj]:
                other = embedding[jj,:]
                dist_squared = np.sum((current-other)**2)
                if dist_squared!=0:
                    grad_coeff = -2*a*b*dist_squared**(b-1)
                    grad_coeff /= (a*dist_squared**b+1)
                else:
                    grad_coeff = 0
                for d in range(2):
                    grad_d = clip(grad_coeff*(current[d]-other[d]))
                    current[d] = current[d]+alpha*grad_d
        for loop in range(loop2):
            j = np.random.randint(n)
            jj = j
            if jj!=j:
                break
            if np.random.rand()<1-A[ii,jj]:
                other = embedding[jj,:]
                dist_squared = np.sum((current-other)**2)
                if dist_squared!=0:
                    grad_coeff = 2*b
                    grad_coeff /= ((0.001+dist_squared)*(a*dist_squared**b+1))
                elif jj == ii:
                    continue
                else:
                    grad_coeff = 0
                for d in range(2):
                    if grad_coeff>0:
                        grad_d = clip(grad_coeff*(current[d]-other[d]))
                    else:
                        grad_d = 4
                    current[d] = current[d]+alpha*grad_d
        embedding[ii,:] = current
    return embedding         
def optimize_embedding(embedding,n_epochs,a,b,A,initial_alpha = 1,loop1 = 20,loop2 = 20):
    alpha = initial_alpha
    for n in range(n_epochs):
        if n%10==0:
            print('Iteration=',n)
        if n==n_epochs-1:
            print('Embedding Finish')
        embedding = optimize_each_point(embedding,A,a,b,alpha,loop1,loop2)
        alpha = initial_alpha*(1-n/n_epochs)
    return embedding

def scAIG_V(Y):
    n_epochs = 500
    min_dist = 0.1
    [num_cluster,idx_cell,grp,Q,C,A,idx_gene] = scAIG_C(Y)
    Q_ = Q/np.linalg.norm(Q,ord = 2,axis = 1,keepdims = True)
    Y = Y[idx_gene,:]
    Y = Y[:,idx_cell]
    Fea = np.concatenate((Y.T,np.mean(Y[Y>0])*Q_),axis=1)
    D = Dissm(Fea.T)[0]
    sigma = mink(D,7)[0]
    sigma = sigma[-1,:]
    sigma = np.reshape(sigma,[len(sigma),1])
    sigma = sigma*sigma
    temp = np.maximum(sigma,sigma.T)+np.spacing(1)
    A = np.exp(-D*D/temp)
    A[np.where(D>2*temp)]=0
    n = A.shape[0]
    Y = vis_rand(C,grp,Q)
    Y = (Y-np.min(Y))/(np.max(Y)-np.min(Y))
    Y = 2*Y  
    D_Y = pdist(Y)
    spread = np.max(D_Y)
    
    a,b = find_ab_params(spread, min_dist)
    Y = optimize_embedding(Y,n_epochs,a,b,A)
    Y = (Y-np.min(Y))/(np.max(Y)-np.min(Y))
    return [Y,grp,idx_cell]