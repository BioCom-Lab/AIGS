function [fea,D,NMI,index_fea] = GeneSele(fea,ng,idg,isdropout)
%Determining which genes to select by self-supervision.
%Pseudolabeling from clustering.

n = size(fea,2);  X = zeros(ng,n); 

D = Dissm(fea);   [~,Q] = GraphEmb(D,5); Q = Q + eps;

NMI = zeros(3*ng,1);

index_fea = zeros(3*ng,1);
for ii = 1:3
        
    grp = kmeans(Q(:,2:ii+2),ii+2,'maxiter',1000,'replicates',50,'Distance','cosine','Options',statset('UseParallel',1));
    
    fea_ = (max(grp)-1)*(fea - min(fea,[],1))./(max(fea,[],1)-min(fea,[],1))+1;
    
    fea_ = round(fea_); d = mypdist2(grp',fea_,@distfun); [~,order] = sort(d,'descend');
    
    X((ii-1)*ng + 1:ii*ng,:) = fea(order(1:ng),:);
    
    index_fea((ii-1)*ng + 1:ii*ng) = order(1:ng);
    
    NMI((ii-1)*ng + 1:ii*ng) = d(order(1:ng));
end

fea = unique(X,'row');

 [index_fea,order] = unique(index_fea);

NMI = NMI(order);