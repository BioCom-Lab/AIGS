function [fea,D,A,idg] = GeneSele(fea,ng,idg,isdropout)
%Determining which genes to select by self-supervision.
%Pseudolabeling from clustering.

n = size(fea,2);  X = zeros(ng,n);  idgs = zeros(ng,1);

D = Dissm(fea);   [A,Q] = GraphEmb(D,5); Q = Q + eps;

for ii = 1:3

    switch isdropout

        case 0
            grp = kmeans(Q(:,2:ii+2),ii+2,'maxiter',1000, ...
                'replicates',50,'Distance','cosine','Options',statset('UseParallel',1));
        case 1
            Q_ = Q(:,2:ii+2)./sqrt(sum(Q(:,2:ii+2).^2,2));
            grp = kmeans(Q_,ii+2,'maxiter',1000, ...
                'replicates',50,'Options',statset('UseParallel',1));

    end
    
    fea_ = (max(grp)-1)*(fea - min(fea,[],1))./(max(fea,[],1)-min(fea,[],1))+1;
    
    fea_ = round(fea_); d = mypdist2(grp',fea_,@distfun); [~,order] = sort(d,'descend');
    
    X((ii-1)*ng + 1:ii*ng,:) = fea(order(1:ng),:); idgs((ii-1)*ng + 1:ii*ng) = idg(order(1:ng));
    
end

fea = unique(X,'row'); idg = unique(idgs,'stable');