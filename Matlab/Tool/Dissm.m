function  [Do,D,m] = Dissm(fea)


D = pdist2(fea',fea',"spearman");

Do = OBPD(D, 10);

D_min = mink(Do,3,2);

sigma_ = D_min(:,3).^2;

Ao = Do<= min(sigma_,sigma_');

bins = conncomp(graph(Ao));

m = max(bins);