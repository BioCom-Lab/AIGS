function [fea,idx_cell] = OutRem(fea,num_knbr,outrate)

D = pdist2(fea',fea',"spearman"); S = exp(-D.^2); n = size(fea,2);

Snk = maxk(S,num_knbr+1); rhos = sum(Snk,1);

outrate = ceil(n*outrate/100);

[~, id_rem] = mink(rhos,outrate); idx_cell = setdiff(1:n,id_rem);

fea(:,id_rem) = [];
