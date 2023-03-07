function [idx,Q,C] = spectral_proj(A,dim)

A = sparse(A);
n = size(A,1);
d = sqrt(sum(A,2));
L = bsxfun(@ldivide,d,A);
L = speye(n)-bsxfun(@rdivide,L,d');
L = (L + L')/2;
L = full(L);
[Q,~] = eigs(L,dim,'smallestabs','SubspaceDimension',min(1000,n));
Q = Q + eps;

%rng default; 

[idx, C] = kmeans(Q,dim,'maxiter',5000,'replicates',10,'Distance','cosine');
end