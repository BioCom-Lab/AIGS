function [A,Q] = GraphEmb(Do,n_class)

warning off;
D_min = mink(Do,7,2); sigma = D_min(:,7).^2;
A = exp(-bsxfun(@rdivide,Do.^2,max(sigma,sigma')+eps));
A = bsxfun(@times,A, Do <= 2*max(sigma,sigma')); 
A = sparse(A); 
n = size(A,1);
d = sqrt(sum(A,2));
L = bsxfun(@ldivide,d,A);
L = bsxfun(@rdivide,L,d');
L = (L + L')/2;

[Q,~] = eigs(L,n_class,'largestabs','SubspaceDimension',min(1000,n));
