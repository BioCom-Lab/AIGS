function DO = OBPD(D,knn)
n = size(D,1); DO = n*ones(n,n);
nknns = 1:knn;
for ii = 1:n
    [~,a] = mink(D(ii,:),knn);
    DO(ii,a) = nknns-1;
end
DO = min(DO,DO');
