function D = distfun(x,X)

m = size(X,1);
D = zeros(m,1);
for ii = 1: m
        D(ii) = calMI(x,X(ii,:));
end
