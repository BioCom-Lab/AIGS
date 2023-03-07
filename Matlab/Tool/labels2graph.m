function A = labels2graph(gnd)
n = length(gnd); A = zeros(n,n);
for ii = 1:n
    for jj = 1:ii-1
        if gnd(ii)==gnd(jj)
            A(ii,jj)=1;
        end
    end
end
A = A + A';
