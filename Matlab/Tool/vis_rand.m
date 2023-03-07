function Qe = vis_rand(C,grp,Q)

n = size(C,1); m = length(grp); Qe = zeros(m,2);

for kk = 1:n
    C(kk,:) = mean(Q(grp==kk,:),1);
end

theta = 0:2*pi/n:2*pi-2*pi/n;
T = [cos(theta'),sin(theta')]; 
dist = 0.01*norm(T(1,:)-T(2,:));
for ii = 1:n
    Dist = Q(grp == ii,:) - C(ii,:);  
    weight = (sum(Dist.^2,2)).^0.1; weight = dist*weight/max(weight);
    Qe(grp == ii,:) = T(ii,:) + weight.*randn(sum(grp == ii),2);
end
Qe = 2*Qe;