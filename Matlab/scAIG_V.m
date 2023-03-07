function [Y,grp,idx_cell] = scAIG_V(Y,options)

n_epochs = 500;
min_dist = 0.1;

out = scAIG_C(Y,options); grp = out.grp;

idx_cell = out.idx_cell; C = out.C;

Q = out.Q; Q_ = Q./sqrt(sum(Q.^2,2));

Y = Y(out.idx_gene,:);
Y = Y(:,out.idx_cell);

Fea = [Y',mean(Y(Y>0))*Q_]; 

D = Dissm(Fea');

D_min = mink(D,7,2); sigma = D_min(:,7).^2;
A = exp(-bsxfun(@rdivide,D.^2,max(sigma,sigma')+eps));
A = bsxfun(@times,A, D <= 2*max(sigma,sigma'));
A = sparse(A); 

%% Initial 2D Points
Y = vis_rand(C,grp,Q);
Y = (Y-min(Y,[],'all'))./(max(Y,[],'all')-min(Y,[],'all'));
Y = 2*Y;

%fit
D_Y = pdist2(Y,Y);
spread = max(D_Y,[],'all');
phi = fittype('1/(1+a*x^(2*b))','independent','x','coefficients',{'a','b'});
x_fit = linspace(0,spread*3,300);
y_fit = exp((-x_fit+min_dist));
y_fit(x_fit<=min_dist) = 1;
cfun = fit(x_fit',y_fit',phi,'StartPoint',[1 1],'Algorithm', 'Levenberg-Marquardt');
a = cfun.a; b = cfun.b;

%embedding
Y = optimize_layout_euclidean(Y,n_epochs,a,b,A);
Y = (Y-min(Y))./(max(Y)-min(Y));




