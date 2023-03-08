function out = scAIG_C(fea,options)
%% load options
% Default options
filter_thr = [log2(3) 3/2];
num_gene_min = 100;
num_knbr = 10;
outrate  = 5;
isdropout = 0;

% Predicted options
if isfield(options,'num_gene_min'), num_gene_min = options.num_gene_min;     end
if isfield(options,'filter_thr'), filter_thr = options.filter_thr;     end
if isfield(options,'num_knbr'), num_knbr = options.num_knbr;     end
if isfield(options,'outrate'), outrate = options.outrate;     end


%% step 1: initial gene filtering 

idg = 1:size(fea ,1); thr_1 = filter_thr(1); thr_2 = filter_thr(2);

absmax = max(fea,[],2); fea = fea(absmax>thr_1,:); idg = idg(absmax>thr_1);

absvar = var(fea,[],2); [pod, i_pod] = maxk(absvar,2*num_gene_min);

if pod(end)>thr_2
   fea = fea(absvar>thr_2,:); idg = idg(absvar>thr_2);
else    
   fea = fea(i_pod,:); idg = idg(i_pod);
end

%% step 2: gene selection.
[fea,D_,A_,idg] = GeneSele(fea,num_gene_min,idg);


%% step 3: Outlier Removing.
[fea,idc] = OutRem(fea,num_knbr,outrate);
out.fea = fea; 
%% step 4：Clustering.
[D,~,nc] = Dissm(fea); [A,Q] = GraphEmb(D,nc + 3); Q = Q + eps;
ratio_max = 0; 

grps = zeros(size(fea,2),4);
for rep = 1:4
    num_cluster = nc + rep - 1;
    Q_ = Q(:,2:num_cluster)./sqrt(sum(Q(:,2:num_cluster).^2,2));
    [grp,C] = kmeans(Q_,num_cluster,'maxiter',1000,'replicates',10,'Options',statset('UseParallel',1));
    grps(:,rep) = grp;
    B = labels2graph(grp);
    As = sign(A);
    ratio = sqrt(num_cluster)*(sum(sum(As.*B))/sum(B(:)))/(sum(sum(As.*~B))/sum(~B(:)));
    if ratio > ratio_max
        grp_final = grp;
        num_cluster_final = num_cluster;
        Q_final = Q(:,2:num_cluster);
        ratio_max = ratio;
        C_final = C;
    end    
end

%% 
out.num_cluster = num_cluster_final;
out.idx_cell = idc; 
out.idx_gene = idg; 
out.nc = nc;
out.grp = grp_final;
out.grps = grps; 
out.Q = Q_final;
out.C = C_final;
out.D_ = D_;
out.A_ = A_;
out.D = D;
out.A = A;
