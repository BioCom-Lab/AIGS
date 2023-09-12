function [idx_marker_gene,idx_marker] = AIGS_M(fea,marker_num,options)
%% load options
% Default options
filter_thr = [log2(3) 3/2];
num_gene_min = 100;

% Predicted options
if isfield(options,'num_gene_min'), num_gene_min = options.num_gene_min;     end
if isfield(options,'filter_thr'), filter_thr = options.filter_thr;     end


%% step 1: initial gene filtering

idg = 1:size(fea,1); thr_1 = filter_thr(1); thr_2 = filter_thr(2);

absmax = max(fea,[],2); fea = fea(absmax>thr_1,:); idg = idg(absmax>thr_1);

absvar = var(fea,[],2); [pod, i_pod] = maxk(absvar,2*num_gene_min);

if pod(end)>thr_2
    fea = fea(absvar>thr_2,:); idg = idg(absvar>thr_2);
else
    fea = fea(i_pod,:); idg = idg(i_pod);
end
%% step 2: cluster and finding marker genes.
out = AIGS_C(fea,[]);
grp = out.grp;
n_class = max(grp);
idx_marker = cell(n_class,1);
fea = fea(:,out.idx_cell);
n = size(fea,2);
for i_class = 1:n_class
    grp_temp = zeros(n,1);
    grp_temp(grp==i_class) = 1;
    Marker = findMarkerGenesp(fea',grp_temp,marker_num);
    idx_marker{i_class} = Marker;
end
idx_marker_gene = zeros(n_class*marker_num,1);
for ii = 1:n_class
    idx_marker{ii} = idg(idx_marker{ii});
    idx_marker_gene((ii-1)*marker_num+1:ii*marker_num) = idx_marker{ii};
end

end