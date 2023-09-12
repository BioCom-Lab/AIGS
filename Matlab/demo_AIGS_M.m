clear;
rng(1);
Datasets = {...
    'Yan';
    };
n_datasets = length(Datasets);
options = [];
%The number of marker genes
marker_num = 30;
for id = 1
    data = Datasets{id,1}; disp(data); eval(['load ' data]); options = [];
    [idx_marker_gene,idx_marker] = AIGS_M(fea,marker_num,options);
end