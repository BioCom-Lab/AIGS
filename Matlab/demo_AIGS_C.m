%% Load data.
Datasets = {...
    'Yan';
    };
n_datasets = length(Datasets);
for id = 1
    data = Datasets{id,1}; disp(data); eval(['load ' data]); 
    options = [];
        tic;
        out = AIGS_C(fea,options);
        t_run = toc;
        %% compute ARI.
        labels = gnd - min(gnd) + 1; 
        labels = labels(out.idx_cell);     
        ARI = calARI(out.grp,labels);
        fprintf('ARI = %2.4f.\n',ARI);
        fprintf('CPU Time = %2.4f.\n',t_run);
end



