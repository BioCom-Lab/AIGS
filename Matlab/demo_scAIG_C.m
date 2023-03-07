%% Load data.
rng(1);
Datasets = {...
'yann';
'goolamn';
'dengn'
'darmanisn'; 
'usoskinn';
'xinn';
'muraron';
'laken';
};
n_datasets = length(Datasets);
ACC_a = zeros(8,1);
NMI_a = zeros(8,1);
ARI_a = zeros(8,1);
res = zeros(2,8);
for id = 1:8
    data = Datasets{id,1}; disp(data); eval(['load ' data]); options = [];
    data_name = data(1:end-1);
    file_name = strcat(data_name,'.scAIG.mat');
    for rep = 1:10
    
        tic;
        out = scAIG_C(fea,options);
        t_run = toc;

        %% compute ARI.
        labels = gnd - min(gnd) + 1; 
        labels = labels(out.idx_cell);
        ARI = calARI(out.grp,labels);
        
        fprintf('CPU Time = %2.4f.\n',t_run);
        ACC = calAC(out.grp,labels);
        NMI = calMI_1(out.grp,labels);
        fprintf('ACC = %2.4f,NMI=%2.4f,ARI = %2.4f.\n',ACC,NMI,ARI);
        ACC_a(id) = ACC;
        NMI_a(id) = NMI;
        ARI_a(id) = ARI;
    end
end



