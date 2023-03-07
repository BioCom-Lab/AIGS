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
    for rep = 1%1:10
        
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
    %% Two-dimensional visualization.
    A = out.A;
    Q = out.Q;
    fea = out.fea;
    idx_cell = out.idx_cell;
    D = pdist2(fea',fea',"spearman");
    Do = OBPD(D,size(A,1));
    D = D-diag(diag(D));
    Do = Do-diag(diag(Do));
    SI_old = silhouette(fea',gnd - min(gnd) + 1,squareform(D));
    SI_new = silhouette(fea',gnd - min(gnd) + 1,squareform(Do));
    boxplot([SI_old,SI_new])
    res(1,id) = mean(SI_old);
    res(2,id) = mean(SI_new);
    fprintf('=%3.2f\n',mean(SI_old));
    fprintf('=%3.2f\n',mean(SI_new));
    C = cell(size(fea,2)+1,2);
    C{1,1} = 'spearman';
    C{1,2} = 'order';
    for ii = 1:size(fea,2)
        C{ii+1,1} = SI_old(ii);
        C{ii+1,2} = SI_new(ii);
    end
    file = strcat(data_name,'.xls');
    writecell(C,file)
end



