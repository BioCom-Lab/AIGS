%% Load data.
load('Colors.mat');
seed = 1;
rng(seed);
Datasets = {...
    'yann';
    'Goolam';
    'Deng';
    'Darmanis';
    'Usoskin';
    'Xin';
    'Muraro';
    'Lake';
    };
n_datasets = length(Datasets); Res_Visual_AIG = cell(n_datasets,4);

for id = 1:8
    data = Datasets{id,1}; disp(data); eval(['load ' data]); options = [];
    [Y,grp,idx_cell] = scAIG_V(fea,options);
    labels = gnd - min(gnd) + 1;
    labels = labels(idx_cell);
    ARI = calARI(grp,labels);
    fprintf('ARI = %2.4f.\n',ARI);
    
    subplot(2,4,id);
    for ii = 1: length(unique(labels))
        hold on;
        scatter(Y(labels == ii,1), Y(labels == ii,2),10,'filled','MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:)); pause(0.1),
    end
    
    s = silhouette(Y,labels); measilhouette = mean(s);
    title([Datasets{id},', Silhouette = ', num2str(round(100*mean(s))/100)]);
    
    xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
    ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);
    
    Res_Visual_AIG{id,1} = data;    
    Res_Visual_AIG{id,2} = Y;
    Res_Visual_AIG{id,3} = labels;
    Res_Visual_AIG{id,4} = measilhouette;
    pause(0.1),
end



