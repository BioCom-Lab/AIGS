clear;
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
%The number of marker genes
marker_number = 5;
for id = 1%1:8
    data = Datasets{id,1}; disp(data); eval(['load ' data]); options = [];
    out = scAIG_C(fea,options);
    labels = gnd - min(gnd) + 1;
    labels = labels(out.idx_cell);
    grp = out.grp;
    [grp,idx] = sort(grp);
    idx_cell = out.idx_cell;
    idx_gene_gf = out.idg_gf;
    fea = out.fea_gf;
    fea = fea(:,idx_cell);
    fea = fea(:,idx);
    [n_gene,n_cell] = size(fea);
    n_class = max(grp);
    idx_marker_grp = cell(n_class,1);
    idx_marker_gene_name = cell(n_class,1);
    for i_class = 1:n_class
        grp_temp = zeros(n_cell,1);
        grp_temp(grp==i_class) = 1;
        Marker = findMarkerGenes(fea',grp_temp,marker_number);
        idx_marker_grp{i_class} = Marker;
        idx_marker_gene_name{i_class} = idx_gene_gf(Marker);
    end
    %plot the heatmap of marker gene and our clustering result.
    subplot(n_class+1,1,1);
    imagesc(grp');
    for i_class = 1:n_class
        subplot(n_class+1,1,i_class+1);
        heatmap(fea(idx_marker_grp{i_class},:)); colormap('parula')
    end
end
