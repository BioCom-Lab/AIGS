clear;clc;
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
for id = 6%1:n_datasets
    load(Datasets{id});
    options = [];
    out = scAIG(fea,options,gnd);
    idx_cell = out.idx_cell;
    grp = out.grp;
    Q = out.Q;
    C = out.C;
    idx_genesele = out.idx_genesele;
    NMI = out.NMI;
    labels = gnd(idx_cell);
    [labels,index] = sort(labels);
    fea = fea(idx_genesele,:);
    fea = fea(:,index);
    [m,n] = size(fea);
    n_class = max(gnd);
    %marker
    res_marker_class = cell(n_class,1);
    mean_fea = zeros(m,n_class);
    index_marker = [];
    for ii = 1:m
        temp = fea(ii,:);
        for jj = 1:n_class
            mean_fea(ii,jj) = mean(temp(labels==jj));
        end

    end
    for ii = 1:m
        gap = zeros(n_class,n_class);
        for kk = 1:n_class
            for ll = 1:n_class
                gap(kk,ll) = mean_fea(ii,kk)-mean_fea(ii,ll);
            end
        end
        gap = abs(gap); gap = mean(gap);
        max_gene = max(mean_fea(ii,:));
        for jj = 1:n_class
            if gap(jj)>max_gene*0.7
                res_marker_class{jj} = [res_marker_class{jj},ii];
                index_marker = [index_marker,ii];
                break;
            end
        end

    end
    NMI_marker = NMI(index_marker);

    for ii = 1:length(res_marker_class)
        res_marker_class{ii} = res_marker_class{ii}';
    end
        num_marker = zeros(n_class,1);
        fea = log10(fea+eps);
        for ii = 1:length(index_marker)
            disp(index_marker(ii))
            subplot(2,1,1);
            imagesc(fea(index_marker(ii),:));axis off;
            subplot(2,1,2);
            imagesc(labels'); axis off;
        end
end
% fea = fea(:,idx_cell);
% fea = log10(fea+1e-5);
% for ii = 1:m
%     clf;
%     X = labels;
%     Y = fea(ii,:);
%     Hdl=violinChart(gca,X,Y,[0 0.447 0.741]);
% end