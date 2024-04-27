clear;clc;close all
warning off; id = 13;

Jobs = {
    1,  'Fig 2 a'
    2,  'Fig 2 b'
    3,  'Fig 2 c'
    4,  'Fig 2 d'
    5,  'Fig 3 a'
    6,  'Fig 3 b'
    7,  'Fig 3 c'
    8,  'Fig 4 a'
    9,  'Fig 4 b'
    10, 'Fig 5 a'
    11, 'Fig 5 b'
    12  'Fig 5 c'
    13  'Fig 5 d'
    };

rng default;
Datasets = {...
    'Yan';
    'Goolam';
    'Deng'
    'Darmanis';
    'Usoskin';
    'Xin';
    'Muraro';
    'Lake';
    };
filter_thr = [log2(3) 3/2];
num_gene_min = 100;
num_knbr = 10;
outrate  = 5;
task = Jobs{id,2};
switch task
    case 'Fig 2 a'
        for fakeloop = 1
            load('Res_AIG_one.mat');

            colors =[0.8,0.9,1.0;
                0.6784,0.8471,0.9020;
                0.2117,0.3921,0.5450];

            for ii = 1:2
                subplot(1,2,ii);
                out = Res{ii};
                gene_number = out.num_gene;
                gene_number(3) = gene_number(3)*10;
                hold on
                for jj = 1:3
                    b = bar(jj,gene_number(jj),0.75,'stacked');
                    set(b(1),'FaceColor',colors(jj,:));
                end
                title(Datasets{ii})
                box off
                ylabel('Number of Genes','fontsize',10)
                set(gca,'XTickLabel','');
                xlabel('Stage','fontsize',10)
                legend({'Raw','Filtered','Selected'},'fontsize',10,'NumColumns',1); hold on
                legend('box', 'off');
            end
        end
    case 'Fig 2 b'
        for fakeloop = 1
            point_size = 5;
            load('Res_AIG_one.mat');
            for id = 1:2
                load(Datasets{id});
                n_class = length(unique(gnd));
                num_cell_class = zeros(n_class,1);

                for ii = 1:n_class
                    num_cell_class(ii) = sum(gnd==ii);
                end
                [xaxis,yaxis] = plot_graph(fea,num_cell_class,gnd);
                subplot(1,2,id);
                f1 = scatter(xaxis,yaxis,point_size*ones(size(xaxis)),[.75,.75,.75],'filled');
                hold on

                out = Res{id};
                idg = out.idx_gene;
                fea_new = fea(idg,:);
                [xaxis,yaxis] = plot_graph(fea_new,num_cell_class,gnd);
                f2 = scatter(xaxis,yaxis,point_size*ones(size(xaxis)),'red','filled');
                set(gca,'xticklabel',[])
                set(gca,'yticklabel',[])
                set(gca,'xtick',[],'xticklabel',[])
                set(gca,'ytick',[],'yticklabel',[])
                xlabel('intra-class variance')
                ylabel('inter-calss variance')
                title(Datasets{id});
                lgd = legend([f2],'Selected genes','FontSize',10,'NumColumns',1);
                legend('box','off')
            end
        end
    case 'Fig 2 c'
        for fakeloop = 1
            res_sil = cell(3,2);
            for data_ii =1:2
                data = Datasets{data_ii,1}; disp(data); eval(['load ' data]); options = [];
                idg = 1:size(fea ,1); thr_1 = filter_thr(1); thr_2 = filter_thr(2);
                absmax = max(fea,[],2); fea = fea(absmax>thr_1,:); idg = idg(absmax>thr_1);
                absvar = var(fea,[],2); [pod, i_pod] = maxk(absvar,2*num_gene_min);
                if pod(end)>thr_2
                    fea = fea(absvar>thr_2,:); idg = idg(absvar>thr_2);
                else
                    fea = fea(i_pod,:); idg = idg(i_pod);
                end
                %                 fea = GeneSele(fea,100,idg,'');
                [Do1,D,m] = Dissm(fea);
                res1 = sil(fea,gnd,D);
                res2 = sil(fea,gnd,Do1);
                res3 =sil(fea,gnd,pdist2(fea',fea','euclidean'));
                res_sil{1,data_ii} = res1;
                res_sil{2,data_ii} = res2;
                res_sil{3,data_ii} = res3;
            end
            res1 = res_sil{1,1};
            res2 = res_sil{2,1};
            res3 = res_sil{3,1};
            res1 = sort(res1);
            res2 = sort(res2);
            res3 = sort(res3);
            res = [res1';res2';res3'];
            figure(1);
            color = [0,0,1;1,0,0;1,165/255,0];
            for ii = 1:3
                plot([1:90],res(ii,:),'Color',color(ii,:),LineWidth=2);
                hold on
                temp = mean(res(ii,:));
                [~,temp_index] = min(abs(res(ii,:)-temp));
                plot(temp_index,temp,'o','Color',color(ii,:),'MarkerSize',10,LineWidth=2);
            end
            title(Datasets{1})
            xlabel('Ranked Cell')
            ylabel('Silhouette Coefficient')
            box off
            figure(2);
            res1 = res_sil{1,2};
            res2 = res_sil{2,2};
            res3 = res_sil{3,2};
            res1 = sort(res1);
            res2 = sort(res2);
            res3 = sort(res3);
            res = [res1';res2';res3'];
            for ii = 1:3
                plot([1:124],res(ii,:),'Color',color(ii,:),LineWidth=2);
                hold on
                temp = mean(res(ii,:));
                [~,temp_index] = min(abs(res(ii,:)-temp));
                plot(temp_index,temp,'o','Color',color(ii,:),'MarkerSize',10,LineWidth=2);
            end
            xlim([1,124]);
            yticks(-0.8:0.6:1);
            title(Datasets{2})
            xlabel('Ranked Cell')
            ylabel('Silhouette Coefficient')
            box off
        end
    case 'Fig 2 d'
        for fakeloop = 1
            load('Res_AIG.mat');
            subplot(1,2,1);
            A =  Res{1,1}.A;
            load('Yan_label.mat'); label_m = label;
            imagesc(A); colorbar; colormap summer;
            title('Yan, AIGS','fontsize',10 );
            plot_labelnet(label_m);
            xlabel('Cell');
            ylabel('Cell');
            subplot(1,2,2);
            load('yan_SIMLR.mat');
            load('Yan.mat');
            imagesc(S); colorbar; colormap summer;
            title('Yan, SIMLR','fontsize',10 );
            plot_labelnet(gnd);
            xlabel('Cell');
            ylabel('Cell');
        end
    case 'Fig 3 a'
        for fakeloop = 1
            Datasets_New = {...
                'Yan';
                'Goolam';
                'Deng'
                'Darmanis';
                'Usoskin';
                'Xin';
                'Muraro';
                'Lake';
                'Average';
                };
            load('Res_measure.mat');
            colors = [
                0 0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.9290 0.6940 0.1250;
                0.4940 0.1840 0.5560;
                0.4660 0.6740 0.1880;
                0.3010 0.7450 0.9330
                ];

            for i_data = 1:4
                subplot(2,5,i_data);
                ARI1 = ARI(:,i_data);
                hold on
                for ii = 1:6
                    b = bar(ii,ARI1(ii),0.75,'stacked');
                    set(b(1),'FaceColor',colors(ii,:));
                end
                xlabel(Datasets_New{i_data})
                set(gca,'XTick','')
                ylabel('ARI','fontsize',10 ); box on;  ylim([0.25,1.2]);
                if i_data == 4
                    legend('AIGS','scDHA','Seurat','dropClust','SIMLR','CIDR','fontsize',10,'NumColumns',1, 'Location', 'eastoutside'); hold on
                    legend('box', 'off');
                end
            end

            for i_data = 6:10
                subplot(2,5,i_data);
                ARI1 = ARI(:,i_data-1);
                hold on
                for ii = 1:6
                    b = bar(ii,ARI1(ii),0.75,'stacked');
                    set(b(1),'FaceColor',colors(ii,:));
                end
                xlabel(Datasets_New{i_data-1})
                set(gca,'XTick','')
                ylabel('ARI','fontsize',10 ); box on;  ylim([0.25,1.2]);
            end

        end
    case 'Fig 3 b'
        for fakeloop = 1
            load('ARI_NOR.mat'); ARI = ARI_NOR; load('colors.mat'); hold on; box on;
            plot(flip(ARI(:,1)), 1:9, '.','Color',colors(6,:),'MarkerSize',40);
            plot(flip(ARI(:,2)), 1:9, '.','Color',colors(1,:),'MarkerSize',40);

            legend('With Outlier Removing','Without Outlier Removing','fontsize',10,'NumColumns',1,'Location', 'northoutside');
            legend('box','off')
            set(gca,'YTick',1:9); ylim([0,10]);
            set(gca,'YTickLabel',{'Average','Lake','Muraro','Xin','Usoskin','Wang','Deng','Goolam','Yan'},'fontsize',10 );
            xlabel('ARI','fontsize',10 );
            xlim([0.6,1])
            set(gca, 'YDir', 'reverse'); % 旋转纵轴
            box off; % 去掉边框
            grid on; % 添加网格
        end
    case 'Fig 3 c'
        for fakeloop = 1
            load('Res_stability.mat');
            yrange = [0.56,1;0.5,1;0.6,1;0.6,1;0.73,1;0.78,1;0.81,1;0.25,1];
            figure(2);
            for ii = 1:8
                subplot(2,4,ii); hold on; box on;
                a = ARIAIG9T(:,ii);
                b = ARIDHA9T(:,ii);
                data = [a;b];
                class = [ones(size(a));2*ones(size(b))];
                boxchart(data,'GroupByColor',class);

                ylabel('ARI');
                set(gca,'XTickLabel','')
                xlabel(Datasets{ii})
                ylim(yrange(ii,:));
                if ii == 1
                    legend('AIGS','scDHA','Location', 'northoutside','NumColumns',2);
                    legend('box','off')
                end
            end
        end
    case 'Fig 4 a'
        for fakeloop = 1
            load('Res_AIG_one.mat')
            for ii = 1:8
                data = Datasets{ii,1};
                eval(['load ' data]);
                subplot(2,4,ii); hold on; box on;
                nc = Res{ii,1}.nc; ne = Res{ii,1}.num_cluster;
                idx_cell = Res{ii,1}.idx_cell;
                grps = Res{ii,1}.grps; aris = zeros(1,4);
                for jj = 1:4
                    label = gnd - min(gnd) + 1;
                    label = label(idx_cell);
                    ARI = calARI(grps(:,jj),label);
                    aris(jj) = ARI;
                end
                plot(nc:nc+3, aris,'b.','MarkerSize',20);
                plot([ne,ne],[min(aris)-0.3 aris(ne-nc+1)],'r','LineWidth',2);
                xlim([nc-0.5,nc+3.5]);
                ylim([min(aris)-0.3, 1.03]);
                ylabel('ARI','fontsize',10);
                set(gca,'linewidth',1.5);
                title(data,'fontsize',10);
            end
        end
    case 'Fig 4 b'
        for fakeloop = 1
            load('Res_Dropout');
            sps = [0.0200    0.0500    0.0800    0.1100    0.1400    0.1700    0.2000];
            c = repmat([1 2 3 4 5 6 7],1,20);
            a = Res(:,:,1);
            a = reshape(a,1,size(a,1)*size(a,2)); y = a;
            catx=repmat({'2%' '5%' '8%' '11%' '14%' '17%' '20%'},1,20);
            g(1,1)=gramm('x',catx,'y',y,'color',c);
            g(1,1).stat_boxplot('width',0.5,'dodge',0);
            g(1,1).set_title(Datasets{1});
            g(1,1).axe_property('ylim',[0 1]);
            g(1,1).set_names('x','Drop rate','y','ARI','color','Species');
            g(1,1).no_legend();
            g.draw();
        end
    case 'Fig 5 a'
        for fakeloop = 1
            load('Res_Visual_AIG.mat'); load('Ccolors.mat');load('Res_other_method.mat');
            for i_id = 1:2
                id = i_id;
                data = Res_Visual_AIG{id,1}; eval(['load ' data]);
                Y = Res_Visual_AIG{id,2};
                labels = Res_Visual_AIG{id,3};
                measilhouette = Res_Visual_AIG{id,4};
                subplot(2,5,5*(i_id-1)+1);
                for ii = 1: length(unique(labels))
                    hold on; box on;
                    scatter(Y(labels == ii,1), Y(labels == ii,2),15,'filled','MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:)); pause(0.1),
                end

                xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
                ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);

                title([Datasets{id},', Silhouette = ', num2str(round(100*measilhouette)/100)]);
                set(gca,'xticklabel',{[]});
                set(gca,'yticklabel',{[]});
                xlabel('AIGS1');
                ylabel('AIGS2');

                xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
                ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);

                for im = 1:4
                    Y = Res_other_method{id,im};
                    subplot(2,5,5*(i_id-1) + im+1);
                    for ii = 1: length(unique(gnd))
                        hold on;
                        scatter(Y(gnd == ii,1), Y(gnd == ii,2),6,'filled','MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:)); pause(0.1),
                    end
                    s = silhouette(Y,gnd); measilhouette = mean(s);
                    title([Datasets{id},', Silhouette = ', num2str(round(100*measilhouette)/100)]);

                    set(gca,'xticklabel',{[]});
                    set(gca,'yticklabel',{[]});
                    if im == 1
                        xlabel('CIDR1');
                        ylabel('CIDR2');
                    elseif im == 2
                        xlabel('scDHA1');
                        ylabel('scDHA2');
                    elseif im == 3
                        xlabel('Seurat1');
                        ylabel('Seurat2');
                    elseif im == 4
                        xlabel('SILMR1');
                        ylabel('SILMR2');
                    end
                end
            end
        end
    case 'Fig 5 b'
        for fakeloop = 1
            load("Usoskin_label_number.mat");
            load("Usoskin_label.mat");
            load('Usoskin_visualization.mat');
            load('Usoskin.mat');
            load('Ccolors.mat');
            load('Res_AIG.mat');
            out = Res{5,1};
            idx_cell = out.idx_cell;
            m = 1;
            label = Usoskin_label_name(:,m);
            label_sub = gnd_new(:,m);
            label = label(idx_cell,:);
            label_sub = label_sub(idx_cell);
            for ii = 1: length(unique(label_sub))
                hold on;
                scatter(Y(label_sub == ii,1), Y(label_sub == ii,2),6,'filled','MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:)); pause(0.1),
            end
            legend(unique(label));
            legend('box','off')
            set(gca,'xticklabel',{[]},'FontSize',6);
            set(gca,'yticklabel',{[]},'FontSize',6);
            xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
            ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);
            xlabel('AIGS1')
            ylabel('AIGS2')
            title('Subcellular level')
        end
    case 'Fig 5 c'
        for fakeloop = 1
            load("Usoskin_label_number.mat");
            load("Usoskin_label.mat");
            load('Usoskin_visualization.mat');
            load('Usoskin.mat');
            load('Ccolors.mat');
            load('Res_AIG.mat');
            out = Res{5,1};
            idx_cell = out.idx_cell;
            m = 2;
            label = Usoskin_label_name(:,m);
            label_sub = gnd_new(:,m);
            label = label(idx_cell,:);
            label_sub = label_sub(idx_cell);
            for ii = 1: length(unique(label_sub))
                hold on;
                scatter(Y(label_sub == ii,1), Y(label_sub == ii,2),6,'filled','MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:)); pause(0.1),
            end
            legend(unique(label));
            legend('box','off')
            set(gca,'xticklabel',{[]},'FontSize',6);
            set(gca,'yticklabel',{[]},'FontSize',6);
            xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
            ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);
            xlabel('AIGS1')
            ylabel('AIGS2')
            title('Subcellular level')
        end
    case 'Fig 5 d'
        for fakeloop = 1
            load('Usoskin.mat');
            load('Usoskin_visualization.mat');
            load('Usoskin_gene_name.mat');
            out = AIGS_C(fea,[]);

            idg = 1:size(fea ,1); thr_1 = filter_thr(1); thr_2 = filter_thr(2);
            absmax = max(fea,[],2); fea = fea(absmax>thr_1,:); idg = idg(absmax>thr_1);
            absvar = var(fea,[],2); [pod, i_pod] = maxk(absvar,2*num_gene_min);
            if pod(end)>thr_2
                fea = fea(absvar>thr_2,:); idg = idg(absvar>thr_2);
            else
                fea = fea(i_pod,:); idg = idg(i_pod);
            end
            label = zeros(590,1);
            neighbour = knnsearch(Y,Y,'k',10);
            fea = fea(:,out.idx_cell);
            label_name = label_name(idg,:);
            A = zeros(590);
            for ii = 1:590
                A(ii,neighbour(ii,:)) = 1;
            end
            A = max(A,A');
            bins = conncomp(graph(A));
            n_class = max(bins);
            marker_gene_name = cell(n_class,1);
            grp = zeros(590,1);
            for ii = 1:n_class
                grp(bins==ii) = ii;
            end
            idx_marker = cell(n_class,1);
            marker_num = 3;
            for i_class = 1:n_class
                grp_temp = zeros(590,1);
                grp_temp(grp==i_class) = 1;
                Marker = findMarkerGenesp(fea',grp_temp,marker_num);
                idx_marker{i_class} = Marker;
                marker_gene_name{i_class} = label_name(Marker,:);
            end
            idx_marker_gene = zeros(n_class*3,1);
            marker_gene_name_all = cell(n_class*3,1);
            for ii = 1:n_class
                idx_marker_gene((ii-1)*3+1:ii*3) = idx_marker{ii};
                marker_gene_name_all((ii-1)*3+1:ii*3) = marker_gene_name{ii};
            end
            marker_gene = {'Mir704','Gas2','Wdfy3','Prg2'};
            fig_gene_index = zeros(4,1);
            for ii = 1:27
                for jj = 1:4
                    if strcmp(marker_gene_name_all{ii},marker_gene{jj})
                        fig_gene_index(jj) = ii;
                    end
                end
            end
            load("mycolor.mat");
            for ii = 1:4
                subplot(2,2,ii);
                idx_gene = fig_gene_index(ii);
                gene = fea(idx_marker_gene(idx_gene),:);
                scatter(Y(:,1),Y(:,2),20,gene,'filled')
                brighten(-0.8)

                colormap(mycolor);
                c = colorbar;
                set(c,'Ticks',[0,4,8])
                set(c,'YTickLabel',{'0','4','8'},'FontSize',14)
                set(c,'TickDirection','out')
                c.Label.String = 'Gene Expression';
                c.Label.FontSize = 10; 
                c.LineWidth = 0.001;
                c.Position(1) = c.Position(1)+0.1;
                c.Position(3) = 0.5*c.Position(3);
                set(gca,'xtick',[],'ytick',[]);
                xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
                ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);
                set(gca,'xticklabel',{[]});
                set(gca,'yticklabel',{[]});
                set(gca,'xtick',[],'ytick',[]);
                xlabel('AIGS1')
                ylabel('AIGS2')
                title(marker_gene{ii});

            end
        end
end

function [xaxis,yaxis] = plot_graph(fea,num_cell_class,gnd)
[m,n] = size(fea);
n_class = length(num_cell_class);
var_gene = zeros(m,n_class);
for ii = 1:n_class
    var_gene(:,ii) = var(fea(:,gnd==ii),0,2);
end
xaxis = zeros(m,1);
for ii = 1:m
    xaxis(ii) = var_gene(ii,:)*num_cell_class;
end
xaxis = xaxis/n;

var_intra_class = zeros(m,n_class);
for ii = 1:m
    for jj = 1:n_class
        temp1 = find(gnd==jj);
        temp2 = setdiff([1:n],temp1);
        mean1 = mean(fea(ii,temp1));
        mean2 = mean(fea(ii,temp2));
        var_intra_class(ii,jj) = (mean1-mean2)^2;
    end
end

yaxis = zeros(m,1);
for ii = 1:m
    yaxis(ii) = var_intra_class(ii,:)*num_cell_class;
end
yaxis = yaxis/n;
end

function res = sil(fea,gnd,D)
n_class = max(gnd);
n = size(fea,2);
a = zeros(n,1);
b = zeros(n,1);
accurate_class = cell(n_class,1);
for ii = 1:n_class
    accurate_class{ii} = find(gnd == ii);
end
for ii = 1:n
    i_class = gnd(ii);
    a(ii) = mean(D(ii,accurate_class{i_class}));
    temp1 = setdiff([1:n_class],i_class);
    temp2 = zeros(n_class-1,1);
    for jj = 1:n_class-1
        temp2(jj) = mean(D(ii,accurate_class{temp1(jj)}));
    end
    b(ii) = min(temp2);
end
res = (b-a)./max(a,b);
end

function a = plot_labelnet(labels)
n_class = length(unique(labels));
a = 1; hold on;
for ii = 1:n_class
    plot(sum(labels<ii)+1,sum(labels<ii)+1:sum(labels<=ii),'r.','Markersize',6); hold on;
    plot(sum(labels<=ii),sum(labels<ii)+1:sum(labels<=ii),'r.','Markersize',6); hold on;
    plot(sum(labels<ii)+1:sum(labels<=ii),sum(labels<ii)+1,'r.','Markersize',6); hold on;
    plot(sum(labels<ii)+1:sum(labels<=ii),sum(labels<=ii),'r.','Markersize',6); hold on;
end
end