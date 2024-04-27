clear;clc;close all
warning off; id = 1;

Jobs = {
    1,  'Fig S1'
    2,  'Fig S2'
    3,  'Fig S3'
    4,  'Fig S4'
    5,  'Fig S5'
    6,  'Fig S6'
    7,  'Fig S7'
    8,  'Fig S8-15'
    9,  'Fig S16'
    10  'Fig S17'
    11  'Fig S18'
    12  'Fig S19'
    13  'Fig S20'
    14  'Fig S21b'
    15  'Fig S22'
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
    case 'Fig S1'
        for fakeloop = 1
            load('Res_AIG_one.mat');
            colors =[0.8,0.9,1.0;
                0.6784,0.8471,0.9020;
                0.2117,0.3921,0.5450];
            text_xlabel = {'Raw','Filtered','Selected'};
            for ii = 3:8
                subplot(2,3,ii-2);
                out = Res{ii};
                gene_number = out.num_gene;
                gene_number(3) = gene_number(3)*10;
                hold on
                for jj = 1:3
                    b = bar(jj,gene_number(jj),0.75,'stacked');
                    set(b(1),'FaceColor',colors(jj,:));
                    text(jj, gene_number(jj), num2str(gene_number(jj)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                    text(jj, -4000, text_xlabel{jj}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                end
                title(Datasets{ii})
                box off
                ylabel('Number of Genes','fontsize',10)
                set(gca,'XTickLabel','');

            end
        end
    case 'Fig S2'
        for fakeloop = 1
            point_size = 5;
            load('Res_AIG_one.mat');
            for id = 3:8
                load(Datasets{id});
                n_class = length(unique(gnd));
                num_cell_class = zeros(n_class,1);

                for ii = 1:n_class
                    num_cell_class(ii) = sum(gnd==ii);
                end
                [xaxis,yaxis] = plot_graph(fea,num_cell_class,gnd);
                subplot(2,3,id-2);
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
    case 'Fig S3'
        for fakeloop = 1
            load('Res_AIG_one.mat')
            for ii = 1:4
                load(Datasets{ii})
                fea = full(fea);
                out = Res{ii};
                idx_gene = out.idx_gene;
                Y_ag = tsne(fea(idx_gene,:)');
                sil_ag = silhouette(Y_ag,gnd);
                idg = 1:size(fea,1); thr_1 = filter_thr(1); thr_2 = filter_thr(2);
                absmax = max(fea,[],2); fea = fea(absmax>thr_1,:);
                absvar = var(fea,[],2); [pod, i_pod] = maxk(absvar,2*num_gene_min);
                if pod(end)>thr_2
                    fea = fea(absvar>thr_2,:);
                else
                    fea = fea(i_pod,:);
                end
                Y_bg = tsne(fea');
                sil_bg = silhouette(Y_bg,gnd);
                subplot(4,2,2*(ii-1)+1)
                gscatter(Y_bg(:,1),Y_bg(:,2),gnd);legend off
                xlabel('tsne1');
                ylabel('tsne2');
                temp1 = num2str(100*mean(sil_bg)/100);
                temp1 = temp1(1:4);
                title_bg = [Datasets{ii},' Before gene selection, Silhouette=', temp1];
                title(title_bg);
                subplot(4,2,2*(ii-1)+2)
                gscatter(Y_ag(:,1),Y_ag(:,2),gnd);legend off
                xlabel('tsne1');
                ylabel('tsne2');
                temp2 = num2str(100*mean(sil_ag)/100);
                temp2 = temp2(1:4);
                title_ag = [Datasets{ii},' After gene selection, Silhouette=', temp2];
                title(title_ag);
            end
        end
    case 'Fig S4'
        for fakeloop = 1
            load('Res_AIG_one.mat')
            for ii = 5:8
                load(Datasets{ii})
                fea = full(fea);
                out = Res{ii};
                idx_gene = out.idx_gene;
                Y_ag = tsne(fea(idx_gene,:)');
                sil_ag = silhouette(Y_ag,gnd);
                idg = 1:size(fea,1); thr_1 = filter_thr(1); thr_2 = filter_thr(2);
                absmax = max(fea,[],2); fea = fea(absmax>thr_1,:);
                absvar = var(fea,[],2); [pod, i_pod] = maxk(absvar,2*num_gene_min);
                if pod(end)>thr_2
                    fea = fea(absvar>thr_2,:);
                else
                    fea = fea(i_pod,:);
                end
                Y_bg = tsne(fea');
                sil_bg = silhouette(Y_bg,gnd);
                subplot(4,2,2*(ii-5)+1)
                gscatter(Y_bg(:,1),Y_bg(:,2),gnd);legend off
                xlabel('tsne1');
                ylabel('tsne2');
                temp1 = num2str(100*mean(sil_bg)/100);
                temp1 = temp1(1:4);
                title_bg = [Datasets{ii},' Before gene selection, Silhouette=', temp1];
                title(title_bg);
                subplot(4,2,2*(ii-5)+2)
                gscatter(Y_ag(:,1),Y_ag(:,2),gnd);legend off
                xlabel('tsne1');
                ylabel('tsne2');
                temp2 = num2str(100*mean(sil_ag)/100);
                temp2 = temp2(1:4);
                title_ag = [Datasets{ii},' After gene selection, Silhouette=', temp2];
                title(title_ag);
            end
        end
    case 'Fig S5'
        for fakeloop = 1
            load('Res_AIG_one.mat');
            load('Goolam.mat');
            idx_cell = Res{2,1}.idx_cell;
            label = gnd(idx_cell);
            subplot(3,2,1);
            A =  Res{2,1}.A;
            load('Goolam_label.mat'); label_m = label; [label_m,sd] = sort(label_m);
            imagesc(A(sd,sd)); colorbar; colormap summer;
            title('Goolam, AIGS','fontsize',10 );
            plot_labelnet(label_m);

            subplot(3,2,2);
            load('goolam_SIMLR.mat');
            [label,sd] = sort(gnd);
            imagesc(S(sd,sd)); colorbar; colormap summer;
            title('Goolam, SIMLR','fontsize',10 );
            plot_labelnet(gnd);

            subplot(3,2,3);
            load('Darmanis.mat');
            idx_cell = Res{4,1}.idx_cell;
            label = gnd(idx_cell);
            A =  Res{4,1}.A;
            load('Darmanis_label.mat'); label_m = label; [label_m,sd] = sort(label_m);
            imagesc(A(sd,sd)); colorbar; colormap summer;
            title('Darmanis, AIGS','fontsize',10 );
            plot_labelnet(label_m);

            subplot(3,2,4);
            load('darmanis_SIMLR.mat');
            [label,sd] = sort(gnd);
            imagesc(S(sd,sd)); colorbar; colormap summer;
            title('Darmanis, SIMLR','fontsize',10);
            plot_labelnet(gnd);

            subplot(3,2,5);
            load('Usoskin.mat');
            A =  Res{5,1}.A;
            idx_cell = Res{5,1}.idx_cell;
            label = gnd(idx_cell);
            label_m = label; [label_m,sd] = sort(label_m);
            imagesc(A(sd,sd)); colorbar; colormap summer;
            title('Usoskin, AIGS','fontsize',10 );
            plot_labelnet(label_m);

            subplot(3,2,6);
            load('usoskin_SIMLR.mat');
            [label,sd] = sort(gnd);
            imagesc(S(sd,sd)); colorbar; colormap summer;
            title('Usoskin, SIMLR','fontsize',10);
            plot_labelnet(gnd);
        end
    case 'Fig S6'
        for fakeloop = 1
            load('Res_measure.mat');
            load('colors.mat');
            subplot(6,2,1:2);
            b = bar(ACC'); set(b,'edgecolor','none');
            for i = 1:6
                set(b(i),'FaceColor',colors(i,:));
            end
            set(gca,'XTickLabel',{'Yan','Goolam','Deng','Darmanis','Usoskin','Xin','Muraro','Lake','Average'},'fontsize',10);
            legend('AIGS','scDHA','Seurat','dropClust','SIMLR','CIDR','fontsize',10,'NumColumns',6); hold on
            ylabel('ACC','fontsize',10); box on;  ylim([0.25,1.02]);
            set(gca,'linewidth',1.5);
            subplot(6,2,3:4);
            b = bar(NMI'); set(b,'edgecolor','none');
            for i = 1:6
                set(b(i),'FaceColor',colors(i,:));
            end
            set(gca,'XTickLabel',{'Yan','Goolam','Deng','Darmanis','Usoskin','Xin','Muraro','Lake','Average'},'fontsize',10);
            ylabel('NMI','fontsize',10); box on;  ylim([0.25,1.02]);
            set(gca,'linewidth',1.5);
            subplot(6,2,5:6);
            b = bar(Jaccard'); set(b,'edgecolor','none');
            for i = 1:6
                set(b(i),'FaceColor',colors(i,:));
            end
            set(gca,'XTickLabel',{'Yan','Goolam','Deng','Darmanis','Usoskin','Xin','Muraro','Lake','Average'},'fontsize',10);
            ylabel('Jaccard','fontsize',10); box on;  ylim([0.25,1.02]);
            set(gca,'linewidth',1.5);
            subplot(6,2,7:8);
            b = bar(Fmeasure'); set(b,'edgecolor','none');
            for i = 1:6
                set(b(i),'FaceColor',colors(i,:));
            end
            set(gca,'XTickLabel',{'Yan','Goolam','Deng','Darmanis','Usoskin','Xin','Muraro','Lake','Average'},'fontsize',10);
            ylabel('Fmeasure','fontsize',10); box on;  ylim([0.25,1.02]);
            set(gca,'linewidth',1.5);
            load('Res_measure.mat');
            yvalues = {'AIGS','scDHA','Seurat','dropClust','SIMLR','CIDR'};
            xvalues = {'Yan','Goolam','Deng','Darmanis','Usoskin','Xin','Muraro','Lake','Average'};
            subplot(6,2,9);
            heatmap(xvalues,yvalues,ACCORDER); colorbar off; colormap summer; grid off;
            title('ACC Rank');
            subplot(6,2,10);
            heatmap(xvalues,yvalues,NMIORDER); colorbar off; colormap summer; grid off;
            title('NMI Rank');
            subplot(6,2,11);
            heatmap(xvalues,yvalues,JaccardORDER); colorbar off; colormap summer; grid off;
            title('Jaccard Rank');
            subplot(6,2,12);
            heatmap(xvalues,yvalues,FmeasureORDER); colorbar off; colormap summer; grid off;
            title('Fmeasure Rank');
        end
    case 'Fig S7'
        for fakeloop = 1
            load('Res_WGS.mat'); load('colors.mat');

            subplot(2,1,1);hold on;
            ARI = [AT_WGS(:,1),AT(:,1)]; ARI = round(100*ARI)/100;
            b = bar(ARI);set(b,'edgecolor','none');
            for i = 1:2
                set(b(i),'FaceColor',colors(i,:));
            end
            xtips1 = b.XEndPoints;
            ytips1 = b.YEndPoints;
            labels1 = string(b(1).YData);
            text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','fontsize',10)
            xtips1 = b(2).XEndPoints;
            ytips1 = b(2).YEndPoints;
            labels2 = string(b(2).YData);
            text(xtips1,ytips1,labels2,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','fontsize',10)
            set(gca,'XTickLabel',{'Yan','Goolam','Deng','Darmanis','Usoskin','Xin','Muraro','Lake','Average'},'fontsize',10 );
            legend('Without Gene Selection','With Gene Selection','fontsize',10,'NumColumns',6 );
            ylabel('ARI','fontsize',10 ); box on;  ylim([0.25,1.1]);
            set(gca,'linewidth',1.5);

            subplot(2,1,2);hold on; ylim([0,20]);
            x = [1,2]; ARI = [AT_WGS(:,2),AT(:,2)]; ARI = round(10*ARI)/10;
            b = bar(ARI);set(b,'edgecolor','none');
            for i = 1:2
                set(b(i),'FaceColor',colors(i,:));
            end
            xtips1 = b.XEndPoints;
            ytips1 = b.YEndPoints;
            labels1 = string(b(1).YData);
            text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','fontsize',10)
            xtips1 = b(2).XEndPoints;
            ytips1 = b(2).YEndPoints;
            labels2 = string(b(2).YData);
            text(xtips1,ytips1,labels2,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','fontsize',10)
            set(gca,'XTickLabel',{'Yan','Goolam','Deng','Darmanis','Usoskin','Xin','Muraro','Lake','Average'},'fontsize',10 );
            legend('Without Gene Selection','With Gene Selection','fontsize',10,'NumColumns',6 );
            ylabel('CPU Time (s)','fontsize',10 ); box on;
            set(gca,'linewidth',1.5);

            exportgraphics(figure(1),'FigA8.pdf','ContentType','vector');
        end
    case 'Fig S8-15'
        for fakeloop = 1
            figure(1);
            load('Yan_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Yan');

            figure(2);
            load('Goolam_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Goolam');

            figure(3);
            load('Deng_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Deng');

            figure(4);
            load('Darmanis_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Darmanis');

            figure(5);
            load('Usoskin_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Usoskin');

            figure(6);
            load('Xin_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Xin');

            figure(7);
            load('Muraro_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Muraro');

            figure(8);
            load('Lake_label.mat');
            List = label2sankey(label,grp);
            colorList=[0.4600    0.5400    0.4600
                0.5400    0.6800    0.4600
                0.4100    0.4900    0.3600
                0.3800    0.5300    0.8400
                0.4400    0.5900    0.8700
                0.5800    0.7900    0.9300
                0.6500    0.6400    0.8400
                0.6300    0.6300    0.8000
                0.5600    0.5300    0.6700
                0.7600    0.8100    0.4300
                0.5600    0.8600    0.9700
                0.7800    0.5900    0.6500
                0.8900    0.9100    0.5300
                0.9300    0.5600    0.2500];
            sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
            axis off; box on;
            title('Lake');

        end
    case 'Fig S16'
        for fakeloop = 1
            load('Res_Dropout');
            sps = [0.0200    0.0500    0.0800    0.1100    0.1400    0.1700    0.2000];
            for id = 1:8
                c = repmat([1 2 3 4 5 6 7],1,20);
                a = Res(:,:,id);
                a = reshape(a,1,size(a,1)*size(a,2)); y = a;
                catx=repmat({'2%' '5%' '8%' '11%' '14%' '17%' '20%'},1,20);
                if mod(id,2)==1
                    g(ceil(id/2),1)=gramm('x',catx,'y',y,'color',c);
                    g(ceil(id/2),1).stat_boxplot('width',0.5,'dodge',0);
                    g(ceil(id/2),1).set_title(Datasets{id}(1:end-1));
                    g(ceil(id/2),1).axe_property('ylim',[0 1])
                    g(ceil(id/2),1).set_names('x','Drop rate','y','ARI','color','Species');
                    g(ceil(id/2),1).no_legend();
                else
                    g(ceil(id/2),2)=gramm('x',catx,'y',y,'color',c);
                    g(ceil(id/2),2).stat_boxplot('width',0.5,'dodge',0);
                    g(ceil(id/2),2).set_title(Datasets{id}(1:end-1));
                    g(ceil(id/2),2).axe_property('ylim',[0 1])
                    g(ceil(id/2),2).set_names('x','Drop rate','y','ARI','color','Species');
                    g(ceil(id/2),2).no_legend();
                end
            end
            g.draw();
        end
    case 'Fig S17'
        for fakeloop = 1
            load('Res_Visual_AIG.mat'); load('Ccolors.mat');load('Res_other_method.mat');
            for i_id = 3:8
                id = i_id;
                data = Res_Visual_AIG{id,1}; eval(['load ' data]);
                Y = Res_Visual_AIG{id,2};
                labels = Res_Visual_AIG{id,3};
                measilhouette = Res_Visual_AIG{id,4};
                subplot(6,5,5*(i_id-3)+1);
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
                    subplot(6,5,5*(i_id-3) + im+1);
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
    case 'Fig S18'
        for fakeloop = 1
            load('Res_Visual_AIG.mat'); load('colors.mat'); load('Res_other_method.mat');
            for i_id = 1:8
                id = i_id; Ss = cell(1,5);
                data = Res_Visual_AIG{id,1}; eval(['load ' data]);
                Y = Res_Visual_AIG{id,2};
                labels = Res_Visual_AIG{id,3};
                Ss{1,1} = silhouette(Y,labels);

                for im = 1:4
                    Y = Res_other_method{id,im};
                    s = silhouette(Y,gnd);
                    Ss{1,im+1} = s;

                end
                subplot(4,2,i_id)
                violin(Ss,'facecolor',colors(2,:),'medc','');
                legend off;
                set(gca,'xticklabel',{'scAIG','CIDR', 'scDHA', 'Seurat','SIMLR'});
                ylabel('silhouette');
                title(Datasets{id});
            end
        end
    case 'Fig S19'
        for fakeloop = 1
            load('Res_Visual_AIG.mat');
            load('Ccolors.mat');
            load('Res_embedding.mat')
            for i_id = 1:4
                Y_AIGS = Res_Visual_AIG{i_id,2};
                Res_embedding{i_id,4} = Y_AIGS;
            end
            Title = {'Initial','tsne','UMAP','AIGS'};
            for i_id = 1:4
                labels = Res_Visual_AIG{i_id,3};
                for f_id = 1:4
                    subplot(4,4,(i_id-1)*4+f_id)
                    Y = Res_embedding{i_id,f_id};
                    sil = mean(silhouette(Y,labels));
                    sil = num2str(sil);
                    sil = sil(1:4);
                    for ii = 1: length(unique(labels))
                        hold on; box on;
                        scatter(Y(labels == ii,1), Y(labels == ii,2),15,'filled','MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:)); pause(0.1),
                    end
                    xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
                    ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);
                    title([Datasets{i_id},', Silhouette = ', sil]);
                    set(gca,'xticklabel',{[]});
                    set(gca,'yticklabel',{[]});
                    xlabel([Title{f_id}, '1']);
                    ylabel([Title{f_id}, '2']);
                end
            end
        end
    case 'Fig S20'
        for fakeloop = 1
            load('Res_Visual_AIG.mat');
            load('Ccolors.mat');
            load('Res_embedding.mat')
            for i_id = 5:8
                Y_AIGS = Res_Visual_AIG{i_id,2};
                Res_embedding{i_id,4} = Y_AIGS;
            end
            Title = {'Initial','tsne','UMAP','AIGS'};
            for i_id = 5:8
                labels = Res_Visual_AIG{i_id,3};
                for f_id = 1:4
                    subplot(4,4,(i_id-5)*4+f_id)
                    Y = Res_embedding{i_id,f_id};
                    sil = mean(silhouette(Y,labels));
                    sil = num2str(sil);
                    sil = sil(1:4);
                    for ii = 1: length(unique(labels))
                        hold on; box on;
                        scatter(Y(labels == ii,1), Y(labels == ii,2),15,'filled','MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:)); pause(0.1),
                    end
                    xlim([min(Y(:,1))-0.05,max(Y(:,1))+0.05]);
                    ylim([min(Y(:,2))-0.05,max(Y(:,2))+0.05]);
                    title([Datasets{i_id},', Silhouette = ', sil]);
                    set(gca,'xticklabel',{[]});
                    set(gca,'yticklabel',{[]});
                    xlabel([Title{f_id}, '1']);
                    ylabel([Title{f_id}, '2']);
                end
            end
        end
    case 'Fig S21b'
        for fakeloop = 1
            load('res_for_mean_deng.mat');
            load('Deng_markergene_name.mat');
            load('Deng_label.mat');
            label = cell(268,2);
            for ii = 1:268
                if strcmp(label{ii,1},'mid2cell')
                    label{ii,2} = 'mid2cell';
                elseif strcmp(label{ii,1},'late2cell')
                    label{ii,2} = 'late2cell';
                end
            end
            class2class = {'16cell',1;'4cell',3;'8cell',1;'blast',6;'late2cell',4;'mid2cell',5;'zygote',2};
            Class_name = {'zygote','mid2cell','late2cell','4cell','8cell,16cell','blast'};
            order = [2,5,4,3,1,6];
            grp_new = zeros(254,1);

            for jj = 1:6
                grp_new(out.grp==order(jj))=jj;
            end
            grp = grp_new;
            num_gene = 3;
            Deng_name = cell(num_gene*6,1);

            for ii = 1:6
                Deng_name((ii-1)*num_gene+1:(ii-1)*num_gene+num_gene) = label_name{order(ii)};
            end
            for ii = 1:6
                A = res{order(ii)};
                A = A(:,order);
                A_all((ii-1)*num_gene+1:(ii-1)*num_gene+num_gene,:) = A;
            end
            A_all = A_all';
            x_axis = [1:2:20];
            y_axis = [num_gene*6:-1:1];
            color = zeros(4,3);
            color(1,:) = [0.6350 0.0780 0.1840];
            color(2,:) = [0.8500 0.3250 0.0980];
            color(3,:) = [0.9290 0.6940 0.1250];
            color(4,:) = [0.3010 0.7450 0.9330];
            color(5,:) = [0,0,0];
            temp = A_all(:);
            temp = sort(temp,'descend');
            a1 = [19,40,60,80,108];
            % a1 = a1(1:end);
            temp = temp(a1);
            color_matrix = zeros(size(A_all));
            size_matrix = zeros(size(A_all));
            for ii = 1:6
                for jj = 1:num_gene*6
                    if A_all(ii,jj)>=temp(1)
                        color_matrix(ii,jj) = 1;
                        size_matrix(ii,jj) = 50;
                    elseif A_all(ii,jj)>=temp(2)
                        color_matrix(ii,jj) = 2;
                        size_matrix(ii,jj) = 40;
                    elseif A_all(ii,jj)>=temp(3)
                        color_matrix(ii,jj) = 3;
                        size_matrix(ii,jj) = 30;
                    elseif A_all(ii,jj)>=temp(4)
                        color_matrix(ii,jj) = 4;
                        size_matrix(ii,jj) = 20;
                    else
                        color_matrix(ii,jj) = 5;
                        size_matrix(ii,jj) = 10;
                    end
                end
            end
            for ii = 1:6
                y = linspace(-0.5,1+num_gene*6,1000);
                plot(x_axis(ii)*ones(size(y)),y,'LineWidth',2,'Color',0.8275*ones(3,1));
                hold on
            end
            for ii = 1:num_gene*6
                x = linspace(0.5,11.5,1000);
                plot(x,ii*ones(size(x)),'LineWidth',2,'Color',0.8275*ones(3,1));
            end
            for ii = 1:6
                for jj = 1:num_gene*6
                    if color_matrix(ii,jj) == 1
                        p1 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    elseif color_matrix(ii,jj) == 2
                        p2 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    elseif color_matrix(ii,jj) == 3
                        p3 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    elseif color_matrix(ii,jj) == 4
                        p4 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    else
                        p5 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    end
                end
            end
            for ii = 1:num_gene*6
              text(-2,1+num_gene*6-ii,Deng_name{ii},'fontsize',10)
            end
            for ii = 1:6
              text(x_axis(ii)-0.5,20,Class_name{ii},'fontsize',10)
            end
            legend([p5,p4,p3,p2,p1],'0%-20%','20%-40%','40%-60%','60%-80%','80%-100%','NumColumns',5,'FontName','Times New Roman','FontSize',15)
            legend box off
            box off
            axis off
        end
    case 'Fig S22'
        for fakeloop = 1
            load('res_for_mean_yan.mat');
            load('Yan_markergene_name.mat');
            load('Yan.mat')
            Yan_gene_name = cell(15,1);
            for ii = 1:5
                for jj = 1:3
                    Yan_gene_name{(ii-1)*3+jj} = label_name{ii}{jj};
                end
            end
            class2class = {'zygote',1;'2cell',1;'4cell',5;'8cell',4;'16cell',2;'blast',3};
            Class_name = {'zygote,2cell','4cell','8cell','16cell','blast'};
            order = [1,5,4,2,3];
            grp_new = zeros(85,1);
            for jj = 1:5
                grp_new(grp==order(jj))=jj;
            end
            grp = grp_new;
            num_gene = 3;
            Deng_name = cell(num_gene*5,1);
            for ii = 1:5
                Deng_name((ii-1)*num_gene+1:(ii-1)*num_gene+num_gene) = label_name{order(ii)};
            end
            for ii = 1:5
                A = res{order(ii)};
                A = A(:,order);
                A_all((ii-1)*num_gene+1:(ii-1)*num_gene+num_gene,:) = A;
            end
            A_all = A_all';
            x_axis = [1:2:20];
            y_axis = [num_gene*5:-1:1];
            color = zeros(4,3);
            color(1,:) = [0.6350 0.0780 0.1840];
            color(2,:) = [0.8500 0.3250 0.0980];
            color(3,:) = [0.9290 0.6940 0.1250];
            color(4,:) = [0.3010 0.7450 0.9330];
            color(5,:) = [0,0,0];
            temp = A_all(:);
            temp = sort(temp,'descend');
            a1 = [15,30,45,60,75];

            temp = temp(a1);
            color_matrix = zeros(size(A_all));
            size_matrix = zeros(size(A_all));
            for ii = 1:5
                for jj = 1:num_gene*5
                    if A_all(ii,jj)>=temp(1)
                        color_matrix(ii,jj) = 1;
                        size_matrix(ii,jj) = 50;
                    elseif A_all(ii,jj)>=temp(2)
                        color_matrix(ii,jj) = 2;
                        size_matrix(ii,jj) = 40;
                    elseif A_all(ii,jj)>=temp(3)
                        color_matrix(ii,jj) = 3;
                        size_matrix(ii,jj) = 30;
                    elseif A_all(ii,jj)>=temp(4)
                        color_matrix(ii,jj) = 4;
                        size_matrix(ii,jj) = 20;
                    else
                        color_matrix(ii,jj) = 5;
                        size_matrix(ii,jj) = 10;
                    end
                end
            end
            for ii = 1:5
                y = linspace(-0.5,1+num_gene*5,1000);
                plot(x_axis(ii)*ones(size(y)),y,'LineWidth',2,'Color',0.8275*ones(3,1));
                hold on
            end
            for ii = 1:num_gene*5
                x = linspace(0,9.5,1000);
                plot(x,ii*ones(size(x)),'LineWidth',2,'Color',0.8275*ones(3,1));
            end
            for ii = 1:5
                for jj = 1:num_gene*5
                    if color_matrix(ii,jj) == 1
                        p1 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    elseif color_matrix(ii,jj) == 2
                        p2 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    elseif color_matrix(ii,jj) == 3
                        p3 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    elseif color_matrix(ii,jj) == 4
                        p4 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    else
                        p5 = plot(x_axis(ii),y_axis(jj),'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                    end
                    %         plot(x_axis(ii),jj,'.','MarkerSize',size_matrix(ii,jj),'Color',color(color_matrix(ii,jj),:),'LineWidth',5);
                end
            end
            for ii = 1:num_gene*5
                text(0,1+num_gene*5-ii,Yan_gene_name{ii},'fontsize',10)
            end
            for ii = 1:5
              text(x_axis(ii)-0.5,17,Class_name{ii},'fontsize',10)
            end
            legend([p5,p4,p3,p2,p1],'0%-20%','20%-40%','40%-60%','60%-80%','80%-100%','NumColumns',5)
            legend box off
            box off
            axis off
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