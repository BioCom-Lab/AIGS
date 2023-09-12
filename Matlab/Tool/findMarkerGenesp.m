function markerGenes = findMarkerGenesp(data, labels,marker_num)

numGenes = size(data, 2); 
pValues = zeros(numGenes, 1); 

for i = 1:numGenes
    geneExp = data(:, i); 
    [isMarkerGene] = isMarker(geneExp, labels); 
    if isMarkerGene && sum(geneExp==0)>0.3*length(labels)
        pValues(i) = anova1(geneExp, labels, 'off'); 
    else
        pValues(i) = Inf; 
    end

end

[~, idx] = sort(pValues, 'ascend'); 

 markerGenes = idx(1:marker_num); 

end


