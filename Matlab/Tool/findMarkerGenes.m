function markerGenes = findMarkerGenes(data, labels,marker_num)

    numGenes = size(data, 2); 
    pValues = zeros(numGenes, 1);
    
    for i = 1:numGenes
        geneExp = data(:, i); 
        isMarkerGene = isMarker(geneExp, labels); 
        if isMarkerGene && sum(geneExp==0)>0.3*length(labels)
            pValues(i) = anova1(geneExp, labels, 'off'); 
        else
            pValues(i) = Inf; 
        end
    end
    
    [~, idx] = sort(pValues, 'ascend'); 
    markerGenes = idx(1:marker_num); 
end

function isMarkerGene = isMarker(geneExpression, geneClass)

    group1 = geneExpression(geneClass == 0);
    group2 = geneExpression(geneClass == 1);
    
    groupVar = [var(group1), var(group2)];
    
    meanDiff = abs(mean(group1) - mean(group2));
    
    p = anova1(geneExpression, geneClass, 'off');
    if p < 0.05 && meanDiff > 2*sqrt(max(groupVar))
        isMarkerGene = true;
    else
        isMarkerGene = false;
    end
end
