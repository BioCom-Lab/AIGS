function [isMarkerGene] = isMarker(geneExpression, geneClass)
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