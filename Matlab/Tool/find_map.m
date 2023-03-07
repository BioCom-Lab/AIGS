function [map] = find_map(idx,labels)
    m = max(idx); n = max(labels);
    dist_idx = cell(m,1); dist_labels = cell(n,1);
    for ii = 1:m
        dist_idx{ii} = find(idx==ii);
    end
    for ii = 1:n
        dist_labels{ii} = find(labels==ii);
    end
    A = zeros(m,n);
    for ii = 1:m
        for jj = 1:n
            A(ii,jj) = length(intersect(dist_idx{ii},dist_labels{jj}));
        end
    end
    [~,map] = max(A,[],2);
end