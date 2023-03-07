function embedding = optimize_layout_euclidean(embedding,n_epochs,a,b,A)
initial_alpha = 1;
alpha = initial_alpha;

for n = 1:n_epochs
    [embedding] = UMAP_each_point(embedding,A,a,b,alpha);
    alpha = initial_alpha*(1-n/n_epochs);
end
end