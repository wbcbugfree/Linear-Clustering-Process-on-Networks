function element_scores = element_sim(cluster1, cluster2)

assert(numel(cluster1) == numel(cluster2));
n = numel(cluster1);
cluster1 = reshape(cluster1,1,n);
cluster2 = reshape(cluster2,1,n);
alpha = 0.9;
ppr_1 = ppr_partition(cluster1, alpha);
ppr_2 = ppr_partition(cluster2, alpha);
node_scores = 1.0 - 1.0/(2.0 * alpha) * sum(abs(ppr_1 - ppr_2));
element_scores = mean(node_scores);
end


function ppr = ppr_partition(cluster, alpha)
ppr = zeros(length(cluster),length(cluster));
for i = 1:length(unique(cluster))
    clusterlist = find(cluster == i);
    Csize = length(clusterlist);
    ppr_result = alpha/Csize * ones(Csize, Csize) + eye(Csize) * (1.0 - alpha);
    ppr(clusterlist,clusterlist) = ppr_result;
end
end