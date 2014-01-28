function z = WardClustering(D, adj_list, n_clust, vox_to_clust)
if (nargin == 4)
    unclust = ~ismember(1:size(D,1),vox_to_clust);
    adj_list{unclust} = [];
    n_clust = n_clust + sum(unclust);
end


z = cluster(LinkageConstrained(D, adj_list),'maxclust',n_clust);
end