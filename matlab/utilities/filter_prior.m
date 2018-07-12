function [filtered_prior] = filter_prior(adj_list,indices,prior)

indices = sort(indices);
idx_place = 1:length(indices);

mapping = containers.Map(indices,idx_place);

filtered_prior = zeros(length(indices),length(indices));

for k = 1:length(adj_list)

    if ismember(k,indices)
        
        neighbors = adj_list{k};
        
        tokeep = [];
        mapped = [];
        
        for n = neighbors
            if ismember(n,indices)
                tokeep = [tokeep,n];
                mapped = [mapped,mapping(n)];
            end
        end
        
        filtered_prior(mapping(k),mapped) = prior(k,tokeep);
    end

end

filtered_prior = sparse(filtered_prior);

end