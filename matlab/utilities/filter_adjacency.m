function [region_list] = filter_adjacency(adj_list,indices)

indices = sort(indices);
idx_place = 1:length(indices);

mapping = containers.Map(indices,idx_place);
region_list = cell(length(indices),1);

for k = 1:length(adj_list)

    if ismember(k,indices)
        
        neighbors = adj_list{k};
        
        mapped = [];
        
        for n = neighbors
            if ismember(n,indices)
                mapped = [mapped;mapping(n)];
            end
        end
        
        region_list{mapping(k)} = mapped;
    end

end

end