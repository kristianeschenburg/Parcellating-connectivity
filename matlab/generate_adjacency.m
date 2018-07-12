function [adjList] = generate_adjacency(Surf,identity)

faces = Surf.faces;

adjList=cell(length(Surf.vertices),1);

for J = 1:length(Surf.vertices)
    
    inds = find(faces == J);
    
    adj = [];
    for k = 1:length(inds)
        [r,~] = ind2sub(size(faces),inds(k));
        adj = [adj,faces(r,:)];
    end
    
    if identity
        adj = unique(adj);
    else
        adj = setdiff(unique(adj),J);
    end
    
    adjList{J} = adj;

end
