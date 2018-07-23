% Given a sparse adjacency matrix G, returns the number of (undirected)
%   connected components K, the component labels z, and a cell list of
%   arrays with the indices of elements in each component

function [K, z, parcels] = ConnectedComp(G)
    
    [K, z] = graphconncomp(G, 'Weak', true);
    [sorted_z, sorted_i] = sort(z);
    parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (K+1)]))));

end