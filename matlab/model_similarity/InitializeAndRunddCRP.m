% Initializes our method using a Ward clustering linkage matrix
%   Z, a (normalized) connectivity matrix D_norm, and adjacency list defining
%   spatial adjacency, the possible numbers of parcels to consider for
%   initialization ("sizes"), hyperparameters alpha, kappa, nu, and sigsq, the
%   number of passes over the dataset MCMC should be run, the ground truth
%   parcellation gt_z (if known, empty otherwise), and a verbose flag which
%   determines whether update information is printed every 1000 iterations.
%   Returns the MAP parcellation map_z, as well as a stats objects with
%   information about the iterations of the model.

function [map_z, stats,initial_parc] = InitializeAndRunddCRP(Z, D_norm, adj_list, sizes, alpha, kappa, ... 
    nu, sigsq, pass_limit, gt_z, verbose, varargin)

% We add an optional parameter to receive the edge_prior matrix.
p = inputParser;
errorMsg = 'Value must be numeric or boolean.'; 
validationFcn = @(x) assert(isnumeric(x) || islogical(x),errorMsg);
p.addParameter('edge_prior',false,validationFcn);

p.parse(varargin{:})
edge_prior = p.Results.edge_prior;

%%% Do not touch below here. %%%

% Standard alpha = 10, kappa = 0.0001, nu = 1

% Find highest-probability greedy parcellation to initialize ddCRP
% sizes can be a list [size_1,size_2...]
logp = LogProbWC(D_norm, Z, sizes, alpha, kappa, nu, sigsq);
[~,max_i] = max(logp);
z = cluster(Z, 'maxclust', sizes(max_i));

% Construct a spanning tree within each cluster as initialization for c
% c is a vector, where each index c(i) is another index in 1:N, indicating
% parent of vertex i
c = ClusterSpanningTrees(z, adj_list);


[map_z,stats,initial_parc] = ddCRP(D_norm, adj_list, [], gt_z, ...
                  pass_limit, alpha, kappa, nu, sigsq, ...
                  5, verbose, 'edge_prior', edge_prior);
end