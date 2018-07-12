% Initializes our method using a Ward clustering linkage matrix
%   Z, a (normalized) connectivity matrix D_norm, and adjacency list defining
%   spatial adjacency, the possible numbers of parcels to consider for
%   initialization ("sizes"), hyperparameters alpha, kappa, nu, and sigsq, the
%   number of passes over the dataset MCMC should be run, the ground truth
%   parcellation gt_z (if known, empty otherwise), and a verbose flag which
%   determines whether update information is printed every 1000 iterations.
%   Returns the MAP parcellation map_z, as well as a stats objects with
%   information about the iterations of the model
<<<<<<< HEAD
function [map_z, stats] = InitializeAndRunddCRP(Z, D_norm, adj_list, sizes, alpha, kappa, ... 
    nu, sigsq, pass_limit, gt_z, verbose, varargin)

% We add an optional parameter to receive the edge_prior matrix.
p = inputParser;
errorMsg = 'Value must be numeric or boolean.'; 
validationFcn = @(x) assert(isnumeric(x) || islogical(x),errorMsg);
p.addParameter('edge_prior',false,validationFcn);

p.parse(varargin{:})
edge_prior = p.Results.edge_prior;
=======
function [map_z, stats] = InitializeAndRunddCRP(Z, D_norm, adj_list, sizes, alpha, kappa, nu, sigsq, pass_limit, gt_z, verbose)
>>>>>>> 503340db92dcf942a0b27eae846dcd176ecdc1c5

% Standard alpha = 10, kappa = 0.0001, nu = 1

% Find highest-probability greedy parcellation to initialize
logp = LogProbWC(D_norm, Z, sizes, alpha, kappa, nu, sigsq);
[~,max_i] = max(logp);
z = cluster(Z, 'maxclust', sizes(max_i));

% Construct a spanning tree within each cluster as initialization for c
c = ClusterSpanningTrees(z, adj_list);
[map_z,stats] = ddCRP(D_norm, adj_list, c, gt_z, ...
                  pass_limit, alpha, kappa, nu, sigsq, ...
                  1000, verbose, 'edge_prior',edge_prior);
end