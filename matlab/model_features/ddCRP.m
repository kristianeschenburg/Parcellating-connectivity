% Main function: Fits our model, given a connectivity matrix D and spatial
%   adjacency specified by adj_list.  An initialization of the voxel links
%   init_c and a ground truth parcellation gt_z (for comparison) can optionally
%   be provided. MCMC will be run for mcmc_passes over the dataset, with
%   hyperparameters alpha, kappa, nu, and sigsq. Diagnostic information will
%   be saved every stats_interval iterations, and will be printed to the console
%   if verbose is True. The MAP parcellation and diagnostic information is
%   returned.
function [map_z,stats] = ddCRP(features, adj_list, init_c, gt_z, mcmc_passes, ...
                          alpha, kappa, nu, sigsq, stats_interval, ... 
                          verbose, varargin)
                      
% We add an optional parameter to receive the edge_prior matrix.
p = inputParser;
errorMsg = 'Value must be numeric or boolean.'; 
validationFcn = @(x) assert(isnumeric(x) || islogical(x),errorMsg);
p.addParameter('edge_prior',false,validationFcn);

p.parse(varargin{:})
edge_prior = p.Results.edge_prior;

hyp = [kappa,nu,sigsq];
nvox = length(adj_list);

% Generate random initialization if not specified
% c is the linkage vector i.e. for connected components 
% c(i) is parent vertex of vertex (i)
if (isempty(init_c))
    c = zeros(nvox, 1);
    for i = 1:nvox
        neighbors = [adj_list{i}; i];
        c(i) = neighbors(randi(length(neighbors)));
    end
else
    c = init_c;
end

% Initialize spatial connection matrix
% G has ones at indices ( i,c(i) )
G = sparse(1:nvox,c,1,nvox,nvox);

% Compute the connected components of G
% K = number of parcels
% z = cluster ID for each vertex
% parcels = index assignments to each cluster
[K, z, parcels] = ConnectedComp(G);

curr_lp = FullProbabilityddCRP_Data(hyp,parcels,features);


stats = struct('times',[],'lp',[],'NMI',[],'K',[], ...
                'z', zeros(0,nvox), 'c', zeros(0,nvox));
max_lp = -Inf;
t0 = tic;
steps = 0;
for pass = 1:mcmc_passes
    % Visit elements randomly
    order = randperm(nvox);
    
    % loop over randomized list of vertices
    for i = order

        if (curr_lp > max_lp)
            max_lp = curr_lp;
            map_z = z;
        end
        
        if (mod(steps, stats_interval) == 0)
            stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, ...
                                       map_z, verbose);
        end
        
        %%% Compute change in log-prob when removing the edge c_i %%%

        % remove index i's parent in sparse matrix G
        G(i,c(i)) = 0;

        %%% %%% %%% %%% %%% CONFUSED ABOUT THIS PART %%% %%% %%% %%% %%%

        % if the vertex pointed to itself
        if (c(i) == i)
            % Removing self-loop, parcellation won't change
            rem_delta_lp = -log(alpha);
            z_rem = z; parcels_rem = parcels;

       	% if vertex points to a different parent than itself
        else

        	% compute connected components, and component IDs for each vertex (z_rem), due to proposal
        	% (remember, G has the link from i-->c(i) removed)
        	% K_rem = new number of clusters after parent removal
        	% z_rem = new cluster assignment after parent removal
        	% parcels_rem = new index groups after parent removal
            [K_rem, z_rem, parcels_rem] = ConnectedComp(G);

            % compute change in likelihood due to proposal
            % removing the link might have split a cluster resulting in 2 new clusters
            if (K_rem ~= K)
                % We split a cluster, compute change in likelihood
                rem_delta_lp = -LogProbDiff(hyp,parcels_rem,z_rem(i),z_rem(c(i)),features);

            % otherwise, if we don't split a cluster
            else
                rem_delta_lp = 0;
            end
        end
        
        %%%%% Compute change in log-prob for each possible edge c_i %%%%%

        % given the current state of the parcellation:
        % K_rem, z_rem, parcels_rem (where rem = 'removed')
        % compute change in likelihood for each neighbor of vertex of interest

        % get neighbors of vertex i %%%
        adj_list_i = adj_list{i};

        % initialize log-probability vector of neighbors + 1 (for self)
        lp = zeros(length(adj_list_i)+1, 1);

        % assign self choice a probability of log(alpha)
        % here basically saying (choose to sit at table by itself)
        % vertex doesn't "see" beyond it's set of neighbors, so 
        % doesn't know if it's being pointed to
        lp(end) = log(alpha);
        cached_merge = zeros(length(adj_list_i),1);

        % loop over neighbors
        for n_ind = 1:length(adj_list_i) 

        	% get neighbor ID
            n = adj_list_i(n_ind);

            % if component ID of neighbor is same as original parent before removal
            if (z_rem(n) == z_rem(c(i)))
                % Just undoing edge removal
                lp(n_ind) = -rem_delta_lp - (c(i)==i)*log(alpha);

            % otherwise, if component ID of neighbor is different (i)
            elseif (z_rem(n) ~= z_rem(i)) 
                % Proposing merge
                % First check cache to see if this is already computed
                prev_lp = find(cached_merge == z_rem(n),1);
                if (~isempty(prev_lp))
                    lp(n_ind) = lp(prev_lp);
                else
                    % This is a novel merge, compute change in likelihood
                    lp(n_ind) = LogProbDiff(hyp,parcels_rem,z_rem(i),z_rem(n),features);
                    cached_merge(n_ind) = z_rem(n);
                end
            end
        end

        %%% %%% %%% %%% %%%
        lp = abs(lp);
        
        % Pick new edge proportional to probability
        % ChooseFromLP selects random index with probability based on lp
        new_neighbor = ChooseFromLP(lp);
        if (new_neighbor <= length(adj_list_i))
            c(i) = adj_list_i(new_neighbor);
        else
            c(i) = i;
        end
        
        % Update likelihood and parcellation
        curr_lp = curr_lp + rem_delta_lp + lp(new_neighbor);
        G(i,c(i)) = 1;
        [K, z, parcels] = ConnectedComp(G);
        steps = steps + 1;
    end
end

stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, verbose);

end

% Update diagnostic stats with time (reported relative to start time t0),current
%   log-probability, current number of clusters, current parcellation z, current
%   voxel links c, number of steps, ground truth (if available), best
%   parcellation so far (map_z). If verbose=True, also print to console.
function stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, verbose)
    stats.lp = [stats.lp curr_lp];
    stats.K = [stats.K K];
    stats.z = [stats.z; z];
    elapsed = toc(t0);
    stats.times = [stats.times elapsed];
    stats.c = [stats.c; c'];
    if (~isempty(gt_z))
        stats.NMI = [stats.NMI CalcNMI(gt_z, map_z)];
    end
    if (verbose)
        if (~isempty(gt_z))
            disp(['Step: ' num2str(steps) ...
              '  Time: ' num2str(elapsed) ...
              '  LP: ' num2str(curr_lp) ...
              '  K: ' num2str(K) ...
              ' NMI: ' num2str(stats.NMI(end))]);
        else
            disp(['Step: ' num2str(steps) ...
                  '  Time: ' num2str(elapsed) ...
                  '  LP: ' num2str(curr_lp) ...
                  '  K: ' num2str(K)]);
        end
    end
end