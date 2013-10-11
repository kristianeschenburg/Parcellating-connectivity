function [map_z stats] = ddCRP(subject, experiment, num_passes, ...
                               alpha, kappa, nu, sigsq, plot_interval)

hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq);
stats = struct('times',[],'lp',[],'NMI',[],'K',[], ...
               'conn_diff', zeros(0,4), 'gt_lp', []);

loaded = load(['../data/' subject '/' experiment '.mat']);
D = loaded.D;
adj_list = loaded.adj_list;
coords = loaded.coords;
nvox = size(coords, 1);

if (strcmp(experiment, 'PPA'))
    const_c = loaded.const_c;
    labels = loaded.labels;
    bold = loaded.bold;
else
    const_c = zeros(nvox, 1);
end
if (strcmp(subject, 'synth'))
    gt_z = loaded.z;
else
    gt_z = [];
end
clear loaded;


c = const_c;
for i = find(const_c==0)'
    c(i) = adj_list{i}(randi(length(adj_list{i})));
end

G = sparse(1:nvox,c,1,nvox,nvox);
[K, z, parcels] = ConnectedComp(G);
curr_lp = FullProbabilityddCRP(D, c, parcels, alpha, hyp);

max_lp = -Inf;
t0 = cputime;
for pass = 1:num_passes
    nonconst_vox = find(const_c==0);
    order = nonconst_vox(randperm(length(nonconst_vox)))';
    
    for i = order
        if (curr_lp > max_lp)
            max_lp = curr_lp;
            map_z = z;
        end
        
        stats.times = [stats.times (cputime-t0)];
        stats.lp = [stats.lp curr_lp];
        stats.K = [stats.K K];
        if (~isempty(gt_z))
            stats.NMI = [stats.NMI CalcNMI(gt_z, map_z)];
            if (abs(stats.NMI(end)-1)<10^(-8))
                %stats.gt_lp = 
                return;
            end
        end
        
        if (strcmp(experiment,'PPA'))
            stats.conn_diff = [conn_diff; ...
                               CalcPPAConnDiff(z, labels, coords, bold)];
        end
        
        if (plot_interval > 0 && ...
                              mod(length(stats.times), plot_interval) == 1)
            if (~isempty(gt_z))
                subplot(3,1,1);
                plot(stats.times,stats.lp); title('Log Prob');
                subplot(3,1,2);
                plot(stats.times,stats.K); title('Num Supervoxels');
                subplot(3,1,3);
                plot(stats.times,stats.NMI); title('NMI');
            elseif (strcmp(experiment,'PPA'))
                subplot(3,1,1);
                plot(stats.times,stats.lp); title('Log Prob');
                subplot(3,1,2);
                plot(stats.times,stats.K); title('Num Supervoxels');
                subplot(3,1,3);
                plot(stats.times,stats.conn_diff); title('Conn Diff');
                legend('LOC','TOS','RSC','IPL');
            else
                subplot(2,1,1);
                plot(stats.times,stats.lp); title('Log Prob');
                subplot(2,1,2);
                plot(stats.times,stats.K); title('Num Supervoxels');
            end
            pause(0.1);
                
        end
        if (c(i) == i)
            rem_delta_lp = -log(alpha);
            K_rem = K; z_rem = z; parcels_rem = parcels;
        else
            G(i,c(i)) = 0;
            [K_rem, z_rem, parcels_rem] = ConnectedComp(G);
            if (K_rem ~= K)
                rem_delta_lp = -LikelihoodDiff(D, ...
                                     parcels_rem, z_rem(i), z_rem(c(i)),...
                                     parcels, z(i), hyp);
            else
                rem_delta_lp = 0;
            end
        end
        
        adj_list_i = adj_list{i};
        lp = zeros(length(adj_list_i)+1, 1);
        lp(end) = log(alpha);
        for n_ind = 1:length(adj_list_i)
            n = adj_list_i(n_ind);
            c_n = c;
            c_n(i) = n;
            G_n = sparse(1:nvox,c_n,1,nvox,nvox);
            [K_n, z_n, parcels_n] = ConnectedComp(G_n);
            if (K_n == K_rem)
                continue;
            end
            if (all(z == z_n))
                % We added an edge equivalent to the one removed
                lp(n_ind) = -rem_delta_lp;
                continue;
            end
            lp(n_ind) = LikelihoodDiff(D, ...
                                       parcels_rem, z_rem(i), z_rem(n), ...
                                       parcels_n, z_n(i), hyp);
        end
        
        new_neighbor = ChooseFromLP(lp);
        if (new_neighbor <= length(adj_list_i))
            c(i) = adj_list_i(new_neighbor);
        else
            c(i) = i;
        end
        curr_lp = curr_lp + rem_delta_lp + lp(new_neighbor);
        G(i,c(i)) = 1;
        [K, z, parcels] = ConnectedComp(G);
    end
end
end

function [K, z, parcels] = ConnectedComp(G)
    [K, z] = graphconncomp(G, 'Weak', true);
    [sorted_z, sorted_i] = sort(z);
    parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (K+1)]))));
end

function ld = LikelihoodDiff(D, parcels_split, split_i1, split_i2, ...
                                parcels_merge, merge_i, hyp)
    
    split_stats = zeros(2*length(parcels_split)-1, 3);
    j = 1;
    for i = 1:length(parcels_split)
        samples = D(parcels_split{i}, parcels_split{split_i1});
        if (i == split_i1)
            samples = samples(logical(triu(ones(size(samples)),1)));
            if (isempty(samples))
                continue;
            end
        else
            samples = samples(:);
        end
        split_stats(j,1) = length(samples);
        split_stats(j,2) = sum(samples)/split_stats(j,1);
        split_stats(j,3) = sum((samples-split_stats(j,2)).^2);
        j = j+1;
    end
    for i = 1:length(parcels_split)
        samples = D(parcels_split{i}, parcels_split{split_i2});
        if (i == split_i1)
            continue;
        end
        if (i == split_i2)
            samples = samples(logical(triu(ones(size(samples)),1)));
            if (isempty(samples))
                continue;
            end
        else
            samples = samples(:);
        end
        split_stats(j,1) = length(samples);
        split_stats(j,2) = sum(samples)/split_stats(j,1);
        split_stats(j,3) = sum((samples-split_stats(j,2)).^2);
        j = j+1;
    end
    split_ll = LogLikelihood(split_stats, hyp);
    
    merge_stats = zeros(length(parcels_merge), 3);
    for i = 1:length(parcels_merge)
        samples = D(parcels_merge{i}, parcels_merge{merge_i});
        if (i == merge_i)
            samples = samples(logical(triu(ones(size(samples)),1)));
            if (isempty(samples))
                continue;
            end
        else
            samples = samples(:);
        end
        merge_stats(i,1) = length(samples);
        merge_stats(i,2) = sum(samples)/merge_stats(i,1);
        merge_stats(i,3) = sum((samples-merge_stats(i,2)).^2);
    end
    merge_ll = LogLikelihood(merge_stats, hyp);
    
    ld = merge_ll - split_ll;
end
