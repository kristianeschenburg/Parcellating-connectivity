addpath('/Users/kristianeschenburg/Documents/Code/Parcellating-connectivity/matlab/model_similarity/');
addpath('/Users/kristianeschenburg/Documents/Code/Parcellating-connectivity/matlab/utilities/');
addpath('/Users/kristianeschenburg/Documents/Code/Parcellating-connectivity/matlab/model_features/');
%%

datadir = char('/Users/kristianeschenburg/Desktop/Research/Data/');
surfacefile = sprintf('%sSurfaces/285345.L.inflated.32k_fs_LR.surf.gii',datadir);
labelfile = sprintf('%sLabels/Desikan/285345.L.aparc.32k_fs_LR.label.gii',datadir);

dataext = char('285345.L.Cortical.Regionalized.noAdj.ProbTrackX2.LogTransformed.Single.aparc.a2009s.mat');
datafile = sprintf('%sCorticalRegionalization/Destrieux/ProbTrackX2/%s',datadir,dataext);

%%

surface = gifti(surfacefile);
label = gifti(labelfile);

temp = load(datafile);
fn = fieldnames(temp);
data = temp.(fn{1});

%%

fprintf('Building adjacency list.\n');
adjacency_list = generate_adjacency(surface,false);

regions = {'L_inferiorparietal'; 'L_supramarginal'};
labelmap = containers.Map(label.labels.name,label.labels.key);
indices = [];

regmap = containers.Map();
regdata = containers.Map();

for r = regions'
    temp = find(label.cdata == labelmap(char(r)));
    regmap(r{1}) = sort(temp);
    regdata(r{1}) = data(sort(temp),:);
    indices = [indices;find(label.cdata == labelmap(char(r)))];
end
indices = sort(indices);

fprintf('Filtering adjacency list.\n');
adj_list = filter_adjacency(adjacency_list,indices);

downsampled = data(indices,:);
    
D = corr(downsampled');
D = D - diag(diag(D));

D_norm = normr(D);

%%
alpha = 0.1;
kappa = 1;
nu = 2;
sigsq = 5;

hyp = [kappa,nu,sigsq];
%%

data_ipl = data(regmap('L_inferiorparietal'),:);
[n_ipl,~] = size(data_ipl);

data_spg = data(regmap('L_supramarginal'),:);
[n_spg,~] = size(data_spg);

suff_ipl = suff_stats(data_ipl);
suff_spg = suff_stats(data_spg);

phyp_ipl = marginal_parameters(hyp,suff_ipl);
phyp_spg = marginal_parameters(hyp,suff_spg);

a = LogProbCluster(hyp,phyp_ipl,n_ipl) + LogProbCluster(hyp,phyp_spg,n_spg);
b = LogProbDiff(hyp,parcel_split,1,2,data);
c = log(alpha);

lp = [a,b,c];


%%

[z,Z] = WardClustering(D_norm,adj_list,7);

%%

c = ClusterSpanningTrees(z, adj_list);

hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq);
nvox = length(adj_list);
%%

G = sparse(1:nvox,c,1,nvox,nvox);
[K_rem, z_rem] = graphconncomp(G, 'Weak', true);
[sorted_z, sorted_i] = sort(z);
parcels_rem = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (K+1)]))));
%%