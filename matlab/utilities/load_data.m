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
kappa = 200;
nu = 1;
sigsq = 200;

sizes = [7];

hyp = [kappa,nu,sigsq];

%%

[z,Z] = WardClustering(D_norm,adj_list,7);

%%

[map_z,stats] = InitializeAndRunddCRP(Z, D_norm, downsampled, adj_list, sizes, ...
    alpha, kappa, nu, sigsq, 50, [], true);
