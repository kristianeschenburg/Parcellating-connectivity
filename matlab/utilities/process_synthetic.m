similarities = 0;

%%

utility_path = char('/mnt/parcellator/parcellation/GitHub/Parcellating-connectivity/matlab/utilities/');
require_path = char('/mnt/parcellator/parcellation/GitHub/Parcellating-connectivity/matlab/requirements/');

addpath(utility_path);
addpath(require_path);

if similarities
    feature_path = char('/mnt/parcellator/parcellation/GitHub/Parcellating-connectivity/matlab/model_similarity/');
else
    feature_path = char('/mnt/parcellator/parcellation/GitHub/Parcellating-connectivity/matlab/model_features/');
end

addpath(feature_path);

%%

alpha=0.1;
mu=0;
kappa=0.1;
nu=1;
sigsq =1;
d=5;
sizes=9;

[~, adj_list, z, coords] = GenerateSynthData('square', 0.1);


%%
[parcels,mus,sigmas,features,data] = sample_synthetic(z,d,mu,kappa,nu,sigsq);

imagesc(corr(cell2mat(features)'));
colorbar();
%%
datafile = sprintf('%stest_features.mu.0.kappa.1.nu.1.sigma.7.mat',utility_path);
save(datafile,'data','-v7.3');

%%
datafile = sprintf('%stest_features.mu.0.kappa.1.nu.1.sigma.7.mat',utility_path);
temp = load(datafile);
fn = fieldnames(temp);
features = temp.(fn{1});

%%

mu_noise = 0;
sigma_noise = 0.2';

noise = normrnd(mu_noise,sigma_noise,[size(features)]);
features_noise = features + noise;

D = corr(features');
D_norm = normr(D);

[wc,Z] = WardClustering(D_norm,adj_list,sizes);

imagesc(reshape(wc,[18,18]))

%%

alpha=20;
mu=0;
kappa=1;
nu=1;
sigsq = 0.04;
d=5;
sizes=9;

if similarities
    [map_z, stats,initial_parc] = InitializeAndRunddCRP(Z,D_norm,adj_list,sizes,alpha,kappa,nu,sigsq,30,[],true);
else
    [map_z, stats,initial_parc] = ddCRP(features,adj_list,[],[],30,alpha,kappa,nu,sigsq,30,true);
end

subplot(2,3,1);
plot(1:length(stats.lp),stats.lp)
title('Log Probability');

subplot(2,3,2)
plot(1:length(stats.K),stats.K);
title('Cluster Count');

subplot(2,3,3)
plot(1:length(stats.max_lp),stats.max_lp);
title('Max-LP');
ylim([min(stats.max_lp)-10,max(stats.max_lp)+10])

subplot(2,3,4)
imagesc(reshape(initial_parc,[18,18]));
tl = {'Initial Parcellation: '; sprintf('%i Clusters',length(unique(initial_parc)))};
title(tl)

subplot(2,3,5)
imagesc(reshape(z,[18,18]))
title('Ground Truth');

subplot(2,3,6)
imagesc(reshape(map_z,[18,18]))
tl = {'Max-LP Map: '; sprintf('%i Clusters',length(unique(map_z)))};
title(tl);
colormap(jet)


sp = sprintf('Alpha: %.1f, Kappa: %i, Nu: %.1f, Sigma: %i',alpha,kappa,nu,sigsq);
suptitle(sp)

figure(2)
labels = zeros(18,18,size(stats.z,1));
for i = 1:size(stats.z,1)
    labels(:,:,i) = reshape(stats.z(i,:),[18,18]);
end

for i = 1:30
    imagesc(labels(:,:,i));
    title(sprintf('Iter: %i, Clusters: %i',i,length(unique(stats.z(i,:)))))
    pause(0.1)
    colormap(jet)
end