addpath('/mnt/parcellator/parcellation/GitHub/Parcellating-connectivity/matlab/model_features/');
addpath('/mnt/parcellator/parcellation/GitHub/Parcellating-connectivity/matlab/utilities/');

alpha=0.01;
mu=0;
kappa=5;
nu=1;
sigma = 7;
d=5;

sizes=9;

[~, adj_list, z, coords] = GenerateSynthData('square', 0.1);
[parcels,mus,sigmas,features,data] = sample_synthetic(z,d,mu,kappa,nu,sigma);

D = corr(data');
D_norm = normr(D);


[wc,Z] = WardClustering(D_norm,adj_list,sizes);
[map_z, stats] = InitializeAndRunddCRP(Z, D_norm, data, adj_list, ...
                    sizes, alpha, kappa, nu, sigma, ...
                    200, [], true);

subplot(2,2,1);
plot(1:length(stats.lp),stats.lp)
title('Log Probability');

subplot(2,2,2)
plot(1:length(stats.K),stats.K);
title('Cluster Count');

subplot(2,2,3)
imagesc(reshape(z,[18,18]))
title('Ground Truth');

subplot(2,2,4)
imagesc(reshape(map_z,[18,18]))
tl = sprintf('Max-LP Map: %i Clusters',length(unique(map_z)));
title(tl);
colormap(jet)

figtitle = sprintf('alpha.%.2f.kappa.%i.nu.%i.sigma.%i.d.%i.fig',alpha,kappa,nu,sigma,d);
savefig(figtitle);

%%