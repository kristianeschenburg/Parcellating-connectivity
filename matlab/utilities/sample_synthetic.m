function [parcels,mus,sigmas,features,data] = sample_synthetic(z,d,mu_0,kappa_0,nu_0,sigma_0)

labels = unique(z);
parcels = cell(length(labels),1);
mus = cell(length(labels),1);
sigmas = cell(length(labels),1);

features = cell(length(labels),1);

for l = labels
    parcels{l} = find(z == l);
    
    [mu_p,sigma_p] = sample_priors(mu_0,kappa_0,nu_0,sigma_0,d);
    
    sigmas{l} = sigma_p;
    mus{l} = mu_p;
    
    samples = mvnrnd(mus{l},diag(sigmas{l}),length(parcels{l}));
    samples = normc(samples);
    features{l} = samples;
    
end

disp(size(parcels))
disp(size(features))

data = zeros(length(z),d);

for p = 1:length(parcels)
    data(parcels{p},:) = features{p};
end

end

function [mu,sigma_p] = sample_priors(mu_0,kappa_0,nu_0,sigma_0,size)

x = chi2rnd(nu_0,[size,1]);
sigma_p = (nu_0*sigma_0)./x;
mu = normrnd(mu_0,sigma_p./kappa_0);

end