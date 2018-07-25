function [parcels,mus,sigmas,features,data] = sample_synthetic(z,d,mu_0,kappa_0,nu_0,sigma_0)

labels = unique(z);
parcels = cell(length(labels),1);
mus = cell(length(labels),1);
sigmas = cell(length(labels),1);

features = cell(length(labels),1);

for l = labels
    parcels{l} = find(z == l);
    
    [mu,sigma] = sample_priors(mu_0,kappa_0,nu_0,sigma_0,d);
    
    sigmas{l} = sigma;
    mus{l} = mu;
    
    samples = mvnrnd(mus{l},diag(sigmas{l}),length(parcels{l}));
    samples = normc(samples);
    features{l} = samples;
    
end

data = zeros(length([parcels{:}]),d);

for l = labels
    data(parcels{l},:) = features{l};
end

end

function [mu,sigma] = sample_priors(mu_0,kappa_0,nu_0,sigma_0,size)

x = chi2rnd(nu_0,[size,1]);
sigma = (nu_0*sigma_0)./x;
mu = normrnd(mu_0,sigma./kappa_0);

end