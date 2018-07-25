% Compute the lop-probability of a single cluster
% hyp : hyperparameters
% p_hyp : posterior parameters after incorporating data
% n : number of samples
function [lp] = LogProbCluster(hyp,p_hyp,n)

kappa_0 = hyp(1);
nu_0 = hyp(2);
sigma_0 = hyp(3);

kappa_n = p_hyp{1};
nu_n = p_hyp{2};
sigma_n = p_hyp{3};

[~,p] = size(sigma_n);

% ratio of gamma functions
gam = gammaln(nu_n/2) - gammaln(nu_0/2);

% inner log computations
inner = (1/2) * (log(kappa_0) + nu_0*log(nu_0*sigma_0) ...
		- log(kappa_n) - nu_n*log(nu_n) - n*log(pi));

% compute sum of sigma_n for each feature
outer = sum(log(1./sigma_n));

lp = p*(gam + inner) + outer;

end