%%%
% Compute the marginal likelihood of a single cluster.
% hyp : prior hyperparameters [kappa_0,nu_0,sigma_0]
% p_hyp : hyperparameters updated with data {kappa_n; nu_n; sigma_n}, 
% length(sigma_n) = dimensionality

% of data
% n : data samples

%%%

function [z] = marginal_Z(hyp,p_hyp,n)

% ratio of gamma functions
gam = gamma(p_hyp{2}/2) / gamma(hyp(2)/2);

% ratio of mean variances
kappa = sqrt(hyp(1)/p_hyp{1});

% ratio of variance means
nu = ((hyp(2)*hyp(3))^(hyp(2)/2))/(p_hyp{2}^(p_hyp{2}/2));

z = (gam*kappa*nu/(pi^(n/2))) * prod(p_hyp{3}(3:end));

end