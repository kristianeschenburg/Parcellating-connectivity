%%%
% Computes the posterior parameters for a single cluster data.

% hyp : vector of prior parameters [kappa_0, nu_0, sigma_0]
% stats : cell array of sufficient statistics for data
%%%

function [p_hyp] = marginal_parameters(hyp,stats)

p_hyp = cell(3,1);

n = stats{1};
mu = stats{2};
ssq = stats{3};

% update kappa and nu with data
kappa_n = hyp(1) + n;
nu_n = hyp(2) + n;

% update sigma with data
% original update formula is n*kappa_0 / (n+kappa_0) * (mu_0 - x_bar)
% here we assume that mu_0 is 0
% deviation = ((n*hyp(2))/(n+hyp(2))) * (hyp(1)-mu);

deviation = (n*hyp(1) / (n+hyp(1))) * (-mu).^2;
sigma_n = (1/nu_n) * (hyp(2)*hyp(3) + ssq + deviation);

% return parameters updated with data
p_hyp{1} = kappa_n;
p_hyp{2} = nu_n;
p_hyp{3} = sigma_n;

end