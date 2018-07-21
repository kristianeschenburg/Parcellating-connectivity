function [stats] = suff_stats(data)

stats = cell(3,1);

% mean of each feature
mu = mean(data,1);
% sum of squares of each feature
ssq = sum(((data - repmat(mu,size(data,1),1)).^2),1);

stats{1} = size(data,1);
stats{2} = mu;
stats{3} = ssq;

end