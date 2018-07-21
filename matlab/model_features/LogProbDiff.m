% Compute change in log-probability when comsidering a new merge
% hyp : hyperparameters
% parcels_split : indices of each parcel in the split map
% split_l1 : label of split region 1
% split_l2 : label of split region 2
% data : input data matrix

function [ld] = LogProbDiff(hyp,parcel_split,split_l1,split_l2,data)

	merge_indices = [parcel_split{split_l1}, parcel_split{split_l2}];
	merge_suff = suff_stats(data(merge_indices,:));
	merge_phyp = marginal_parameters(hyp,merge_suff);

	merge_ll = LogProbCluster(hyp,merge_phyp,merge_suff{1});
	split_ll = LogProbSplit(hyp,parcel_split,split_l1,split_l2,data);

	ld = merge_ll - split_ll;

end

% Compute the log-probability of splitting a parcel into two
% hyp : hyperparameters
% parcels_split : indices of each parcel in the split map
% split_l1 : label of split region 1
% split_l2 : label of split region 2
% data : input data matrix

function [split_ll] = LogProbSplit(hyp,parcel_split,split_l1,split_l2,data)

	suff1 = suff_stats(data(parcel_split{split_l1},:));
	suff2 = suff_stats(data(parcel_split{split_l2},:));

	phyp1 = marginal_parameters(hyp,suff1);
	phyp2 = marginal_parameters(hyp,suff2);

	split_ll = LogProbCluster(hyp,phyp1,suff1{1}) + LogProbCluster(hyp,phyp2,suff2{1});

end