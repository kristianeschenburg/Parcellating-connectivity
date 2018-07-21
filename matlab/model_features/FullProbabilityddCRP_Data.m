function [lp] = FullProbabilityddCRP_Data(hyp,parcels,data)

lp = 0;

for j = 1:length(parcels)
    
    sufficient = suff_stats(data(parcels{j},:));
    p_hyp = marginal_parameters(hyp,sufficient);
    lp = lp + LogProbCluster(hyp,p_hyp,sufficient{1});
    
end

end