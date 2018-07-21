% Compute full probability of a given parcellation of D, specified both in terms
%   of voxel links c and list of arrays of element indices "parcels".
%   Hyperparmeters as specified as alpha and vectorized hyp, and whether D is
%   symmetric is given by the boolean sym.
%   Note that this is very slow for large matrices, and should only be used
%   during initialization - the likelihood is updated incrementally during inference
function logp = FullProbabilityddCRP(D, c, parcels, alpha, hyp, sym)

    % if the input similarity matrix is symmetric
    % only need to compute sufficient statistics for half of matrix
    if (sym)
        stats = zeros(length(parcels)*(length(parcels)+1)/2,3);
        j = 1;
        for c1 = 1:length(parcels)
            for c2 = c1:length(parcels)
                samples = D(parcels{c1},parcels{c2});
                if (c1 == c2)
                    samples = samples(logical(triu(ones(size(samples)),1)));
                    if (isempty(samples))
                        continue;
                    end
                else
                    samples = samples(:);
                end
                stats(j,1) = length(samples);
                stats(j,2) = sum(samples)/stats(j,1);
                stats(j,3) = sum((samples-stats(j,2)).^2);
                j = j+1;
            end
        end
        logp = log(alpha) * sum(c' == 1:length(c)) + LogLikelihood(stats, hyp);
        
    % if the input similarity matrix is asymmetric
    else
        
        % initialize pairwise parcel statistics
        stats = zeros(length(parcels)*length(parcels),3);
        j = 1;
        
        % loop over each pair of parcels
        for c1 = 1:length(parcels)
            for c2 = 1:length(parcels)
                
                % get data matrix corresonponding to pair of parcels
                samples = D(parcels{c1},parcels{c2});
                
                % if its the diagonal elements
                if (c1 == c2)
                    
                    % create matrix of all ones
                    off_diags = true(size(samples));
                    % set diagonal elements to 0 / false
                    off_diags(1:(size(samples,1)+1):end) = false;
                    % get off diagonal elements
                    samples = samples(off_diags);
                    if (isempty(samples))
                        continue;
                    end
                else
                    samples = samples(:);
                end
                
                %%% Compute suff. stats. of pairwise similarities %%%
              
                % get number of samples
                stats(j,1) = length(samples);
                % get mean of samples
                stats(j,2) = sum(samples)/stats(j,1);
                % compute sum of squared differences
                stats(j,3) = sum((samples-stats(j,2)).^2);
                j = j+1;
            end
        end
        logp = log(alpha) * sum(c' == 1:length(c)) + LogLikelihood(stats, hyp);
    end
end
