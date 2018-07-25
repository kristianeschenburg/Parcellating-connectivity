import numpy as np
from scipy.special import gammaln
from cachetools import cached, TTLCache 
cache = TTLCache(maxsize=100, ttl=300)

import time

from subgraphs import ConnectedComponents, ClusterSpanningTrees

class ddCRP(object):
    
    """
    Class to implement the distance-dependent Chinese Restaurant Process.
    
    Parameters:
    - - - - -
        alpha : concentration parameter of CRP prior
        mu_0, kappa_0 : hyperparameters on feature mean prior
        nu_0, sigma_0 : hyperparameters on feature variance prior

        mcmc_passes : number of MCMC passes to apply to data
        stats_interval : number of passes to run before recording statistics
        verbose : boolean to print statistics every stats_interval
        
    """
    
    def __init__(self,alpha,mu_0,kappa_0,nu_0,sigma_0,mcmc_passes=100,
        stats_interval=500,verbose=true):
        
        """
        Initialize ddCRP object.
        """
        
        self.alpha = float(alpha)
        self.mu0 = float(mu_0)
        self.kappa0 = float(kappa_0)
        self.nu0 = float(nu_0)
        self.sigma0 = float(sigma_0)
        self.mcmc_passes = int(mcmc_passes)
        self.stats_interval = int(stats_interval)
        self.verbose = verbose

    def fit(self,features,adj_list,init_c=None,gt_z=None,edge_prior=None):

        """
        Main function to fit the distance-dependent Chinese Restaurant Process.Restaurant

        Parameters:
        - - - - -
            features : data array of features for each sample
            adj_list : adjacency list of samples
            init_c : initialized cortical map, default = []
            gt_z : ground truth map for computing normalized mutual information
            edge_prior : nested dictionary, probability of neighboring vertices beloning
                        to same parcels

        """

        stats = {'times': [],
                    'lp': [],
                    'NMI': [],
                    'K': [],
                    'z': np.zeros((self.mcmc_passes,nvoc)),
                    'c': np.zeros((self.mcmc_passes,nvoc))}

        nvox = len(adj_list)

        # initialize parent vector, if not provided
        if not init_c:
            c = np.zeros((nvox,))
            for i in np.arange(nvox):
                neighbors = adj_list[i] + [i]
                c[i] = neighbors[np.random.ra]
                c(i) = neighbors[np.random.randint(low=0,high=len(neighbors))]
        else:
            c = init_c;
        c = c.astype(np.int32)

        # initialize sparse linkage matrix
        G = sparse.csc_matrix((np.ones((nvox,)),(np.arange(nvox),c)), shape=(nvox,nvox))
        G = G.tolil()
        
        # compute initial parcel count and parcel assignments
        [K, z, parcels] = ConnectedComp(G)

        # compute log-likelihood of initial cortical map
        curr_lp = self._fullProbabilityDDCRP(parcels,features)

        max_lp = -1.*np.inf;
        t0 = time.time()
        steps = 0;

        order = np.arange(nvox)

        # perform mcmc_passes of over all samples
        for pas in np.arange(self.mcmc_passes):

            np.random.shuffle(order)

            for i in order:

                if curr_lp > max_lp:
                    max_lp = curr_lp
                    map_z = z

                # remove current link to parent
                G[i,c[i]] = 0

                # if link was self-link
                if c[i] == i:
                    # Removing self-loop, parcellation won't change
                    rem_delta_lp, z_rem, parcels_rem = -np.log(alpha), z, parcels
                else:
                    # otherwise compute new connected components
                    K_rem, z_rem, parcels_rem = ConnectedComp(G)

                    # if number of components changed
                    if K_rem != K:
                        # We split a cluster, compute change in likelihood
                        rem_delta_lp = -LikelihoodDiff(D, parcels_rem, z_rem[i],
                            z_rem[c[i]], hyp, sym)
                    else:
                        rem_delta_lp = 0


    @cached(cache)
    def _fullProbabilityDDCRP(self,parcels,features):
        
        """
        Compute the full log-likelihood of the clustering.
        
        Parameters:
        - - - - -
            parcels : dictionary mapping cluster IDs to data indices
            features : array of features for full dataset
        """
        
        lp = 0
        
        for parc,idx in parcels.items():
            
            sufficient = self._sufficient_statistics(features[idx,:])
            [kappan,nun,sigman] = self._marginal_parameters(sufficient)

            lp += self._LikelihoodCluster(kappan,nun,sigman,sufficient[0])

        return lp
    
    def _LikelihoodCluster(self,kappa,nu,sigma,n):
        
        """
        Compute the log-likelihood of a single cluster.
        
        Parameters:
        - - - - -
            kappa : kappa of cluster, marginalized of mu0,sigma0
            nu : nu of cluster, marginalized of mu0,sigma0
            sigma : sigma of cluster, marginalized of mu0,sigma0
            n : sample size of cluster
        """
        
        p = len(sigma)
        
        # ratio of gamma functions
        gam = gammaln(nu/2) - gammaln(self.nu0/2);
        
        # terms with square roots in likelihood function
        inner = (1/2) * (np.log(self.kappa0) + self.nu0*np.log(self.nu0*self.sigma0) - 
                         np.log(kappa) - nu*np.log(nu) - n*np.log(np.pi));
        
        # sum of sigma_n for each feature
        outer = np.log((1./sigma)).sum()
        
        lp = p*(gam + inner) + outer;
        
        return lp

    def _LogProbDiff(self,parcel_split,split_l1,split_l2,features):

        """
        Compute change in log-likelihood when considering a merge.

        Parameters:
        - - - - -
            parcel_split : indices of each component in split map
            split_l1 , split_l2 : label values of components to merge
            features : input feature array
        """

        merged_indices = np.concatenate([parcel_split[split_l1],
            parcel_split[split_l2]])

        # compute sufficient statistics and marginalized parameters of merged parcels
        stats = self._sufficient_statistics(features[merged_indices,:]);
        phyp = self._marginal_parameters(stats);

        # compute likelihood of merged parcels
        merge_ll = self._LikelihoodCluster(phyp,stats[0]);

        # compute likelihood of split parcels
        split_ll = self._LogProbSplit(parcel_split,split_l1,split_l2,data);

        ld = merge_ll - split_ll;

        return ld

    def _LogProbSplit(self,parcel_split,split_l1,split_l2,features):

        """
        Compute change in log-likelihood when consiering a split.

        Parameters:
        - - - - -
            parcel_split : indices of each component in split map
            split_l1 , split_l2 : label values of components to merge
            features : input feature array
        """

        idx1 = parcel_split[split_l1]
        idx2 = parcel_split[split_l2]

        suff1 = self._sufficient_statistics(features[idx1,:]);
        suff2 = self._sufficient_statistics(features[idx2,:]);

        phyp1 = self._marginal_parameters(suff1);
        phyp2 = self._marginal_parameters(suff2);

        split_ll = self._LikelihoodCluster(phyp1,suff1[0]) + self._LikelihoodCluster(phyp2,suff2[0]);

        return split_ll
 
    def _sufficient_statistics(self,cluster_features):
        
        """
        Compute sufficient statistics for data.
        
        Parameters:
        - - - - -
            cluster_features : data array for single cluster 
        """
        
        # n samples
        [n,_] = cluster_features.shape
        # feature means
        mu = cluster_features.mean(0)
        # feature sum of squares
        ssq = ((cluster_features-mu[None,:])**2).sum(0)

        return [n,mu,ssq]
    
    def _marginal_parameters(self,suff_stats):
        
        """
        Compute cluster-specific marginal hyperparameters after collapsing over mu and sigma.
        
        Parameters:
        - - - - -
            suff_stats : sufficient statistics for single cluster
        """
        
        n = suff_stats[0]
        mu = suff_stats[1]
        ssq = suff_stats[2]

        # update kappa and nu
        kappaN = self.kappa0 + n
        nuN = self.nu0 + n
        
        deviation = ((n*self.kappa0) / (n+self.kappa0)) * (self.mu0 - mu);
        sigmaN = (1/nuN) * (self.nu0*self.sigma0 + ssq + deviation);
        
        return [kappaN,nuN,sigmaN]
    