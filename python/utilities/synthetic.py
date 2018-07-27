import numpy as np
from scipy import stats

def synthetic_features(z,d,mu_0,kappa_0,nu_0,sigma_0):

	"""
	Sample synthetic features for a given parcellation.def

	Parameters:
	- - - - -
		z : parcellation
		d : feature dimensions
		mu_0, kappa_0 : hyperparameters on prior mean
		nu_0, sigma_0 : hyperparameters on prior variance
	"""

	parcels = {k : np.where(z == k)[0] for k in set(z)}
	params = {k: {'mu': None, 'std': None} for k in set(z)}

	parcel_features = {}.fromkeys(parcels.keys())

	for parc,idx in parcels.items():

		m,s = sample_priors(mu_0,kappa_0,nu_0,sigma_0,d)
		params[parc]['mu'] = m
		params[parc]['std'] = s

		parcel_features[parc] = sample_features(m,s,len(parcels[parc]))

	feature_array = np.zeros((len(z),d))
	for parc,idx in parcels.items():
		feature_array[idx,:] = parcel_features[parc]


	return [parcels,params,parcel_features,feature_array]

def sample_priors(mu_0,kappa_0,nu_0,sigma_0,size):

	"""
	Sample prior mean and variance of synthetic data.

	Parameters:
	- - - - -
		mu_0, kappa_0 : hyperparameters on prior mean
		nu_0, sigma_0 : hyperparameters on prior variance
	"""

	x = stats.chi2.rvs(nu_0,size=[size,])
	sigma_p = (nu_0*sigma_0)/x;
	mu_p = stats.norm.rvs(mu_0,(sigma_p/kappa_0))

	return [mu_p,sigma_p]

def sample_features(mu,sigma,d):

	"""
	Sample from the prior distribution.

	Parameters:
	- - - - -
		mu : prior mean
		sigma : priovar (co)variance
		d : number of samples
	"""

	samples = stats.multivariate_normal.rvs(mean=mu,cov=np.diag(sigma),size=[d,])

	return samples