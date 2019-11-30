''' Calculate some analytics'''
import numpy as np
import scipy.stats

def bic_loglike(log_likelihood, k_free, n_data):
	# natural log
	return -2.* log_likelihood + k_free* np.log(n_data)

def aic_loglike(log_likelihood, k_free):	
	# natural log
	return 2.* k_free -2.* log_likelihood

def aic_c_loglike(log_likelihood, k_free, n_data):	
	# natural log
	aic= aic_loglike(log_likelihood, k_free)
	cfactor= (2.* k_free**2. + 2.*k_free) / (n_data- k_free - 1.)
	return  aic * cfactor

def bic_rss(rss, k_free, n_data):	
	return n_data*np.log(rss/n_data) + k_free *np.log(n_data)

def aic_rss(rss, k_free, n_data):	
	return 2.* k_free + n_data*np.log(rss) # + some constant?

def aic_c_rss(rss, k_free, n_data):	
	aic= aic_rss(rss, k_free, n_data)
	cfactor= (2.* k_free**2. + 2.*k_free) / (n_data- k_free - 1.)
	return  aic * cfactor
