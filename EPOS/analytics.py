''' Calculate some analytics'''
import numpy as np
import scipy.stats

def pfmodel(epos, sma_cut=1., Verbose=False):
	assert hasattr(epos, 'pfm')
	pfm= epos.pfm
	
	#func= scipy.stats.mstats.gmean
	func=np.median
	
	inc= func(pfm['inc'])
	mass= func(pfm['M'])
	Pin= func(pfm['Pin'])
	dP= func(pfm['dP'][pfm['dP']>0])
	
	_ , counts= np.unique(epos.pfm['ID'],return_counts=True)
	Nk= np.mean(counts)
	
	# How many systems fit in this range with kepler spacing?
	dP_Kep= 2.0 # kepler spacing
	Nobs= np.log10(365.24/Pin) / np.log10(dP_Kep)
	
	s_inc= 'inc= {:.2g} degree'.format(inc)
	s_mass= 'mass= {:.2g} M_earth'.format(mass)
	s_inner= 'inner= {:.2g} days'.format(Pin)
	s_ratio= 'ratio= {:.2g}'.format(dP)
	s_multi= 'Nk= {}'.format(Nk)
	#s_multi= 'Nk= {} ({:.2g})'.format(Nk, Nobs)
	
	s_tab= '${:.2g}$ '.format(mass)
	s_tab+= '& ${:.2g}$ '.format(inc)
	s_tab+= '& ${:.2g}$ '.format(Pin) 
	s_tab+= '& ${:.2g}$ '.format(dP)
	s_tab+= '& ${:.1f}$ '.format(Nk)
	#s_tab+= '& ${:.1f} ~({:.1f})$ '.format(Nk, Nobs)
	s_tab+= r'\\'
	
	if Verbose:
		print s_inc, s_mass, s_inner, s_ratio, s_multi
	
	#return s_inc, s_mass, s_inner, s_ratio, s_multi
	return s_tab

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
