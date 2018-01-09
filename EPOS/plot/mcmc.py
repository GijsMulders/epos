import warnings
import numpy as np
import matplotlib.pyplot as plt

import helpers
import parametric, periodradius, multi

try:
	import corner
except ImportError:
	print '\nWarning: corner.py not imported'
	warnings.warn('corner.py not imported',ImportWarning)
except:
	print '??'

def all(epos):
	if hasattr(epos, 'chain'):
		print '\nPlotting chain...'
		chain(epos)
		try:
			corners(epos)
		except NameError:
			print '  (skipping corner plot)'
		
		if epos.populationtype is 'parametric':
			parametric.oneD(epos, MCMC=True)
			parametric.twoD(epos, MCMC=True)
		
		# plot sample from posterior
		if hasattr(epos, 'plotsample'):
			periodradius.panels(epos, MCMC=True)
			if epos.Multi:
				multi.multiplicity(epos, MCMC=True)
				multi.periodratio(epos, MCMC=True)
				multi.periodinner(epos, MCMC=True)
			
	else:
		print '\nNo chain to plot, did you run EPOS.run.mcmc()? \n'
	
def chain(epos):
	nwalker, nstep, npara= epos.chain.shape

	# plot multiplicity CDF
	f, axlist = plt.subplots((npara+1)/2, 2, sharex=True)
	f.set_size_inches(8,8)

	#axlis.set_title('')
	for axl in axlist[-1]:
		axl.set_xlabel('Step number')
	
	# loop over existing axes only
	for k, (xlabel, ax) in enumerate(zip(epos.fitpars.keysfit,axlist.flatten())):
		ax.set_ylabel(xlabel)
	
		#ax.set_xlim(1, 10)
	
		#plot each walker
		for i in range(nwalker):
			ax.plot(epos.chain[i,:,k],color='r',alpha=10./nwalker)
	
		ax.axvline(epos.burnin, ls='--', color='k')
	
	f.subplots_adjust(hspace=0.15, wspace=0.4)
			
	helpers.save(plt, '{}mcmc/chain'.format(epos.plotdir))	

def corners(epos):
	fig = corner.corner(epos.samples, labels=epos.fitpars.keysfit,
                      truths=epos.fitpars.getfit(Init=True), 
                      quantiles=[0.16, 0.5, 0.84], show_titles=True)
	fig.savefig('{}mcmc/triangle.png'.format(epos.plotdir))
		
