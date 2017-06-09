import numpy as np
import matplotlib.pyplot as plt
import helpers
import corner
import parametric

def all(epos):
	chain(epos)
	corners(epos)
	parametric.oneD(epos, MCMC=True)
	parametric.twoD(epos, MCMC=True)
	

def chain(epos):
	nwalker, nstep, npara= epos.chain.shape

	# plot multiplicity CDF
	f, axlist = plt.subplots((npara+1)/2, 2, sharex=True)
	f.set_size_inches(8,8)

	#axlis.set_title('')
	for axl in axlist[-1]:
		axl.set_xlabel('Step number')
	
	# loop over existing axes only
	for k, (xlabel, ax) in enumerate(zip(epos.pname,axlist.flatten())):
		ax.set_ylabel(xlabel)
	
		#ax.set_xlim(1, 10)
	
		#plot each walker
		for i in range(nwalker):
			ax.plot(epos.chain[i,:,k],color='r',alpha=10./nwalker)
	
		ax.axvline(epos.burnin, ls='--', color='k')
	
	f.subplots_adjust(hspace=0.15, wspace=0.4)
			
	helpers.save(plt, '{}mcmc/chain'.format(epos.plotdir))	

def corners(epos):
	fig = corner.corner(epos.samples, labels=epos.pname,
                      truths=epos.p0, 
                      quantiles=[0.16, 0.5, 0.84], show_titles=True)
	fig.savefig('{}mcmc/triangle.png'.format(epos.plotdir))
		