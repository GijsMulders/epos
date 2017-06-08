import numpy as np
import matplotlib.pyplot as plt
import helpers
import corner

def chain(epos):
	nwalker, nstep, npara= epos.chain.shape
	print nwalker, nstep, npara

	# plot multiplicity CDF
	f, axlist = plt.subplots((npara+1)/2, 2, sharex=True)
	f.set_size_inches(8,6)

	#axlis.set_title('')
	for axl in axlist:
		axl[-1].set_xlabel('Step number')
	
	# loop over existing axes only
	for k, (xlabel, ax) in enumerate(zip(epos.spara,axlist.flatten())):
		ax.set_ylabel(xlabel)
	
		#ax.set_xlim(1, 10)
	
		#plot each walker
		for i in range(nwalker):
			ax.plot(epos.chain[i,:,k],color='k')	
	
	f.subplots_adjust(hspace=0)
			
	helpers.save(plt, '{}mcmc/chain'.format(epos.plotdir))	

def corners(epos, nburn=0):
	chain1d= epos.chain[:, nburn:, :].reshape((-1, len(epos.spara)))
	fig = corner.corner(chain1d, labels=epos.spara,
                      truths=epos.para_in, 
                      quantiles=[0.16, 0.5, 0.84], show_titles=True)
	fig.savefig('{}mcmc/triangle.png'.format(epos.plotdir))
	
