import matplotlib.pyplot as plt
import numpy as np
import helpers

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def readme():
	print '\nMakes all the plots'
	print 
	print 'observed:'
	print 'model:'
	print 'result:'


''' parameter study '''
def multiD(epos):
	print '\nPlotting parameter study...'
	
	# identify indices to dimensions with arrays
	k_dim= []
	for k, dim in enumerate(epos.para['grid']): 
		if len(dim)>1: k_dim.append(k) 
	
	all_dim= epos.para['prob'].shape
	ndim= len(k_dim) 
	indices= range(len(all_dim)) 
				
	print '  {} dimensional; {}'.format(ndim, k_dim)
	if ndim>1:
		for dim in k_dim:
			# sum over all axes _except_ dim :S
			index= tuple(indices[x] for x in indices if (x is not dim))
			print dim, index
			p_1d= np.sum(epos.para['prob'], axis=index)
			print p_1d.shape
			para_1d(epos, p_1d, epos.para['grid'][dim], 'dim{}'.format(dim))

def oneD(epos, pdf, xgrid, fname):
	f, ax = plt.subplots()
	
	ax.set_title('Marginalized PDF')
	ax.set_xlabel(fname)
	ax.set_ylabel('Probability')

	ax.plot(xgrid, pdf, ls='-', marker='+', color='b')
	
	aux.save(plt, '{}grid/1d.{}'.format(epos.plotdir,fname))		