import numpy as np

def readme():
	print '\nThis module contains helper functions for multiplanet systems:'
	print
	print 'indices(): returns the indices to single/multi planet systems'
	print '  ID:  array of star identifiers'
	print
	print 'frequency(): returns the frequency of single/double/triple/etc systems'
	print '  ID:  array of star identifiers'
	print
	print 'frequency(): returns the period ratios of adjecent planets'
	print '  ID:  array of star identifiers'
	print '  P:  array of orbital periods'
	print

def indices(ID, Verbose=False):
	IDsys, toplanet, counts= np.unique(ID, return_inverse=True,return_counts=True)
	if Verbose:
		print '  {} single systems, {} multis'.format(np.sum(counts==1),np.sum(counts>1))

	return counts==1, counts>1

def frequency(ID, Verbose=False):
	smulti= ['single','double','triple','quad','quint','sext','sept','oct','nint']
	#bc= np.bincount(np.bincount(ID)) # only for int
	_ , counts= np.unique(ID,return_counts=True)
	bincounts= np.bincount(counts)
	assert bincounts[0]==0, "unique items can't have frequency zero"
	if Verbose:
		for nmulti, text in zip(bincounts[1:], smulti):
			print '  - {}: {}'.format(text, nmulti)
	
	return np.arange(1,bincounts.size), bincounts[1:]

def periodratio(ID, P, Verbose=False):
	IDsys, toplanet, counts= np.unique(ID, return_inverse=True,return_counts=True)
	ismulti= (counts>1) # systems that are multi
	Pmulti= P[ismulti[toplanet]] # Periods of _all_ planets in multis
	assert np.all(IDsys[toplanet] == ID)
	
	# get index of first planet in each system
	di= np.roll(counts,1)
	di[0]=0
	i1= np.cumsum(di) # index to first planet
	assert np.all(IDsys == ID[i1]), ' assumes lexsort( (P,ID) )'
	
	Pratio= []
	for i in range(2,len(np.bincount(counts)) ):
		#print '\nmultiplicity: {}'.format(i)
		im= i1[counts>=i] # multis with 
		dP= P[im+(i-1)]/P[im+(i-2)]
		#for k, _dP in zip(im,Pratio):
		#	print ' ID {}, P={}, dP= {}'.format(ID[k], P[ID==ID[k]], _dP)
		Pratio.extend(dP)
	
	return Pratio