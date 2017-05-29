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