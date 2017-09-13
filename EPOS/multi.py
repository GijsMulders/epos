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
	print 'cdf(): returns the cdf of planets in multi-systems'
	print '  ID:  array of star identifiers'
	print
	print 'periodratio(): returns the period ratios of adjacent planets'
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

def cdf(ID, Verbose=False):
	bin, count= frequency(ID)
	xlist=[[bin[k]]*count[k]*bin[k] for k in range(len(bin))]
	#if Verbose: print '  multi cdf planets: {}'.format(len(np.concatenate(xlist)))
	return np.concatenate(xlist)

def periodratio(ID, P, N=None, Verbose=False):
	IDsys, toplanet, counts= np.unique(ID, return_inverse=True,return_counts=True)
	ismulti= (counts>1) # systems that are multi
	Pmulti= P[ismulti[toplanet]] # Periods of _all_ planets in multis
	assert np.all(IDsys[toplanet] == ID)
	
	# get index of first planet in each system
	di= np.roll(counts,1)
	di[0]=0
	i1= np.cumsum(di) # index to first planet
	assert np.all(IDsys == ID[i1]), ' assumes lexsort( (P,ID) )'
	
	Pinner= P[i1[counts>1]] # innerplanet in multi
	
	Pratio= []
	for i in range(2,len(np.bincount(counts)) ):
		#print '\nmultiplicity: {}'.format(i)
		im= i1[counts>=i] # multis with 
		dP= P[im+(i-1)]/P[im+(i-2)]
		#for k, _dP in zip(im,Pratio):
		#	print ' ID {}, P={}, dP= {}'.format(ID[k], P[ID==ID[k]], _dP)
		#print 'i={}: {}'.format(i,dP.size)
		Pratio.extend(dP)
	#print 'n dP= {}'.format(len(Pratio))
	
	if N is None:
		return np.array(Pratio), Pinner
	else:
		''' Innermost observed planet is nth planet'''
		PN, dPN= [], []
		for m in range(1,10):
			PN.append(Pinner[N[i1[counts>1]]==m])
			#print m, np.sum(N[i1[counts>1]]==m) # mostly 1 or 2, never 3+
			#print N[i1[counts>1]]
			
			# Period ratio of adjacent planets?
			#print 'm={}'.format(m)
			PratioN= []
			for i in range(2,len(np.bincount(counts)) ):
				im= i1[counts>=i] # multis with 
				idx= N[im+(i-1)]-N[im+(i-2)] == m # adjacent, 1,2,3... inbetween
				#print N[im+(i-1)]- N[im+(i-2)]
				#print i, im.size, 
				#print m, idx.sum()
				dP= P[im[idx]+(i-1)]/P[im[idx]+(i-2)]
				PratioN.extend(dP)
			#print 'm={}, n dPN={}'.format(m, len(PratioN))
			dPN.append(PratioN)

		#return np.array(Pratio), Pinner, PN, dPN
		return PN, dPN
		