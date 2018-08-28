'''
his module contains helper functions for multiplanet systems
'''
import numpy as np

def indices(ID, Verbose=False):
	''' 
	retruns indices to planets that are observed to be single or multi
	
	Args:
		ID(np.array):  array of planet host star identifiers. 
			Planets with the same ID are in the same system.
	
	Returns:
		single(np.array of int): index to observed single systems
		multi(np.array of int): index to observed multis systems
		k(list of int): planetary systems witk k planets
		multis(np.array of int): indices to observed systems witk k planets
	
	'''
	IDsys, toplanet, counts= np.unique(ID, return_inverse=True,return_counts=True)
	if Verbose:
		print '\n  {} singles, {} multis'.format(np.sum(counts==1),np.sum(counts>1))
	
	single=(counts==1)[toplanet]
	multi= (counts>1)[toplanet]
	ksys= range(2,len(np.bincount(counts)))
	multis= [(counts==k)[toplanet] for k in ksys] # tuple for indexing
	
	return single, multi, ksys, multis

def nth_planet(ID, P):
	'''
	Index to observed nth planet in the system
	
	Args:
		ID(np.array):  array of planet host star identifiers. 
			Planets with the same ID are in the same system.
		ID(np.array):  array of planet orbital periods. 
	
	Returns:
		single(np.array of int): index to observed single systems
		multi(np.array of int): index to observed multis systems
		nth(list of int): nth planet in the system
		multis(np.array of int): indices to nth planet
	'''
	IDsys, toplanet, counts= np.unique(ID, return_inverse=True,return_counts=True)
	assert np.all(IDsys[toplanet] == ID)

	# get index of first planet in each system
	i1= _first_planet_in_system(counts)	

	single=(counts==1)[toplanet]
	multi= (counts>1)[toplanet]

	ksys= range(1, len(np.bincount(counts))) # [1,2,3,...n]
	#ksys= range(2, len(np.bincount(counts)))
	
	multis=  [i1[counts>1]] # 1st planet in multi
	for i in ksys[1:]:
	#for i in ksys:
		# 2nd, 3rd, ..
		multis.append(tuple([i1[counts>=i]+(i-1)])) # tuple for indexing
		
	return single, multi, ksys, multis

def frequency(ID, Verbose=False):
	'''
	returns the frequency of single/double/triple/etc systems
	
	Args:
		ID(np.array):  array of planet host star identifiers. 
			Planets with the same ID are in the same system.
	'''
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
	'''
	returns the cdf of planets in multi-systems
	'''
	bin, count= frequency(ID)
	xlist=[[bin[k]]*count[k]*bin[k] for k in range(len(bin))]
	#if Verbose: print '  multi cdf planets: {}'.format(len(np.concatenate(xlist)))
	return np.concatenate(xlist)

def periodratio(ID, P, N=None, R=None, Verbose=False):
	'''
	returns the period ratios of adjacent planets as a list
	'''
	IDsys, toplanet, counts= np.unique(ID, return_inverse=True,return_counts=True)

	ismulti= (counts>1) # systems that are multi
	Pmulti= P[ismulti[toplanet]] # Periods of _all_ planets in multis
	assert np.all(IDsys[toplanet] == ID)

	# get index of first planet in each system
	i1= _first_planet_in_system(counts)	
	assert np.all(IDsys == ID[i1]), ' assumes lexsort( (P,ID) )'
	Pinner= P[i1[counts>1]] # innerplanet in multi
	
	Pratio= []
	Rpair= [] # size of inner planet in pair
	for i in range(2,len(np.bincount(counts)) ):
		#print '\nmultiplicity: {}'.format(i)
		im= i1[counts>=i] # multis with 
		_dP= P[im+(i-1)]/P[im+(i-2)]
		if R is not None:
			_R= R[im+(i-1)] # size of outer planet
		#for k, _dP in zip(im,Pratio):
		#	print ' ID {}, P={}, dP= {}'.format(ID[k], P[ID==ID[k]], _dP)
		#print 'i={}: {}'.format(i,dP.size)
		Pratio.extend(_dP)
		if R is not None:
			Rpair.extend(_R)
	#print 'n dP= {}'.format(len(Pratio))
	
	if R is not None:
		return np.array(Pratio), Pinner, Rpair
	elif N is None:
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

def _first_planet_in_system(counts):
	#_, counts= np.unique(ID, return_counts=True) # slower
	di= np.roll(counts,1)
	di[0]=0
	return np.cumsum(di) # index to first planet
