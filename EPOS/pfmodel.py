import h5py
import glob
import numpy as np
import sys, os

import cgs

''' Helper functions to read in planet formation models'''

#def symba(name='HMSim1', dir='hdf5/Sim1', plts_mass=0, istep=None, Verbose=False):
def symba(name, fname, plts_mass=0, cut=-np.inf, smacut=np.inf, istep=None, 
	Verbose=False, Saved=True):
	''' 
	returns a list of planetary systems
	sma in au
	mass in earth masses
	remove planetesimals (m < plt_mass)
	cut planets that are still forming (m < cut*au)
	'''
	
	dir= 'npz/{}'.format(name)
	if not os.path.exists(dir): os.makedirs(dir)
	fnpz= '{}/{}.npz'.format(dir, name)	
	
	''' Load hdf5 file or npz dictionary for quicker access'''
	if os.path.isfile(fnpz) and Saved:
		print '\nLoading saved status from {}'.format(fnpz)
		npz= np.load(fnpz)
				
		# check if keys present
		for key in ['sma','mass','inc','starID']:
			if not key in npz: 
				raise ValueError('Key {} not present\n{}'.format(key,npz.keys()))
	else:
		print '\nProcessing Symba HDF5 file for {}'.format(name)
		#fname= '{}/{}_set??.h5'.format(dir,name)
		flist= glob.glob(fname)
		if len(flist)==0: 
			raise ValueError('file not found: {}'.format(fname))
		else:
			if Verbose: print '  {} files'.format(len(flist))
	
		sma, mass, inc, ecc, ID= [], [], [], [], []
		sma0, mass0, inc0, ecc0, ID0 =[], [], [], [], []
		symbamassunit= cgs.Msun/cgs.Mearth / (2.*np.pi)**2.
	
		for i,fname in enumerate(flist):
			with h5py.File(fname,'r') as hf:

				# get final steps
				if istep is None:
					nsteps=0
					for particle in hf:
						nsteps=max(nsteps, hf.get(particle).shape[0])
				else:
					nsteps=istep

				if Verbose: 
					print '\n{}, step {}'.format(fname, nsteps)
				else:
					amtDone= float(i)/len(flist)
					print '\r  [{:50s}] {:5.1f}%'.format('#' * int(amtDone * 50), amtDone * 100),
					sys.stdout.flush() 
			
				L_sma, L_mass, L_inc, L_ecc= [], [], [], []
				L_sma0, L_mass0, L_inc0, L_ecc0= [], [], [], []
			
				for particle in hf:
					#if Verbose: print particle, hf.get(particle).shape[0]
					if hf.get(particle).shape[0] == nsteps:
						final_architecture= np.array(hf.get(particle))[-1,:]

						_mass= final_architecture[8]* symbamassunit
						_sma= final_architecture[2]
						if (_mass>plts_mass) & (_mass>cut*_sma**1.5) & (_sma<smacut):
							L_sma.append(_sma)
							L_mass.append(_mass)
							L_inc.append(final_architecture[4])
							L_ecc.append(final_architecture[3])
							ID.append(i)

					# initial conditions
					if hf.get(particle).shape[0]== 1:
						pass
					else:
						#print particle
						initial_architecture= np.array(hf.get(particle))[0,:]
						L_sma0.append(initial_architecture[2])
						L_mass0.append(initial_architecture[8]* symbamassunit)
						L_inc0.append(initial_architecture[4])
						L_ecc0.append(initial_architecture[3])
						ID0.append(i)

						#if not 'time' in sg: sg['time']= final_architecture[0]/1e6
				
					# 			# print '\nLoaded subgroup {} with {} planetary systems'.format(sg['name'], sg['n'])
	# 			for system in sg['system']:
	# 				print 'system has {} planets, {:.1f} Mearth:'.format(
	# 					system['np'],np.sum(system['mass']))
	# 				for a,m in zip(system['sma'], system['mass']):
	# 					print '  {:.2f} au, {:.1f} Mearth'.format(a,m)
					#print sg['all_Pratio']

				
					#initial_condition= np.array(hf.get(particle))[0,:]
					#print hf.get(particle).shape
					#L_sma0.append(initial_condition[2]) ?? crashes at HMSim3
					#L_mass0.append(initial_condition[8]* symbamassunit)
				
					# -Col 1: time in yrs
					# -Col 2: Body identifier
					# -Col 3: Semi major axis
					# -Col 4: Eccentricity
					# -Col 5: Inclination
					# -Col 6: Longitude of ascending node
					# -Col 7: Argument of Periapsis
					# -Col 8: Mean anomaly at epoch
					# -Col 9: Mass in units where 1 solar mass = (2*pi)^2
			
				# print system properties:
				if Verbose:
					print 'System {} has {} planets, {:.1f} Mearth:'.format(
						i, len(L_sma),np.sum(L_mass))
					for a,m in zip(L_sma, L_mass):
						print '  {:.2f} au, {:.1f} Mearth'.format(a,m)
			
				sma.extend(L_sma)
				mass.extend(L_mass)
				inc.extend(L_inc)
				ecc.extend(L_ecc)
				#ID (set in loop)
			
				sma0.extend(L_sma0)
				mass0.extend(L_mass0)
				inc0.extend(L_inc0)
				ecc0.extend(L_ecc0)

			if not Verbose: print '\r  [{:50s}] {:.1f}%'.format('#' * int(1 * 50), 1 * 100),
		
		npz={'sma':np.asarray(sma), 'mass':np.asarray(mass), 
			'inc':np.asarray(inc), 'starID':np.asarray(ID),
			'ecc':np.asarray(ecc),
			'sma0':np.asarray(sma0), 'mass0':np.asarray(mass0), 
			'inc0':np.asarray(inc0), 'starID0':np.asarray(ID0),
			'ecc0':np.asarray(ecc0)}
		
		#if Saved:
		print '\nSaving status in {}'.format(fnpz)
		#np.save(fname, epos.chain)
		# compression slow on loading?
		np.savez_compressed(fnpz, **npz)
		
	return npz
	
def mercury(fname, istep=None, icut=-np.inf, smacut=[0,np.inf], masscut= [0,np.inf],
	Verbose=False):
	print '\nProcessing Mercury file'
	flist= glob.glob(fname)
	if len(flist)==0: raise ValueError('file not found: {}'.format(fname))
	sma, mass, inc, ID= [], [], [], []
	
	for i,fname in enumerate(flist):

		if Verbose: 
			print '  {}: Using planetary system at final step'.format(fname)
		else:
			amtDone= float(i)/len(flist)
			print '\r  [{:50s}] {:5.1f}%'.format('#' * int(amtDone * 50), amtDone * 100),
			sys.stdout.flush() 
		
		try:
			a= np.loadtxt(fname, unpack=True,ndmin=2) # ndmin to guarantee lists

			L_sma= a[1]
			L_mass= a[7] * cgs.Msun/cgs.Mearth
			L_inc= a[3]
		
			# story copy(?) in list of planetary systems
			# cut out systems with low inclinatons (optional)
			if np.median(L_inc) > icut:
				sma.extend(L_sma)
				mass.extend(L_mass)
				inc.extend(L_inc)
				ID.extend([i]*len(L_sma))

		except ValueError:
			print '\n(skipping {})\n'.format(fname)

		except: raise

		if not Verbose: print '\r  [{:50s}] {:.1f}%'.format('#' * int(1 * 50), 1 * 100),

	sma= np.array(sma)
	mass=np.array(mass)
	inc=np.array(inc)
	starID=np.array(ID)

	# cut out certain planets
	include= (smacut[0]<=sma) & (sma<=smacut[-1])
	include&= (masscut[0]<=mass) & (mass<=masscut[-1])

	npz={'sma':sma[include], 'mass':mass[include], 
		'inc':inc[include], 'starID':starID[include]}
		
	return npz

def pa_bert(name='1Dlin', dir='PA_Bert/', Verbose=False):
	fname= '{}/Data{}.out'.format(dir,name)
	#with open(fanem,'r') as f:
	a= np.loadtxt(fname, unpack=True)
	include= np.isfinite(a[4]) & (a[4]>0) # & (a[-1]<-0.2) # metallicity cut
	sma= a[4][include]
	mass= a[5][include]
	FeH= a[-1][include]
	
	print '\nLoad Pebble Accretion from Bertram Bitsch'
	print '  sma:  {:.2e} ... {:.1f}'.format(min(sma), max(sma))
	print '  mass: {:.2e} ... {:.1f}'.format(min(mass), max(mass))
	print '  Fe/H: {:.2f} ... {:.2f}'.format(min(FeH), max(FeH))
	
	return sma, mass, FeH
	
def bern(name='syntheticpop_20emb_983systems.txt', dir='Bern', 
		smacut=[0,np.inf], masscut= [0,np.inf], Rcut=[0,np.inf], 
		Verbose=False, Single=False):
	fname= '{}/{}'.format(dir,name)

	''' Read the header '''
	header= np.genfromtxt(fname, max_rows=1, dtype=str, delimiter=',')
	if Verbose:
		#print header
		print '\nColumns:'
		for k, colname in enumerate(header):
			print '  a[{}]: {}'.format(k, colname)
	#print k

	''' Read the data '''
	a= np.loadtxt(fname, unpack=True, skiprows=1)
	if Verbose:
		print '\nraw data: {} systems, {} planets'.format(np.unique(
			a[0 if Single else 1]).size, a[0].size)

	# cut out certain planets
	include= (smacut[0]<=a[2]) & (a[2]<=smacut[-1])
	include&= (masscut[0]<=a[3]) & (a[3]<=masscut[-1])
	include&= (Rcut[0]<=a[4]) & (a[4]<=Rcut[-1])
	#include&= (cut[0]<=a[]) & (a[]<=cut[-1])

	#planet= a[0][include]
	ID= a[0 if Single else 1][include]
	sma= a[2][include]
	mass= a[3][include]
	radius= a[4][include]
	ecc= a[5][include]
	inc=a[6][include] * (180./np.pi) # rad -> degrees
	FeH= a[7][include]
	Mcore= a[8][include]
	Menv= a[9][include]
	sma0= a[10][include]
	fice= a[11][include]

	print '\nLoad population synthesis model {}'.format(name)
	print '  included {} systems, {} planets'.format(np.unique(ID).size, sma.size)
	print '  sma:  	 {:.2e} ... {:.1f}'.format(min(sma), max(sma))
	print '  sma0: 	 {:.2e} ... {:.1f}'.format(min(sma0), max(sma0))
	print '  mass:   {:.2e} ... {:.1f}'.format(min(mass), max(mass))
	print '  radius: {:.2f} ... {:.1f}'.format(min(radius), max(radius))
	print '  inc:    {:.2e} ... {:.1f}'.format(min(inc), max(inc))
	print '  Fe/H:   {:.2f} ... {:.2f}'.format(min(FeH), max(FeH))
	
	order= np.lexsort((sma,ID)) 
	
	npz={'sma':sma[order], 'mass':mass[order], 'radius':radius[order], 
		'inc':inc[order], 'ecc':ecc[order], 'starID':ID[order], 
		'Mcore':Mcore[order], 'Menv':Menv[order], 'fice':fice[order],
		'tag':FeH[order], 'sma0':sma0[order]}
	
	if Single: npz['inc']=None

	return npz

def morby(name, fname, Verbose=False, Saved=True):
	''' 
	returns a list of planetary systems
	sma in au
	mass in earth masses
	'''
	
	dir= 'npz/{}'.format(name)
	if not os.path.exists(dir): os.makedirs(dir)
	fnpz= '{}/{}.npz'.format(dir, name)	
	
	''' Load hdf5 file or npz dictionary for quicker access'''
	if os.path.isfile(fnpz) and Saved:
		print '\nLoading saved status from {}'.format(fnpz)
		npz= np.load(fnpz)
				
		# check if keys present
		for key in ['sma','mass','inc','starID']:
			if not key in npz: 
				raise ValueError('Key {} not present\n{}'.format(key,npz.keys()))
	else:
		print '\nProcessing Symba HDF5 file for {}'.format(name)
		#fname= '{}/{}_set??.h5'.format(dir,name)
		flist= glob.glob(fname)
		if len(flist)==0: 
			raise ValueError('file not found: {}'.format(fname))
		else:
			if Verbose: print '  {} files'.format(len(flist))
	
		sma, mass, inc, ecc, ID= [], [], [], [], []
	
		for i,fname in enumerate(flist):
			a= np.loadtxt(fname,skiprows=1, usecols=(1,2,3,7), unpack=True)
			#if a.ndim==1: a=[a]
			try:
				if a.ndim==1:
					sma.append([a[0]])
					ecc.append([a[1]])
					inc.append([a[2]])
					mass.append([a[3]])

					ID.append([i])
				elif a.ndim==2:

					sma.append(a[0])
					ecc.append(a[1])
					inc.append(a[2])
					mass.append(a[3])

					ID.append([i]*len(a[0]))
				else:
					print a.shape
					print a
					raise ValueError('Error reading in') 
			except (IndexError, TypeError):
				print a.shape
				print a
				raise

			# print system properties:
			if Verbose:
				print 'System {} has {} planets, {:.1f} Mearth:'.format(
					i, len(sma[-1]),np.sum(mass[-1]))
				for items in zip(sma[-1], mass[-1], inc[-1]):
					print '  {:.2f} au, {:.1f} Mearth, {:.3g} deg'.format(*items)
		
		npz={'sma':np.concatenate(sma), 'mass':np.concatenate(mass), 
			'inc':np.concatenate(inc), 'starID':np.concatenate(ID),
			'ecc':np.concatenate(ecc)}
		
		if Saved:
			print 'Saving status in {}'.format(fnpz)
			#np.save(fname, epos.chain)
			# compression slow on loading?
			np.savez_compressed(fnpz, **npz)
		
	return npz


def mordasini(name='syntheticpopmordasini1MsunJ31', dir='Mordasini', smacut=np.inf,
		Rcut=0, Single=False, Verbose=False):
	fname= '{}/{}.dat'.format(dir,name)
	header= np.genfromtxt(fname, max_rows=1, dtype=str)
	print header
	a= np.loadtxt(fname, unpack=True, skiprows=1, usecols=(0,1,2,3,4,6,7))

	npl= a[2].size
	ns= np.unique(a[0 if Single else 1]).size

	include= (a[2]<smacut) & (a[4]>Rcut)
	ID= a[0 if Single else 1][include]
	sma= a[2][include]
	mass= a[3][include]
	radius= a[4][include]
	inc=a[5][include]
	FeH= a[6][include]

	print '\nLoad population synthesis model {}'.format(name)
	print '  {} stars, {} planets'.format(ns, npl)
	print '  sma:  	 {:.2e} ... {:.1f}'.format(min(sma), max(sma))
	print '  mass:   {:.2e} ... {:.1f}'.format(min(mass), max(mass))
	print '  radius: {:.2f} ... {:.1f}'.format(min(radius), max(radius))
	print '  inc:    {:.2e} ... {:.1f}'.format(min(inc), max(inc))
	print '  Fe/H:   {:.2f} ... {:.2f}'.format(min(FeH), max(FeH))
	
	order= np.lexsort((sma,ID)) 
	
	npz={'sma':sma[order], 'mass':mass[order], 'radius':radius[order], 
		'inc':inc[order], 'starID':ID[order], 'tag':FeH[order]}

	if Single: npz['inc']=None
		
	return npz

def mordasini_ext(name='syntheticpopmordasini1MsunJ31extended', dir='Mordasini', smacut=np.inf,
		Rcut=0, Verbose=False):
	fname= '{}/{}.dat'.format(dir,name)
	header= np.genfromtxt(fname, max_rows=1, dtype=str)
	print header
	a= np.loadtxt(fname, unpack=True, skiprows=1, usecols=(1,2,3,4,6,7,10))

	npl= a[2].size
	ns= np.unique(a[0 if Single else 1])

	include= (a[1]<smacut) & (a[3]>Rcut)
	ID= a[0][include]
	sma= a[1][include]
	mass= a[2][include]
	radius= a[3][include]
	inc=a[4][include]
	FeH= a[5][include]
	sma0= a[6][include]

	print '\nLoad population synthesis model {}'.format(name)
	print '  {} stars, {} planets'.format(ns, npl)
	print '  sma:  	 {:.2e} ... {:.1f}'.format(min(sma), max(sma))
	print '  sma0: 	 {:.2e} ... {:.1f}'.format(min(sma0), max(sma0))
	print '  mass:   {:.2e} ... {:.1f}'.format(min(mass), max(mass))
	print '  radius: {:.2f} ... {:.1f}'.format(min(radius), max(radius))
	print '  inc:    {:.2e} ... {:.1f}'.format(min(inc), max(inc))
	print '  Fe/H:   {:.2f} ... {:.2f}'.format(min(FeH), max(FeH))
	
	order= np.lexsort((sma,ID)) 
	
	npz={'sma':sma[order], 'mass':mass[order], 'radius':radius[order], 
		'inc':inc[order], 'starID':ID[order], 'tag':FeH[order], 'sma0':sma0[order]}
		
	return npz

def combine(pfmlist, tags):
	''' Combine Multiple Planet Formation Model'''
	assert len(tags) == len(pfmlist)

	''' Create list of tags '''
	taglist= []
	for pfm, tag in zip(pfmlist, tags):
		taglist.append(np.full_like(pfm['sma'], tag))

	''' Create merged dictionary'''
	npz={'tag':np.concatenate(taglist)}
	for key in pfm.keys():
		npz[key]= np.concatenate([pfm[key] for pfm in pfmlist])

	return npz