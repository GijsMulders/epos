import h5py
import glob
import numpy as np
import cgs
import sys

#def symba(name='HMSim1', dir='hdf5/Sim1', plts_mass=0, istep=None, Verbose=False):
def symba(fname, plts_mass=0, cut=-np.inf, istep=None, Verbose=False):
	''' 
	returns a list of planetary systems
	sma in au
	mass in earth masses
	remove planetesimals (m < plt_mass)
	cut planets that are still forming (m < cut*au)
	'''
	print '\nProcessing Symba HDF5 file'
	#fname= '{}/{}_set??.h5'.format(dir,name)
	flist= glob.glob(fname)
	if len(flist)==0: raise ValueError('file not found: {}'.format(fname))
	
	sma, mass, inc= [], [], []
	sma0, mass0=[], []
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
				print '  {}: Using planetary system at step {}'.format(fname, nsteps)
			else:
				amtDone= float(i)/len(flist)
				print '\r  [{:50s}] {:5.1f}%'.format('#' * int(amtDone * 50), amtDone * 100),
				sys.stdout.flush() 
			
			L_sma, L_mass, L_inc= [], [], []
			L_sma0, L_mass0= [], []
			
			for particle in hf:
				#if Verbose: print particle, hf.get(particle).shape[0]
				if hf.get(particle).shape[0] == nsteps:
					final_architecture= np.array(hf.get(particle))[-1,:]
					#print final_architecture.shape

 					_mass= final_architecture[8]* symbamassunit
 					_sma= final_architecture[2]
 					if (_mass > plts_mass) & (_mass > cut * _sma**1.5):
						L_sma.append(_sma)
						L_mass.append(_mass)
						L_inc.append(final_architecture[4])
					#if not 'time' in sg: sg['time']= final_architecture[0]/1e6
				
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
			
			# story copy(?) in list of planetary systems
			sma.append(L_sma)
			mass.append(L_mass)
			inc.append(L_inc)
			#sma0.append(L_sma0)
			#mass0.append(L_mass0)
		if not Verbose: print '\r  [{:50s}] {:.1f}%'.format('#' * int(1 * 50), 1 * 100),
	
# 	print sma0
# 	sma_min= min([min(x) for x in sma0])
# 	sma_max= max([max(x) for x in sma0])
# 	print '{}: {:.2f}-{:.1f} au'.format(name, sma_min,sma_max)
#	print '    {:.3f}-{:.3f} Mearth'.format(name, np.min(mass0),np.max(sma))
	 
			
	return sma, mass, inc
	#return {'sma':sma, 'mass':mass, 'inc':inc}
	

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
	
def dace_screengrab(name='CD753', dir='DACE', Verbose=False):
	fname= '{}/{}.txt'.format(dir,name)
	#with open(fanem,'r') as f:
	header= np.genfromtxt(fname, max_rows=1, dtype=str)
	print header
	a= np.loadtxt(fname, unpack=True, skiprows=1, usecols=(2,3,4,5), delimiter='\t')
	sma= a[1]
	mass= a[2]
	radius= a[3]

	print '\nLoad DACE model {}'.format(name)
	print '  sma:  	 {:.2e} ... {:.1f}'.format(min(sma), max(sma))
	print '  mass:   {:.2e} ... {:.1f}'.format(min(mass), max(mass))
	print '  radius: {:.2f} ... {:.1f}'.format(min(radius), max(radius))
	
	return sma, mass, radius