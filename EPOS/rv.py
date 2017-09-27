import numpy as np

'''
NOTE: RV from Howard 2010
'''

def readme():
	#import clearscreen # works only on loading, not on execution
	print '\nThis module contains helper functions to load RV survey data into EPOS:'
	print
	print 'Howard2010(): returns the Howard et al. 2010 exoplanet catalogue as numpy arrays'
	print '  P:   Orbital period in days'
	print '  Rp:  Planet radius in earth radii'
	print '  KID: star identifier (for multis, not implemented)'
	print
	print 'eff_Q16(): returns the Q16 survey detection efficiency as a 2D matrix'
	print '  P:   Orbital period in days'
	print '  Rp:  Planet radius in earth radii'
	print '  eff: 2D matrix of detection efficiency [0-1]'
	print
	print 'dr25(): returns the DR25 planet list and '
	print '        survey detection efficiency as a 2D matrix'
	print '  P:   Orbital period in days'
	print '  Rp:  Planet radius in earth radii'
	print '  eff: 2D matrix of detection efficiency [0-1]'


def Howard2010():
	Msini,P,completeness = \
		np.loadtxt('files/msini_per_completeness.txt',usecols=(1,0,2),unpack=True)
	P_obs, Msini_obs = \
		np.loadtxt('files/planet_star_per_msini_ref.txt',usecols=(2,3),unpack=True)
	ID = np.loadtxt('files/planet_star_per_msini_ref.txt',usecols=(1,),dtype=str)

	''' make 2D grid '''
	P_1D= np.unique(P)
	Msini_1D= np.unique(Msini)

	X, Y = np.meshgrid(P_1D, Msini_1D)
	Z= completeness.reshape(X.shape)
	
	obs= {'xvar':P_obs,
		'yvar':Msini_obs, 
		'starID':ID}
	
	survey= {'xvar':P_1D, 'yvar':Msini_1D, 'eff_2D':Z.T} # 'Mstar':None, 'Rstar':None}
	obs['nstars']= 125 # 
	
	return obs, survey