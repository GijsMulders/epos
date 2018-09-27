import numpy as np
import os.path
import EPOS

fpath= os.path.dirname(EPOS.__file__)

'''
This module contains helper functions to load RV survey data into EPOS:'
'''

def Mayor2011():
	Msini, P, Completeness = np.loadtxt(fpath+'/files/mayor_completeness.txt', 
		usecols = (0,1,2), unpack = True)
	P_obs, Msini_obs = np.loadtxt(fpath+'/files/Mayor2011.csv', usecols=(2,5), 
			unpack = True, skiprows = 1, delimiter = ',')
	ID = np.loadtxt(fpath+'/files/Mayor2011.csv', usecols = (1,), skiprows = 1, 
			dtype = str, delimiter = ',')
	
	P_1D= np.unique(P)
	Msini_1D= np.unique(Msini)

	X, Y = np.meshgrid(P_1D, Msini_1D)
	Z= Completeness.reshape(X.shape)
	
	obs= {'xvar':P_obs, 'yvar':Msini_obs, 'starID':ID}
	
	survey= {'xvar':P_1D, 'yvar':Msini_1D, 'eff_2D':Z.T}
	obs['nstars']= 822
	
	return obs, survey

def _Howard2010():
	# this function is superseded by the Mayor2011 completeness
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