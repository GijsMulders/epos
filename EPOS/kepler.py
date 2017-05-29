import numpydict as npdict
import pickle
import cgs

'''
NOTE: Kepler files from Q16-epos.py
cp ~/Kepler/pickle/completeness.Q16.epos.*.pickle files/
cp ~/Kepler/npdict/Q16.occ.Rp-P.npdict files/Q16.occ.Rp-P.all.npdict
cp ~/Kepler/npdict/Q16.occ.Rp-P.sptype.*.npdict files/

TODO: switch to better format, rename to DR24 etc.
'''

def readme():
	#import clearscreen # works only on loading, not on execution
	print '\nThis module contains helper functions to load kepler survey data into EPOS:'
	print
	print 'obs_Q16(): returns the Q16 exoplanet catalogue as numpy arrays'
	print '  P:   Orbital period in days'
	print '  Rp:  Planet radius in earth radii'
	print '  KID: star identifier (for multis, not implemented)'
	print
	print 'eff_Q16(): returns the Q16 survey detection efficiency as a 2D matrix'
	print '  P:   Orbital period in days'
	print '  Rp:  Planet radius in earth radii'
	print '  eff: 2D matrix of detection efficiency [0-1]'
	print

def obs_Q16(subsample='all'):
	KOI= npdict.load('files/Q16.occ.Rp-P.sptype.{}.npdict'.format(subsample))
	return KOI['P']/cgs.day, KOI['Rp']/cgs.Rearth, KOI['KID']

def eff_Q16(subsample='all'):
	fname= 'files/completeness.Q16.epos.{}.pickle'.format(subsample)
	with open(fname,'r') as f: 
		grid= pickle.load(f)
	
	return grid['x'], grid['y'], grid['completeness'] 