import numpy as np

'''
NOTE: Kepler files from Q16-epos.py
cp ~/Kepler/npz/completeness.Q16.epos.*.npz files/
cp ~/Kepler/npz/KOI.Q16.epos.*.npz files/

#cp ~/Kepler/npdict/Q16.occ.Rp-P.sptype.*.npdict files/
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
	KOI= np.load('files/KOI.Q16.epos.{}.npz'.format(subsample))
	return KOI['P'], KOI['Rp'], KOI['KID']

def eff_Q16(subsample='all'):
	eff= np.load('files/completeness.Q16.epos.{}.npz'.format(subsample))
	return eff['x'], eff['y'], eff['fsnr']
