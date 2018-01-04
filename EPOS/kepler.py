''' 
This module contains helper functions to load kepler survey data into EPOS
'''
import numpy as np
import os.path

try:
	from astropy.table import Table
except ImportError:
	print '\nWarning: Could not import astropy.table'

'''
NOTE: Kepler files from Q16-epos.py
cp ~/Kepler/npz/completeness.Q16.epos.*.npz files/
cp ~/Kepler/npz/KOI.Q16.epos.*.npz files/

#cp ~/Kepler/npdict/Q16.occ.Rp-P.sptype.*.npdict files/
'''


def obs_Q16(subsample='all'):
	'''
	returns the Q16 exoplanet catalogue as numpy arrays
	
	Returns:
		P(np.array):   Orbital period in days
		Rp(np.array):  Planet radius in earth radii'
		KID(np.array): star identifier (for multis)'
	'''
	KOI= np.load('files/KOI.Q16.epos.{}.npz'.format(subsample))
	return KOI['P'], KOI['Rp'], KOI['KID']

def eff_Q16(subsample='all'):
	'''
	returns the Q16 survey detection efficiency as a 2D matrix
	
	Returns:
		P(np.array):   Orbital period in days
		Rp(np.array):  Planet radius in earth radii'
		eff(np.array):2D matrix of detection efficiency [0-1]
	'''
	eff= np.load('files/completeness.Q16.epos.{}.npz'.format(subsample))
	return eff['x'], eff['y'], eff['fsnr']

def dr25(subsample='all', score=0.9, Huber=False):
	'''
	Generates Kepler DR25 planet population and detection efficiency
	
	Args:
		subsample(str):	Subsample, choose from 'all', 'M', 'K', 'G', or 'F'
		
	Returns:
		tuple: two dictionaries
		
		obs(dict):
			xvar(np.array of float): orbital period
			yvar(np.array of float): planet radius
			starID(np.array): stellar ID
			nstars(int): number of stars surveyed
		survey(dict):
			the grid is in xvar,yvar, the detection efficiency in eff_2D 
	'''
	fkoi= 'files/q1_q17_dr25_koi.npz'
	if os.path.isfile(fkoi):
		print '\nLoading planets from {}'.format(fkoi)
		koi= np.load(fkoi)
	else:
		print '\nReading planets from IPAC file' 
		try:
			ipac=Table.read('files/q1_q17_dr25_koi.tbl',format='ipac')
		except NameError:
			raise ImportError('You need to install astropy for this step')
		koi= {key: np.array(ipac[key]) for key in ['kepid','koi_prad','koi_period',
				'koi_steff', 'koi_slogg','koi_pdisposition','koi_score'] }
		np.savez(fkoi, **koi)

	''' Select reliable candidates, remove giant stars'''
	if Huber:
		isdwarf= stars['logg'] > 1./4.671 * \
					np.arctan((stars['Teff']-6300.)/-67.172)+ 3.876
		isgiant= stars['logg'] < np.where(stars['Teff']>5000, 
									13.463-0.00191*stars['Teff'],3.9)
		issubgiant= ~isdwarf & ~isgiant
		suffix='huber'
	else:
		isdwarf= koi['koi_slogg']>4.2
		suffix='logg42'
	
	iscandidate= koi['koi_pdisposition']=='CANDIDATE'
	isreliable= koi['koi_score']>score # removes rolling band, ~ 500 planets
	print '  {}/{} dwarfs'.format(isdwarf.sum(), isdwarf.size)
	print '  {} candidates, {} false positives'.format((isdwarf&iscandidate).sum(), 
				(isdwarf&~iscandidate).sum() )
	print '  {}+{} with score > {:.2f}'.format((isdwarf&iscandidate&isreliable).sum(), 
				(isdwarf&~iscandidate&isreliable).sum(),score )

	isall= islogg42 & isreliable
	#isall= isdwarf & isreliable # reliability score cuts out false positives
#	isall= isdwarf & iscandidate & isreliable

	slice={'all':isall}
	for spT, Tmin, Tmax in zip(['M','K','G','F'],
		[2400, 3865, 5310, 5980], [3865, 5310, 5980, 7320]):
		slice[spT]= isall & (Tmin<koi['koi_steff']) & (koi['koi_steff']<=Tmax)
	
	obs= {'xvar':koi['koi_period'][slice[subsample]],
		'yvar':koi['koi_prad'][slice[subsample]], 
		'starID':koi['kepid'][slice[subsample]]}
	
	# from dr25_epos.py in 
	eff= np.load('files/completeness.dr25.{}.{}.npz'.format(subsample, suffix))
	#eff= np.load('files/det_eff.dr25.{}.npz'.format(subsample))
	survey= {'xvar':eff['P'], 'yvar':eff['Rp'], 'eff_2D':eff['fsnr'], 
			'Mstar': eff['Mst'], 'Rstar':eff['Rst']}
	obs['nstars']= eff['n']
	
	return obs, survey