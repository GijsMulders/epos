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

def dr25(subsample='all', score=0.9, Huber=True, Vetting=False):
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

	''' Remove giant stars'''
	if Huber:
		isdwarf= koi['koi_slogg'] > 1./4.671 * \
					np.arctan((koi['koi_steff']-6300.)/-67.172)+ 3.876
		isgiant= koi['koi_slogg'] < np.where(koi['koi_steff']>5000, 
									13.463-0.00191*koi['koi_steff'],3.9)
		issubgiant= ~isdwarf & ~isgiant
		suffix='huber'
	else:
		isdwarf= koi['koi_slogg']>4.2
		suffix='logg42'
	
	''' Select reliable candidates '''
	iscandidate= koi['koi_pdisposition']=='CANDIDATE'
	isreliable= koi['koi_score']>score # removes rolling band, ~ 500 planets
	print '  {}/{} dwarfs'.format(isdwarf.sum(), isdwarf.size)
	print '  {} candidates, {} false positives'.format((isdwarf&iscandidate).sum(), 
				(isdwarf&~iscandidate).sum() )
	print '  {}+{} with score > {:.2f}'.format((isdwarf&iscandidate&isreliable).sum(), 
				(isdwarf&~iscandidate&isreliable).sum(),score )

	isall= isdwarf & isreliable # reliability score cuts out false positives
#	isall= isdwarf & iscandidate & isreliable

	''' Select a spectral type sub-smaple'''
	if subsample=='all':
		slice=isall
	elif subsample in ['M','K','G','F']:
		Teff={'M':[2400, 3865], 'K':[3865, 5310], 'G':[5310, 5980], 'F':[5980, 7320]}
		slice= isall & (Teff['subsample'][0]<koi['koi_steff']) \
						& (koi['koi_steff']<=Teff['subsample'][1])
	elif subsample[0] is 'T':
		Tmin=int(subsample[1:])-250
		Tmax=int(subsample[1:])+250
		slice= isall & (Tmin<koi['koi_steff']) & (koi['koi_steff']<=Tmax)
	else:
		raise ValueError('Subsample {} not recognized'.format(subsample))
	#if Huber: slice['subgiants']=issubgiant
	
	obs= {'xvar':koi['koi_period'][slice],
		'yvar':koi['koi_prad'][slice], 
		'starID':koi['kepid'][slice]}
	
	''' Load pre-calculated detection efficiencies '''
	# from dr25_epos.py
	sfile= 'dwarfs' if subsample=='all' else subsample
	eff= np.load('files/completeness.dr25.{}.{}.npz'.format(sfile, suffix))
	#eff= np.load('files/det_eff.dr25.{}.npz'.format(subsample))
	survey= {'xvar':eff['P'], 'yvar':eff['Rp'], 'eff_2D':eff['fsnr'], 
			'Mstar': eff['Mst'], 'Rstar':eff['Rst']}
	obs['nstars']= eff['n']
  	
	''' Add vetting completeness '''
	if Vetting:
		if subsample=='all' and score == 0.9:
			X,Y= np.meshgrid(eff['P'], eff['Rp'], indexing='ij')
			#vet_2D= fbpl2d( (X,Y), 0.9, 46., -0.053, -0.37, 5.6, 0.19, -2.3)
			vet_2D= fbpl2d( (X,Y), 0.9, 46., -0.053, -0.37, 5.6, 0.19, 0.19)
			survey['vet_2D']= vet_2D
			assert vet_2D.shape == eff['fsnr'].shape
		else:
			print 'no vetting completeness for {} with score={}'.format(score, subsample)
	
	return obs, survey
	
def fbpl2d((x,y), a, b, c, d, e, f, g):
	bpl= a* (x/b)**np.where(x<b, c, d) * (y/e)**np.where(y<e, f,g)
	return np.maximum(0.2, np.minimum(bpl, 1))