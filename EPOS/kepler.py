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
	print 'dr25(): returns the DR25 planet list and '
	print '        survey detection efficiency as a 2D matrix'
	print '  P:   Orbital period in days'
	print '  Rp:  Planet radius in earth radii'
	print '  eff: 2D matrix of detection efficiency [0-1]'

def obs_Q16(subsample='all'):
	KOI= np.load('files/KOI.Q16.epos.{}.npz'.format(subsample))
	return KOI['P'], KOI['Rp'], KOI['KID']

def eff_Q16(subsample='all'):
	eff= np.load('files/completeness.Q16.epos.{}.npz'.format(subsample))
	return eff['x'], eff['y'], eff['fsnr']

def dr25(subsample='all', score=0.9):

	from astropy.table import Table
	print '\nReading planets from IPAC file' 
	ipac=Table.read('files/q1_q17_dr25_koi.tbl',format='ipac')
	isdwarf= ipac['koi_slogg']>4.2
	iscandidate= ipac['koi_pdisposition']=='CANDIDATE'
	# koi_score
	isreliable= ipac['koi_score']>score # removes rolling band, ~ 500 planets
	print '  {}/{} dwarfs'.format(isdwarf.sum(), isdwarf.size)
	print '  {} candidates, {} false positives'.format((isdwarf&iscandidate).sum(), 
				(isdwarf&~iscandidate).sum() )
	print '  {}+{} with score > {:.2f}'.format((isdwarf&iscandidate&isreliable).sum(), 
				(isdwarf&~iscandidate&isreliable).sum(),score )

	isall= isdwarf & isreliable # reliability score cuts out false positives
#	isall= isdwarf & iscandidate & isreliable

	slice={'all':isall}
	for spT, Tmin, Tmax in zip(['M','K','G','F'],
		[2400, 3865, 5310, 5980], [3865, 5310, 5980, 7320]):
		slice[spT]= isall & (Tmin<ipac['koi_steff']) & (ipac['koi_steff']<=Tmax)
	
	obs= {'xvar':ipac['koi_period'][slice[subsample]],
		'yvar':ipac['koi_prad'][slice[subsample]], 
		'starID':ipac['kepid'][slice[subsample]]}
	
	# from dr25_epos.py in 
	eff= np.load('files/completeness.dr25.{}.npz'.format(subsample))
	#eff= np.load('files/det_eff.dr25.{}.npz'.format(subsample))
	survey= {'xvar':eff['P'], 'yvar':eff['Rp'], 'eff_2D':eff['fsnr'], 
			'Mstar': eff['Mst'], 'Rstar':eff['Rst']}
	obs['nstars']= eff['n']
	
	return obs, survey