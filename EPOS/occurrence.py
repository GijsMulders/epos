import numpy as np
from scipy import interpolate
import multiprocessing
from functools import partial

from EPOS.population import periodradius

def all(epos):
	if hasattr(epos,'occurrence'):
		planets(epos)
		if 'bin' in epos.occurrence:
			binned(epos)
		if 'xzoom' in epos.occurrence and 'yzoom' in epos.occurrence:
			zoomed(epos)
		
		# posterior per bin
		if epos.Prep and epos.populationtype is 'parametric' and (not epos.Multi) \
				and 'bin' in epos.occurrence:
			parametric(epos)
	else:
		print 'No bins for calculating occurrence rate, did you use epos.set_bins() ?'

def planets(epos, Log=False):
	if not epos.Range: epos.set_ranges()
	
	''' Interpolate occurrence for each planet (on log scale) '''
	focc= epos.occurrence
	
	print '\nInterpolating planet occurrence per bin'
	#with np.errstate(divide='ignore'): 
	#print '{}x{}=?={}'.format(epos.eff_xvar.shape, epos.eff_yvar.shape, epos.completeness.shape)
	
	if Log:
		# does not seem to make a big difference
		pl_comp_log= interpolate.RectBivariateSpline(np.log10(epos.eff_xvar), 
			np.log10(epos.eff_yvar), np.log10(epos.completeness))
		completeness= 10.**pl_comp_log(np.log10(epos.obs_xvar), np.log10(epos.obs_yvar), grid=False) 
		#print completeness
	else:
		# if zeros in detection efficiency?
		pl_comp= interpolate.RectBivariateSpline(epos.eff_xvar, 
			epos.eff_yvar, epos.completeness)
		completeness= pl_comp(epos.obs_xvar, epos.obs_yvar, grid=False)
	focc['planet']={}
	focc['planet']['completeness']= completeness
	focc['planet']['occ']= 1./completeness/epos.nstars
	#print epos.planet_occurrence

def binned(epos):	
	focc= epos.occurrence

	# set y bin
	if epos.MassRadius:
		focc['bin']['y']= epos.MR(focc['bin']['y in'])[0]
		for i, ybin in enumerate(focc['bin']['y']):
			if ybin[0]>ybin[-1]: focc['bin']['y'][i]= focc['bin']['y'][i][::-1] 
	else:
		focc['bin']['y']= focc['bin']['y in']

	# calc n, i, occ, err, xc, yc, dlnx, dlny from x & y
	_occ_per_bin(epos, focc['bin'])

def zoomed(epos):
	focc= epos.occurrence

	print '\n  x zoom bins'
	_occ_per_bin(epos, focc['xzoom'])
	print '\n  y zoom bins'
	_occ_per_bin(epos, focc['yzoom'])

def _occ_per_bin(epos, foccbin):	
	_occ, _n, _inbin, _xc, _yc, _dlnx, _dlny=[],[],[],[],[],[],[]
	
	for xbin, ybin in zip(foccbin['x'],foccbin['y']):
		inbinx= (xbin[0]<=epos.obs_xvar) & (epos.obs_xvar<xbin[1])
		inbiny= (ybin[0]<=epos.obs_yvar) & (epos.obs_yvar<ybin[1])
		inbin= inbinx & inbiny
		
		_inbin.append(inbin)
		_n.append(inbin.sum())
		_occ.append(epos.occurrence['planet']['occ'][inbin].sum())
		
		print '  x: [{:.3g},{:.3g}], y: [{:.2g},{:.2g}], n={}, occ={:.2g}'.format(
			xbin[0],xbin[-1], ybin[0],ybin[-1], _n[-1], _occ[-1])
		
		_xc.append(np.sqrt(xbin[0])*np.sqrt(xbin[-1]) )
		_yc.append(np.sqrt(ybin[0])*np.sqrt(ybin[-1]) )
		_dlnx.append(np.log(xbin[-1]/xbin[0]))
		_dlny.append(np.log(ybin[-1]/ybin[0]))

	
	foccbin['n']= np.array(_n)
	foccbin['i']= np.array(_inbin)
	foccbin['occ']= np.array(_occ)
	foccbin['err']= foccbin['occ']/np.where(
							foccbin['n']>0,np.sqrt(foccbin['n']),1.)
	
	foccbin['xc']= np.array(_xc)
	foccbin['yc']= np.array(_yc)
	foccbin['dlnx']= np.array(_dlnx)
	foccbin['dlny']= np.array(_dlny)
	
def parametric(epos):
	assert epos.Prep and epos.populationtype is 'parametric' and (not epos.Multi)
	focc= epos.occurrence
	
	''' loop over all bins '''
	_eta, _gamma, _area= [], [], []
	_pos, _sigp, _sign=[], [], []
	
	print '\n  posterior per bin'
	for xbin, ybin in zip(focc['bin']['x'],focc['bin']['y in']):
		_area.append(np.log(xbin[1]/xbin[0])*np.log(ybin[1]/ybin[0]))
		_, pdf, _, _= periodradius(epos, Init=True, xbin=xbin, ybin=ybin)
		_gamma.append(np.average(pdf))
		_eta.append(_gamma[-1]*_area[-1])

		print '  x: [{:.3g},{:.3g}], y: [{:.2g},{:.2g}], area={:.2f}, eta={:.2g}'.format(
			xbin[0],xbin[-1], ybin[0],ybin[-1], _area[-1], _eta[-1])

		''' Posterior?'''
		if hasattr(epos, 'samples'):
			posterior= []
			if epos.Parallel:
				pool = multiprocessing.Pool()
				one_arg_func= partial(_posterior, epos, xbin=xbin, ybin=ybin)
				posterior= pool.map(one_arg_func, epos.samples)
			else:
				for sample in epos.samples:
					_, pdf, _, _= periodradius(epos, fpara=sample, xbin=xbin, ybin=ybin)
					posterior.append(np.average(pdf))

			pos= np.percentile(posterior, [16, 50, 84])
			_pos.append(pos[1])
			_sigp.append(pos[2]-pos[1])
			_sign.append(pos[1]-pos[0])
		
			print '  gamma= {:.1%} +{:.1%} -{:.1%}'.format(_pos[-1],_sigp[-1],_sign[-1])
			print '  eta= {:.1%} +{:.1%} -{:.1%}'.format(
				_pos[-1]*_area[-1],_sigp[-1]*_area[-1],_sign[-1]*_area[-1])

	focc['bin']['area']= np.array(_area)
	focc['bin']['gamma0']= np.array(_gamma)
	focc['bin']['eta0']= np.array(_eta)	
	
	if hasattr(epos, 'samples'):
		focc['bin']['gamma']= np.array(_pos)
		focc['bin']['gamma+']= np.array(_sigp)
		focc['bin']['gamma-']= np.array(_sign)
		focc['bin']['eta']= focc['bin']['gamma']*focc['bin']['area']
		focc['bin']['eta+']= focc['bin']['gamma+']*focc['bin']['area']
		focc['bin']['eta-']= focc['bin']['gamma-']*focc['bin']['area']

def _posterior(epos, sample, xbin, ybin):
	_, pdf, _, _= periodradius(epos, fpara=sample, xbin=xbin, ybin=ybin)
	return np.average(pdf)
	
	