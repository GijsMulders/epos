import numpy as np
from scipy import interpolate
from run import _pdf

def all(epos):
	planets(epos)
	if epos.Prep and epos.populationtype is 'parametric' and (not epos.Multi) :
		parametric(epos)

def planets(epos, Log=False):
	if not epos.Range: epos.set_ranges()
	
	''' Interpolate occurrence for each planet (on log scale) '''
	focc= epos.occurrence
	
	print '\nInterpolating planet occurrence:'
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
	
	''' Occurrence per bin '''
	_occ, _n, _inbin=[],[],[]
	for xbin, ybin in zip(focc['bin']['x'],focc['bin']['y']):
		inbinx= (xbin[0]<=epos.obs_xvar) & (epos.obs_xvar<xbin[1])
		inbiny= (ybin[0]<=epos.obs_yvar) & (epos.obs_yvar<ybin[1])
		inbin= inbinx & inbiny
		
		_inbin.append(inbin)
		_n.append(inbin.sum())
		_occ.append(focc['planet']['occ'][inbin].sum())
		
		print 'x: {}, y: {}, n={}, occ={:.2g}'.format(xbin, ybin, _n[-1], _occ[-1])
	
	focc['bin']['n']= np.array(_n)
	focc['bin']['i']= np.array(_inbin)
	focc['bin']['occ']= np.array(_occ)
	#print type(np.array(_inbin))

def parametric(epos):
	assert epos.Prep and epos.populationtype is 'parametric' and (not epos.Multi)
	focc= epos.occurrence
	
	''' loop over all bins '''
	_eta, _gamma, _area= [], [], []
	_pos, _sigp, _sign=[], [], []
	
	for xbin, ybin in zip(focc['bin']['x'],focc['bin']['y']):
		_area.append(np.log(xbin[1]/xbin[0])*np.log(ybin[1]/ybin[0]))
		_, pdf, _, _= _pdf(epos, Init=True, xbin=xbin, ybin=ybin)
		_gamma.append(np.average(pdf)*epos.scale)
		_eta.append(_gamma[-1]*_area[-1])
		print 'x: {}, y: {}, area={:.2f}, eta={:.2g}'.format(
			xbin, ybin, _area[-1], _eta[-1])
		
		''' Posterior?'''
		if hasattr(epos, 'samples'):
			posterior= []
			for sample in epos.samples:
				_, pdf, _, _= _pdf(epos, fpara=sample, xbin=xbin, ybin=ybin)
				posterior.append(np.average(pdf)*epos.scale )

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
