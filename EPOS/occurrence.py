import numpy as np
from scipy import interpolate

def planets(epos, Log=False):
	if not epos.Range: epos.set_ranges()
	
	''' Interpolate occurrence for each planet (on log scale) '''
	focc= epos.occurrence
	
	print '\nInterpolating planet occurrence:'
	#with np.errstate(divide='ignore'): 
	print '{}x{}=?={}'.format(epos.eff_xvar.shape, epos.eff_yvar.shape, epos.completeness.shape)
	
	if Log:
		# does not seem to make a big difference
		pl_comp_log= interpolate.RectBivariateSpline(np.log10(epos.eff_xvar), 
			np.log10(epos.eff_yvar), np.log10(epos.completeness))
		completeness= 10.**pl_comp_log(np.log10(epos.obs_xvar), np.log10(epos.obs_yvar), grid=False) 
		print completeness
	else:
		# if zeros in detection efficiency?
		pl_comp= interpolate.RectBivariateSpline(epos.eff_xvar, 
			epos.eff_yvar, epos.completeness)
		completeness= pl_comp(epos.obs_xvar, epos.obs_yvar, grid=False)
	focc['planet']= 1./completeness/epos.nstars
	#print epos.planet_occurrence
	
	''' Occurrence per bin '''
	_occ, _n, _inbin=[],[],[]
	for xbin, ybin in zip(focc['bin']['x'],focc['bin']['y']):
		inbinx= (xbin[0]<=epos.obs_xvar) & (epos.obs_xvar<xbin[1])
		inbiny= (ybin[0]<=epos.obs_yvar) & (epos.obs_yvar<ybin[1])
		inbin= inbinx & inbiny
		
		_inbin.append(inbin)
		_n.append(inbin.sum())
		_occ.append(focc['planet'][inbin].sum())
		
		print 'x: {}, y: {}, n={}, occ={:.2g}'.format(xbin, ybin, _n[-1], _occ[-1])
	
	focc['bin']['n']= np.array(_n)
	focc['bin']['i']= np.array(_inbin)
	focc['bin']['occ']= np.array(_occ)
	
	
	
	