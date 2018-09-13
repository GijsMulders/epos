import numpy as np
from scipy import interpolate
import multiprocessing
from functools import partial

from EPOS.population import periodradius

def all(epos):
	if hasattr(epos,'occurrence'):
		planets(epos)

	if hasattr(epos,'pfm'):
		models(epos)

	if hasattr(epos,'occurrence'):
		if 'bin' in epos.occurrence:
			binned(epos)
		if 'xzoom' in epos.occurrence and 'yzoom' in epos.occurrence:
			zoomed(epos)
		
		# posterior per bin
		if epos.Prep and epos.Parametric and (not epos.Multi) \
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
	focc['planet']['xvar']= epos.obs_xvar
	focc['planet']['yvar']= epos.obs_yvar	
	focc['planet']['completeness']= completeness
	focc['planet']['occ']= 1./completeness/epos.nstars
	#print epos.planet_occurrence

def models(epos, Log=False):
	if not epos.Range: epos.set_ranges()

	if not hasattr(epos,'occurrence'):
		epos.occurrence={}
	
	''' Interpolate occurrence for each planet (on log scale) '''
	focc= epos.occurrence
	print '\nInterpolating model planet occurrence'
	
	pfm=epos.pfm
	if not 'R' in pfm:
		pfm['R'], _= epos.MR(pfm['M'])
	
	if Log:
		# does not seem to make a big difference
		pl_comp_log= interpolate.RectBivariateSpline(np.log10(epos.eff_xvar), 
			np.log10(epos.eff_yvar), np.log10(epos.completeness))
		completeness= 10.**pl_comp_log(np.log10(pfm['P']), np.log10(pfm['R']), grid=False) 
		#print completeness
	else:
		# if zeros in detection efficiency?
		pl_comp= interpolate.RectBivariateSpline(epos.eff_xvar, 
			epos.eff_yvar, epos.completeness)
		completeness= pl_comp(pfm['P'], pfm['R'], grid=False)
		
	focc['model']={}
	focc['model']['completeness']= completeness
	#focc['model']['occ']= 1./completeness/epos.nstars
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

	print '\n  Observed Planets'
	# calc n, i, occ, err, xc, yc, dlnx, dlny from x & y
	_occ_per_bin(epos, focc['bin'])
	
	if 'model' in focc:
		eta= epos.fitpars.getpps() if hasattr(epos, 'fitpars') else 1.
		print '\n  Modeled planets, eta= {:.3g}'.format(eta)
		_model_occ_per_bin(epos, focc['bin'], focc['model'], weights=eta/epos.pfm['ns'])
		focc['model']['eta']= eta

def zoomed(epos):
	''' Occurrence (inverse detection efficiency) along x,y axis.'''
	focc= epos.occurrence

	print '\n  x zoom bins'
	_occ_per_bin(epos, focc['xzoom'])
	print '\n  y zoom bins'
	_occ_per_bin(epos, focc['yzoom'])

def _occ_per_bin(epos, foccbin):
	''' Planet occurrence (inverse detection efficiency) per bin '''	
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

def _model_occ_per_bin(epos, foccbin, foccmodel, weights=None):	
	_occ, _n, _inbin, _xc, _yc, _dlnx, _dlny=[],[],[],[],[],[],[]
	
	for xbin, ybin in zip(foccbin['x'],foccbin['y']):
		inbinx= (xbin[0]<=epos.pfm['P']) & (epos.pfm['P']<xbin[1])
		inbiny= (ybin[0]<=epos.pfm['R']) & (epos.pfm['R']<ybin[1])
		inbin= inbinx & inbiny
		
		_inbin.append(inbin)
		_n.append(inbin.sum())
		_occ.append(weights* inbin.sum())
		
		print '  x: [{:.3g},{:.3g}], y: [{:.2g},{:.2g}], n={}, occ={:.2g}'.format(
			xbin[0],xbin[-1], ybin[0],ybin[-1], _n[-1], _occ[-1])
		
		_xc.append(np.sqrt(xbin[0])*np.sqrt(xbin[-1]) )
		_yc.append(np.sqrt(ybin[0])*np.sqrt(ybin[-1]) )
		_dlnx.append(np.log(xbin[-1]/xbin[0]))
		_dlny.append(np.log(ybin[-1]/ybin[0]))

	_foccbin= foccmodel['bin']= {}
	_foccbin['n']= np.array(_n)
	_foccbin['i']= np.array(_inbin)
	_foccbin['occ']= np.array(_occ)
	_foccbin['err']= _foccbin['occ']/np.where(
							_foccbin['n']>0,np.sqrt(_foccbin['n']),1.)
	
	_foccbin['xc']= np.array(_xc)
	_foccbin['yc']= np.array(_yc)
	_foccbin['dlnx']= np.array(_dlnx)
	_foccbin['dlny']= np.array(_dlny)
	
	_foccbin['x']= foccbin['x']
	_foccbin['y']= foccbin['y']
	
def parametric(epos):
	''' Calculates the occurrence rate per bin from the parametric model, 
	with uncertainties if samples from an MCMC chain are available
	'''
	assert epos.Prep and epos.Parametric and (not epos.Multi)
	focc= epos.occurrence
	
	''' loop over all pre-defined bins '''
	print '\n  posterior per bin'
	xbins= focc['bin']['x']
	ybins= focc['bin']['y in']
	eta, gamma, area, pos, sigp, sign= _posterior_per_bin(epos, xbins, ybins, Verbose=True)

	focc['bin']['area']= np.array(area)
	focc['bin']['gamma0']= np.array(gamma)
	focc['bin']['eta0']= np.array(eta)	
	
	if hasattr(epos, 'samples'):
		focc['bin']['gamma']= np.array(pos)
		focc['bin']['gamma+']= np.array(sigp)
		focc['bin']['gamma-']= np.array(sign)
		focc['bin']['eta']= focc['bin']['gamma']*focc['bin']['area']
		focc['bin']['eta+']= focc['bin']['gamma+']*focc['bin']['area']
		focc['bin']['eta-']= focc['bin']['gamma-']*focc['bin']['area']

	''' bin for normalization '''
	if hasattr(epos.fitpars, 'normkeyx') and hasattr(epos.fitpars, 'normkeyy'):
		print '\n  normalization per unit of ln area'
		dw= 1.01
		xnorm= epos.fitpars.get(epos.fitpars.normkeyx)
		ynorm= epos.fitpars.get(epos.fitpars.normkeyy)
		xbin= [xnorm/dw, xnorm*dw]
		ybin= [ynorm/dw, ynorm*dw]
		eta, gamma, _, gamma_fit, gamma_p, gamma_n=  _posterior_per_bin(epos, [xbin],[ybin])
		print '  x={:.2g}, y={:.2g}, gamma= {:.2g}'.format(xnorm, ynorm, gamma[0])
		if len(gamma_fit)>0:
			print '  gamma= {:.2g} +{:.2g} - {:.2g}'.format(gamma_fit[0], gamma_p[0], gamma_n[0])

def _posterior(epos, sample, xbin, ybin):
	_, pdf, _, _= periodradius(epos, fpara=sample, xbin=xbin, ybin=ybin)
	return np.average(pdf)

def _posterior_per_bin(epos, xbins, ybins, Verbose=True):
	eta, gamma, area= [], [], []
	pos, sigp, sign=[], [], []
	
	for xbin, ybin in zip(xbins,ybins):
		area.append(np.log(xbin[1]/xbin[0])*np.log(ybin[1]/ybin[0]))
		_, pdf, _, _= periodradius(epos, Init=True, xbin=xbin, ybin=ybin)
		gamma.append(np.average(pdf))
		eta.append(gamma[-1]*area[-1])

		if Verbose:
			print '  x: [{:.3g},{:.3g}], y: [{:.2g},{:.2g}], area={:.2f}, eta_0={:.2g}'.format(
				xbin[0],xbin[-1], ybin[0],ybin[-1], area[-1], eta[-1])

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

			#pos= np.percentile(posterior, [16, 50, 84])
			perc= np.percentile(posterior, [2.3, 15.9, 50., 84.1, 97.7])
			pos.append(perc[2])
			sigp.append(perc[3]-perc[2])
			sign.append(perc[2]-perc[1])
			sig2p= (perc[4]-perc[2])
			sig2n= (perc[2]-perc[0])
		
			if Verbose:
				print '  gamma= {:.1%} +{:.1%} -{:.1%}'.format(pos[-1],sigp[-1],sign[-1])
				print '  eta= {:.1%} +{:.1%} -{:.1%}'.format(
					pos[-1]*area[-1],sigp[-1]*area[-1],sign[-1]*area[-1])
				#print '  gamma 2sig= {:.1%} +{:.1%} -{:.1%}'.format(_pos[-1],sig2p,sig2n)

	return eta, gamma, area, pos, sigp, sign
	