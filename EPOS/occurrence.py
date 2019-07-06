import numpy as np
from scipy import interpolate
import multiprocessing
from functools import partial

import shapely.geometry

from EPOS.population import periodradius
import EPOS.analytics

def all(epos, BinnedMetric=False, MLE=False):
	if hasattr(epos,'occurrence'):
		planets(epos)

	if hasattr(epos,'pfm'):
		models(epos)

	if hasattr(epos,'occurrence'):
		if 'bin' in epos.occurrence:
			binned(epos, MLE=MLE)
		if 'poly' in epos.occurrence:
			poly_binned(epos)
		if 'xzoom' in epos.occurrence and 'yzoom' in epos.occurrence:
			zoomed(epos, MLE=MLE)
		
		# posterior per bin
		if epos.Prep and epos.Parametric and (not epos.Multi) \
				and 'bin' in epos.occurrence:
			parametric(epos)

			''' BIC/AIC with chi^2'''
			if BinnedMetric and 'yzoom' in epos.occurrence or 'xzoom' in epos.occurrence:
				print '\n  Binned occurrence rate metrics'
				binned_metric(epos)
	else:
		print 'No bins for calculating occurrence rate, did you use epos.set_bins() ?'
	
def planets(epos, Log=False):
	''' 
	Calculate the intrinsic occurrence of each observed planet by interpolating the survey completeness
	'''
	if not epos.Range: epos.set_ranges()
	
	''' Interpolate occurrence for each planet (on log scale) '''
	focc= epos.occurrence
	
	print '\nInterpolating planet occurrence'
	#with np.errstate(divide='ignore'): 
	#print '{}x{}=?={}'.format(epos.eff_xvar.shape, epos.eff_yvar.shape, epos.completeness.shape)
	
	completeness= _interpolate_occ(epos, epos.obs_xvar, epos.obs_yvar, Log=Log)

	focc['planet']={}
	focc['planet']['xvar']= epos.obs_xvar
	focc['planet']['yvar']= epos.obs_yvar	
	focc['planet']['completeness']= completeness
	focc['planet']['occ']= 1./completeness/epos.nstars
	#print epos.planet_occurrence

def models(epos, Log=False):
	'''
	Interpolate the survey completeness for each modeled planet
	'''
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
	focc['model']['eta']= epos.fitpars.getpps() if hasattr(epos, 'fitpars') else 1.
	#focc['model']['occ']= focc['model']['eta']/completeness/epos.nstars
	#print epos.planet_occurrence

def binned(epos, Log=False, MLE=False):
	'''
	Calculate the planet occurrence rate per bin using inverse detection efficiency.
	(these are the bins from the detection efficiency grid)
	'''	
	focc= epos.occurrence

	# set y bin
	if epos.MassRadius:
		focc['bin']['y']= epos.MR(focc['bin']['y in'])[0]
		for i, ybin in enumerate(focc['bin']['y']):
			if ybin[0]>ybin[-1]: 
				print '\nWarning: observed bins are not increasing in size'
				print '   {} > {}'.format(ybin[0],ybin[-1])
				focc['bin']['y'][i]= focc['bin']['y'][i][::-1] 
	else:
		focc['bin']['y']= focc['bin']['y in']

	print '\n  Observed Planets'
	# calc n, i, occ, err, xc, yc, dlnx, dlny from x & y
	_occ_per_bin(epos, focc['bin'], Log=Log)
	if MLE: _occ_per_bin_MLE(epos, focc['bin'], Log=Log)

	if 'model' in focc:
		eta= focc['model']['eta']
		print '\n  Modeled planets, eta= {:.3g}'.format(eta)
		_model_occ_per_bin(epos, focc['bin'], focc['model'], weights=eta/epos.pfm['ns'])

def poly_binned(epos, Log=False):
	'''
	Calculate the planet occurrence rate per bin using inverse detection efficiency.
	(these are the polygonic bins)
	'''	
	focc= epos.occurrence

	# set y bin
	if epos.MassRadius:
		assert False, 'implement polygons for planet mass'

	print '\n  Observed Planets'
	# calc n, occ, err, dlnxy 
	_occ_per_polygon(epos, focc['poly'], Log=Log)
	
	if 'model' in focc:
		eta= focc['model']['eta']
		print '\n  Modeled planets, eta= {:.3g}'.format(eta)
		_model_occ_per_polygon(epos, focc['poly'], focc['model'], 
			weights=eta/epos.pfm['ns'])

def zoomed(epos, TestScore= True, MLE=False):
	''' Occurrence (inverse detection efficiency) along x,y axis.'''
	focc= epos.occurrence

	print '\n  x zoom bins'
	_occ_per_bin(epos, focc['xzoom'])
	if MLE: _occ_per_bin_MLE(epos, focc['xzoom'])
	print '\n  y zoom bins'
	_occ_per_bin(epos, focc['yzoom'])
	if MLE: _occ_per_bin_MLE(epos, focc['yzoom'])

def _occ_per_bin(epos, foccbin, Log=False):
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

	foccbin['xc']= np.array(_xc)
	foccbin['yc']= np.array(_yc)
	foccbin['dlnx']= np.array(_dlnx)
	foccbin['dlny']= np.array(_dlny)

	foccbin['n']= np.array(_n)
	foccbin['i']= np.array(_inbin)
	foccbin['occ']= np.array(_occ)
	foccbin['err']= foccbin['occ']/np.where(
							foccbin['n']>0,np.sqrt(foccbin['n']),1.)
	isUL= foccbin['n']==0
	if isUL.sum()>0:
		completeness_UL= _interpolate_occ(epos, foccbin['xc'][isUL], foccbin['yc'][isUL], Log=Log)
		foccbin['err'][isUL]= 1./completeness_UL/epos.nstars

	# quick MLE by using completenss at bin centre
	# completeness_MLE= _interpolate_occ(epos, foccbin['xc'], foccbin['yc'], Log=Log)
	# foccbin['occ_MLE']= 1.* foccbin['n']/epos.nstars /completeness_MLE
	# foccbin['err_MLE']= foccbin['occ_MLE']/np.where(
	# 	foccbin['n']>0,np.sqrt(foccbin['n']),1.)

def _occ_per_polygon(epos, foccpoly, Log=False):
	''' Planet occurrence (inverse detection efficiency) per bin '''	
	_occ, _n, _inbin, _xc, _yc, _dlnxy=[],[],[],[], [], []
	
	xyp= np.hstack([epos.obs_xvar[:,None],epos.obs_yvar[:,None]])

	for xycoords in foccpoly['coords']:
		shp=shapely.geometry.polygon.Polygon(xycoords)
		shp_log=shapely.geometry.polygon.Polygon(np.log(xycoords))

		points= map(shapely.geometry.Point, xyp)
		inpoly= map(shp.contains, points)

		logcenter= shp_log.centroid
		
		_inbin.append(inpoly)
		_n.append(np.sum(inpoly))
		_occ.append(epos.occurrence['planet']['occ'][inpoly].sum())
		_dlnxy.append(np.log(shp_log.area))

		_xc.append(np.exp(logcenter.x))
		_yc.append(np.exp(logcenter.y))

		print '  poly, n={}, occ={:.2g}, area={:.2g}'.format(_n[-1], _occ[-1], _dlnxy[-1])
		print '    xc= {:.2g}, yc={:.2g}'.format(_xc[-1], _yc[-1])
	
	foccpoly['xc']= np.array(_xc)
	foccpoly['yc']= np.array(_yc)	
	foccpoly['dlnxy']= np.array(_dlnxy)

	foccpoly['n']= np.array(_n)
	foccpoly['i']= np.array(_inbin)

	foccpoly['occ']= np.array(_occ)
	foccpoly['err']= foccpoly['occ']/np.where(
							foccpoly['n']>0,np.sqrt(foccpoly['n']),1.)
	
	isUL= foccpoly['n']==0
	if isUL.sum()>0:
		completeness_UL= _interpolate_occ(epos, foccpoly['xc'][isUL], foccpoly['yc'][isUL], Log=Log)
		foccpoly['err'][isUL]= 1./completeness_UL/epos.nstars 

def _occ_per_bin_MLE(epos, foccbin, Log=False, nres= 30):
	''' Planet occurrence (maximim likelihood estimate) per bin '''	
	_occ, _n, _nb, _comp, _xc, _yc, _dlnx, _dlny=[],[], [],[],[],[],[],[]
	
	# speedup by calculating completeness on fine grid once instead of per bin?
	# xgrid= 
	# ygrid= 
	# comp_fine= 

	for xbin, ybin in zip(foccbin['x'],foccbin['y']):
		inbinx= (xbin[0]<=epos.obs_xvar) & (epos.obs_xvar<xbin[1])
		inbiny= (ybin[0]<=epos.obs_yvar) & (epos.obs_yvar<ybin[1])
		_n.append(np.sum(inbinx & inbiny)) # planets

		# average completeness over bin
		xgrid= np.geomspace(xbin[0],xbin[-1], nres)
		ygrid= np.geomspace(ybin[0],ybin[-1], nres)
		completeness= _interpolate_occ(epos, xgrid, ygrid, Log=Log)

		#inbinx= (xbin[0]<=epos.X) & (epos.X<xbin[1])
		#inbiny= (ybin[0]<=epos.Y) & (epos.Y<ybin[1])
		#inbin= inbinx & inbiny
		#_nb.append(np.sum(inbinx & inbiny)) # bins
		#_comp.append(np.mean(epos.completeness[inbin]))

		#from scipy.stats.mstats import gmean
		#_comp.append(gmean(completeness))	
		_comp.append(np.mean(completeness))	
		_occ.append(1.*_n[-1]/epos.nstars/_comp[-1])
		
		print '  x: [{:.3g},{:.3g}], y: [{:.2g},{:.2g}], nb={}, occ={:.2g}'.format(
			xbin[0],xbin[-1], ybin[0],ybin[-1], nres*nres, _occ[-1])
		#_nb[-1]
		
		_xc.append(np.sqrt(xbin[0])*np.sqrt(xbin[-1]) )
		_yc.append(np.sqrt(ybin[0])*np.sqrt(ybin[-1]) )
		_dlnx.append(np.log(xbin[-1]/xbin[0]))
		_dlny.append(np.log(ybin[-1]/ybin[0]))

	foccbin['xc']= np.array(_xc)
	foccbin['yc']= np.array(_yc)
	foccbin['dlnx']= np.array(_dlnx)
	foccbin['dlny']= np.array(_dlny)

	foccbin['n']= np.array(_n)
	foccbin['comp']= np.array(_comp)
	foccbin['occ_MLE']= np.array(_occ)
	foccbin['err_MLE']= np.where(foccbin['n']>0,
		foccbin['occ_MLE']/np.sqrt(foccbin['n']),
		1./epos.nstars/ foccbin['comp'])

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

def _model_occ_per_polygon(epos, foccpoly, foccmodel, weights=None):	
	_occ, _n, _inbin, _xc, _yc, _dlnxy=[],[],[],[], [], []
	
	xyp= np.hstack([epos.pfm['P'][:,None],epos.pfm['R'][:,None]])

	for xycoords in foccpoly['coords']:
		shp=shapely.geometry.polygon.Polygon(xycoords)
		shp_log=shapely.geometry.polygon.Polygon(np.log(xycoords))

		points= map(shapely.geometry.Point, xyp)
		inpoly= map(shp.contains, points)

		logcenter= shp_log.centroid
		
		_inbin.append(inpoly)
		_n.append(np.sum(inpoly))
		_occ.append(weights* np.sum(inpoly))
		_dlnxy.append(np.log(shp_log.area))

		_xc.append(np.exp(logcenter.x))
		_yc.append(np.exp(logcenter.y))

		print '  poly, n={}, occ={:.2g}, area={:.2g}'.format(_n[-1], _occ[-1], _dlnxy[-1])
		print '    xc= {:.2g}, yc={:.2g}'.format(_xc[-1], _yc[-1])
		
	_foccpoly= foccmodel['poly']= {}
	_foccpoly['n']= np.array(_n)
	_foccpoly['i']= np.array(_inbin)
	_foccpoly['occ']= np.array(_occ)
	_foccpoly['err']= _foccpoly['occ']/np.where(
							_foccpoly['n']>0,np.sqrt(_foccpoly['n']),1.)
	
	_foccpoly['xc']= np.array(_xc)
	_foccpoly['yc']= np.array(_yc)
	_foccpoly['dlnxy']= np.array(_dlnxy)
	
	_foccpoly['coords']= foccpoly['coords']

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

def _interpolate_occ(epos, out_xvar, out_yvar, Log=False):
	if Log:
		# does not seem to make a big difference
		pl_comp_log= interpolate.RectBivariateSpline(np.log10(epos.eff_xvar), 
			np.log10(epos.eff_yvar), np.log10(epos.completeness))
		out_completeness= 10.**pl_comp_log(np.log10(out_xvar), np.log10(out_yvar), grid=False) 
		#print completeness
	else:
		# if zeros in detection efficiency?
		pl_comp= interpolate.RectBivariateSpline(epos.eff_xvar, epos.eff_yvar, epos.completeness)
		out_completeness= pl_comp(out_xvar, out_yvar, grid=False)

	return out_completeness

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
	

def binned_metric(epos):
	''' occurrence for one y bin, along x axis'''
	if 'yzoom' in epos.occurrence:
		xc= epos.occurrence['yzoom']['xc']
		occ_obs= epos.occurrence['yzoom']['occ'] / epos.occurrence['yzoom']['dlnx']
		occ_err= epos.occurrence['yzoom']['err'] / epos.occurrence['yzoom']['dlnx']

		_, _, pdf_X, _= periodradius(epos, ybin=epos.yzoom, xgrid= epos.occurrence['yzoom']['xc'])
		#for args in zip(xc, occ_obs, occ_err, pdf_X):
			#print 'xc= {:.2g}, occ={:.2g}+-{:.2g}, model= {:.2g}'.format(*args)

		squared_residuals= ((occ_obs-pdf_X)/occ_err)**2.

		zoom= (epos.xzoom[0] < xc) & (xc < epos.xzoom[-1])

		rss= np.sum(squared_residuals[zoom])
		nbins= zoom.sum()
		kfree= epos.fitpars.get_kfree()
		dof= nbins-kfree # kfree includes normalization parameter
		chi2= rss / dof

		bic= EPOS.analytics.bic_rss(rss, kfree, nbins)
		aic= EPOS.analytics.aic_rss(rss, kfree, nbins)
		aic_c= EPOS.analytics.aic_c_rss(rss, kfree, nbins)

		print '  x: (n={}, k={})'.format(nbins, kfree)
		print '    chi^2= {:.1f}, reduced= {:.1f}'.format(rss, chi2)
		print '    bic= {:.1f}'.format(bic)
		print '    aic= {:.1f}, AICc= {:.1f}'.format(aic, aic_c)

	if 'xzoom' in epos.occurrence:
		yc= epos.occurrence['xzoom']['yc']
		occ_obs= epos.occurrence['xzoom']['occ'] / epos.occurrence['xzoom']['dlny']
		occ_err= epos.occurrence['xzoom']['err'] / epos.occurrence['xzoom']['dlny']

		_, _, _, pdf_Y= periodradius(epos, xbin=epos.xzoom, ygrid= epos.occurrence['xzoom']['yc'])
		#for args in zip(yc, occ_obs, occ_err, pdf_Y):
			#print 'yc= {:.2g}, occ={:.2g}+-{:.2g}, model= {:.2g}'.format(*args)

		squared_residuals= ((occ_obs-pdf_Y)/occ_err)**2.

		zoom= (epos.yzoom[0] < yc) & (yc < epos.yzoom[-1])

		rss= np.sum(squared_residuals[zoom])
		nbins= zoom.sum()
		kfree= epos.fitpars.get_kfree()
		dof= nbins-kfree # kfree includes normalization parameter
		chi2= rss / dof

		bic= EPOS.analytics.bic_rss(rss, kfree, nbins)
		aic= EPOS.analytics.aic_rss(rss, kfree, nbins)
		aic_c= EPOS.analytics.aic_c_rss(rss, kfree, nbins)

		print '  y: (n={}, k={})'.format(nbins, kfree)
		print '    chi^2= {:.1f}, reduced= {:.1f}'.format(rss, chi2)
		print '    bic= {:.1f}'.format(bic)
		print '    aic= {:.1f}, AICc= {:.1f}'.format(aic, aic_c)
