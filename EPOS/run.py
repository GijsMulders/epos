import numpy as np
import time
from scipy import interpolate
from scipy.stats import ks_2samp
import os
import cgs
import multi
from functools import partial

def once(epos, fac=1.0):
	'''
	A test run with equal weights
	TODO: throw in a bunch of assertions
	'''
	if not epos.Prep:
		
		print '\nPreparing EPOS run...'
		# prep the input population (this should be after MC run?)
		if epos.populationtype is 'parametric':
			#epos.p= epos.p0 # ??
			pass
		elif epos.populationtype is 'model':
			summedweight= np.sum([sg['weight'] for sg in epos.groups])
			for sg in epos.groups:
				sg['weight']*= fac/summedweight
				print 'set weight {} to {}'.format(sg['name'],sg['weight']) 
			#epos.weights= [fac/len(epos.groups)]* len(epos.groups) # equal weights
			assert epos.RadiusMassConversion, 'set mass-to-radius function'
		else: assert False
		
		# prep the detection efficiency / observations
		if not epos.Range: epos.set_ranges()
		prep_eff(epos) # set the MC grid
		prep_obs(epos) # make pdf, cdf
		epos.Prep=True
	
	''' set weights / parameters '''
	if epos.populationtype is 'parametric':
		fpara= epos.p0 
	elif epos.populationtype is 'model':
		fpara= [sg['weight'] for sg in epos.groups]
	else: assert False
	
	''' Time the first MC run'''
	print '\nStarting the first MC run'
	tstart=time.time()
	MC(epos, fpara, Store=True, Parametric=(epos.populationtype is 'parametric')) # parallel: epos.weights separetely
	tMC= time.time()
	print 'Finished one MC in {:.3f} sec'.format(tMC-tstart)
	epos.tMC= tMC-tstart
	
def mcmc(epos, nMC=500, nwalkers=100, dx=0.1, nburn=50):
	assert epos.Prep
	import emcee
	
	''' set starting parameters '''
	if epos.populationtype is 'parametric':
		fpara= epos.p0
	elif epos.populationtype is 'model':
		fpara= [sg['weight'] for sg in epos.groups]
	else: assert False
	
	''' Load previous chain?'''
	ndim= len(epos.p0)
	shape= (nwalkers, nMC, ndim)
	
	# store, npy is uncompressed, savez returns as dict instead of array
	dir= 'chain/{}'.format(epos.name)
	fname= '{}/{}x{}x{}.npy'.format(dir, nwalkers, nMC, ndim)
	if not os.path.exists(dir): os.makedirs(dir)
	if os.path.isfile(fname):
		print '\nLoading saved status from {}'.format(fname)
		epos.chain= np.load(fname) 
		assert epos.chain.shape == (nwalkers, nMC, ndim)
	else:
		
		''' start the timer '''
		tstart=time.time()
		nsims= nMC*nwalkers
		runtime= (epos.tMC/3600.)*nsims
		if runtime>1:
			print '\nPredicted runtime {:.3f} hours for {} runs at {:.3f} sec'.format(
					runtime, nsims, epos.tMC)
		else:
			print '\nPredicted runtime {:.1f} minutes for {} runs at {:.3f} sec'.format(
					runtime*60., nsims, epos.tMC)
	
		''' Wrap function '''
		lnmc= partial(MC, epos, Parametric=(epos.populationtype is 'parametric'),
						Verbose=False, LogProb= True)
	
		''' Set up the MCMC walkers '''
		p0 = [np.array(fpara)*np.random.uniform(1.-dx,1+dx,len(fpara)) 
				for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, len(fpara), lnmc) # args= ...
	
		''' run the chain '''
		# add progress bar?
		sampler.run_mcmc(p0, nMC)
		print '\nDone running\n'
		print 'Mean acceptance fraction: {0:.3f}'.format(
					np.mean(sampler.acceptance_fraction))

		''' Print run time'''	
		tMC= time.time()
		runtime= tMC-tstart
		if runtime>3600:
			print '  Runtime was {:.3f} hours at {:.3f} sec'.format(
					runtime/3600, (tMC-tstart)/nsims )
		else:
			print '  Runtime was {:.1f} minutes at {:.3f} sec'.format(
					runtime/60., (tMC-tstart)/nsims)
	
		epos.chain= sampler.chain
		print 'Saving status in {}'.format(fname)
		np.save(fname, epos.chain) 
		
	''' the posterior samples after burn-in '''
	epos.samples= epos.chain[:, nburn:, :].reshape((-1, len(epos.p0)))
	epos.burnin= nburn
	fitpars = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(epos.samples, [16, 50, 84],
                                                axis=0)))
	epos.pfit=[p[0] for p in fitpars]
	
	''' estimate #planets/star '''
	#print epos.samples.shape # 30000, 7
	pps= [np.sum(epos.func(epos.X, epos.Y,*para)) for para in epos.samples]
	eta= np.percentile(pps, [16, 50, 84]) 
	print '  eta= {:.3g} +{:.3g} -{:.3g}'.format(eta[1], 
			eta[2]-eta[1], eta[1]-eta[0])
	
	''' estimate eta_earth'''
	# eta_earth
	n_earth= epos.func(365., 1.0, *epos.samples.transpose() )
	xnorm= epos.MC_xvar.size-1 #  == xrange/dx
	ynorm= epos.MC_yvar.size-1 #  == yrange/dy
	xHZ= np.log10((1.67/0.95)**(1.5) )/ np.log10(epos.MC_xvar[-1]/epos.MC_xvar[0])
	yHZ= np.log10((1.5/0.7))/ np.log10(epos.MC_yvar[-1]/epos.MC_yvar[0])

	eta_earth= np.percentile(n_earth, [16, 50, 84]) * (xHZ*xnorm) * (yHZ*ynorm)
	print '  eta_earth= {:.3g} +{:.3g} -{:.3g}'.format(eta_earth[1], 
			eta_earth[2]-eta_earth[1], eta_earth[1]-eta_earth[0])

	''' Best-fit values'''
	print '\nBest-fit values'
	for pname, fpar in zip(epos.pname, fitpars): 
		print '  {}= {:.3g} +{:.3g} -{:.3g}'.format(pname,*fpar)


	print '\nStarting the best-fit MC run'	
	MC(epos, Store=True, Parametric=(epos.populationtype is 'parametric'),
		fpara=epos.pfit)
	
	
def prep_eff(epos):
	#epos.MC_scale= (epos.MC_xvar[-1]/epos.MC_xvar[-2])* (epos.MC_yvar[1]/epos.MC_yvar[0])
	#print 'scale={}'.format(epos.MC_scale)
	
	# MC det eff
	# NOTE: the efficiency factor epos.MC_max can be used to decrease sample size
	epos.Pindex= 2./3.
	epos.MC_eff= epos.eff_trim * epos.X**epos.Pindex
	epos.MC_max= np.max(epos.MC_eff)
	epos.MC_eff/= epos.MC_max
	
	# check transit probability at edge of grid
	# TODO: split transit probability out of detection efficiency for multis
# 	Rstar=1.1
# 	epos.prob_P= Rstar*cgs.Rsun/ (cgs.au /365.24**(2./3.))
# 	print 'P={}, prob={}, R/a= {}'.format(epos.MC_xvar[0],eff_trim[0,-1],
# 						epos.prob_P/epos.MC_xvar[0]**(2./3.))

def prep_obs(epos):
	# occurrence pdf on sma from plot_input_diag?

	# cdf of obs, zoomed only
	z=epos.obs_zoom={}
	ix= (epos.xzoom[0]<=epos.obs_xvar) & (epos.obs_xvar<=epos.xzoom[1])
	iy= (epos.yzoom[0]<=epos.obs_yvar) & (epos.obs_yvar<=epos.yzoom[1])
	
	if np.sum(ix&iy)==0:
		raise ValueError('No planets in the box!')
	
	x= epos.obs_xvar[ix&iy]
	y= epos.obs_yvar[ix&iy]

	for key, var in zip(['x','y'], [x,y]):
		z[key]= np.sort(var) # add (x0,0) and (x1,1)?
		z[key+' cum']= np.arange(z[key].size,dtype=float)
		z[key+' cdf']= z[key+' cum']/z[key+' cum'][-1]
	
	# multis
	z['multi']={}
	z['multi']['bin'], z['multi']['count']= multi.frequency(epos.obs_starID[ix&iy])
	
	z['multi']['Pratio']= multi.periodratio(epos.obs_starID[ix&iy], 
							epos.obs_xvar[ix&iy])

def MC(epos, fpara, Store=False, Verbose=True, Parametric=True, KS=True,
		CoPlanar=False, Isotropic=False, LogProb=False):
	''' construct arrays with P, M for MC 
	NOTE: increase speed by reducing sample size with epos.MC_max
		to discard non-transiting geometries a-priori
		(this doesn't seem necessary)
	'''	
	
	tstart=time.time()
	
	''' construct 1D arrays for P, R
	dimension equal to sample size
	also keeping track off subgroup (SG), inc (I), and period ratio (dP)
	'''
	if Parametric:
		
		''' Check if parameters make sense'''
		assert len(fpara) is 7, '{}'.format(fpara)
		if fpara[0] <= 0:
			#print 'oops: {}'.format(fpara) # rare
			if Store: raise ValueError('normalization needs to be larger than 0')
			return -np.inf
		if (fpara[1]<= 0) or (fpara[4]<= 0):
			#print 'oops II: {}'.format(fpara) # rare
			if Store: raise ValueError('breaks needs to be larger than 0')
			return -np.inf
		
		''' create PDF, CDF'''
		pdf= epos.func(epos.X, epos.Y, *fpara)
		pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
		cum_X, cum_Y= np.cumsum(pdf_X), np.cumsum(pdf_Y)
		pps_x, pps_y=  cum_X[-1], cum_Y[-1]
		planets_per_star= 0.5*(pps_x+pps_y) # should be equal
		tcdf= time.time()
		#print '  CDF in {:.3f} sec'.format(tcdf-tstart) # this is pretty fast :)
		
		ndraw= int(round(planets_per_star*epos.nstars))
		if ndraw < 1: 
			print 'no draws ({}, {}*{})'.format(ndraw, epos.nstars, planets_per_star)
			if Store: raise ValueError('no planets')
			return -np.inf
		allP= np.interp(np.random.uniform(0,pps_x,ndraw), cum_X, epos.MC_xvar)
		allR= np.interp(np.random.uniform(0,pps_y,ndraw), cum_Y, epos.MC_yvar)
		
# 		i,j= -1, 4
# 		print ' eta_earth= {:.2f}'.format(epos.pdf[i,j] *(4./epos.MC_scale) )
# 		print ' P={:.0f}, R={:.1f}'.format(epos.MC_xvar[i], epos.MC_xvar[j])
		
# 		if Store:
# 			epos.pdf= pdf
# 			epos.pdf_X, epos.pdf_Y= pdf_X, pdf_Y
# 			epos.planets_per_star= planets_per_star
		
	else:
		# Draw systems from subgroups, array with length survey size
		assert len(fpara) is len(epos.groups)
		L_P, L_M, L_I= [], [], []
		L_SG=[] # keep track of subgroup
		L_ID=[] # keep track of multis
		L_dP=[] # Period ratios (input)
		IDstart=0
		for i, sg in enumerate(epos.groups):
			ndraw= int(round(1.*epos.nstars*sg['weight']/sg['n']))
			# note: do a proper calc with // and % for low weights
			if Verbose: 
				print '  {} planets in subgroup {} with {} simulations'.format(sg['all_P'].size, sg['name'], sg['n'])
				print '  {} stars in survey, {} draws from subgroup (w={:.3f})'.format(epos.nstars, ndraw, sg['weight'])
			L_P.append(np.tile(sg['all_P'], ndraw))
			L_M.append(np.tile(sg['all_mass'], ndraw))
			L_SG.append(np.full_like(L_M[-1],i))
			for j in range(ndraw):
				L_ID.append(sg['all_ID']+ IDstart)
				IDstart= L_ID[-1][-1]
			
			if 'all_Pratio' in sg: L_dP.append(np.tile(sg['all_Pratio'], ndraw))
			if 'all_inc' in sg: L_I.append(np.tile(sg['all_inc'], ndraw))

		allP=np.concatenate(L_P)
		allM=np.concatenate(L_M)
		allSG=np.concatenate(L_SG)
		allID=np.concatenate(L_ID)
		if 'all_Pratio' in sg: alldP=np.concatenate(L_dP)
		if 'all_inc' in sg: allI=np.concatenate(L_I)
		#print len(np.unique(L_allID))/ndraw # == sg['n']
		
	# remove planets according to geometric transit probability
	# isotropic inclinations:
	#Isotropic=False
	Isotropic= True if Parametric else (not ('all_inc' in sg))
	if Isotropic:
		p_trans= epos.MC_max* allP**(-epos.Pindex)
		itrans= p_trans >= np.random.uniform(0,1,allP.size)
	elif False:
		# test transit probability calculation with same random numbers
		themc= np.random.uniform(0,1,allP.size)
		R_a= cgs.Rsun/ (cgs.au*(allP/365.24)**(2./3.))
		
		# from integrated probability
		itrans= R_a >= themc
		
		# from MC
		mc_trans= np.arcsin(themc) < np.arcsin(R_a)
		print itrans.sum(), mc_trans.sum() # duh
	else:
		# draw same numbers for multi-planet systems
		IDsys, toplanet= np.unique(allID, return_inverse=True) # return_counts=True
		print '  {}/{} systems'.format(IDsys.size, allID.size)
		
		# draw system viewing angle proportionate to sin theta
		inc_sys= np.arcsin(np.random.uniform(0,1,IDsys.size))
		inc_pl= inc_sys[toplanet]
		assert inc_pl.size == allP.size
		
		R_a= cgs.Rsun/ (cgs.au*(allP/365.24)**(2./3.))
		mutual_inc= allI
		#mutual_inc= 0.0 # planar distribution
		#mutual_inc= 1.0 # fit 
		print '  Average mutual inc={:.1f} degrees'.format(np.median(allI))
		delta_inc= mutual_inc *np.cos(np.random.uniform(0,np.pi,allP.size)) * np.pi/180.
		itrans= np.abs(inc_pl+delta_inc) < np.arcsin(R_a)
	
	# Print multi statistics	
	if Verbose and not Parametric: 
		print '\n  {} planets, {} transit their star'.format(itrans.size, itrans.sum())
		multi.frequency(allID[itrans], Verbose=True)

		
	# going to MC with these arrays
	if Parametric:
		MC_P= allP[itrans] 
		MC_R= allR[itrans]
	else:
		MC_M= allM[itrans]	
		MC_P= allP[itrans] * (1.+0.1*np.random.normal(size=MC_M.size) )
		MC_SG= allSG[itrans]
		MC_ID= allID[itrans]
#		if len(alldP)>0: MC_dP= alldP[itrans]	
		if 'all_Pratio' in sg: MC_dP= alldP[itrans]
	
		# convert mass to radius
		mean, dispersion= epos.RM(MC_M)
		MC_R= mean+ dispersion*np.random.normal(size=MC_M.size)	
		#MC_R= epos.RM_mean(MC_M)+ epos.RM_dispersion(MC_M)*np.random.normal(size=MC_M.size)

		
	''' other diagnostics? (multi)'''
	#if Multi
	
	if Store:
		tr= epos.transit={} # store for verification plot.
		tr['R']= MC_R
		tr['P']= MC_P
		if not Parametric:
			tr['M']= MC_M
			tr['i sg']= MC_SG
	
	# remove planets due to SNR ( self.MC_eff )
	f_snr= interpolate.RectBivariateSpline(epos.MC_xvar, epos.MC_yvar, epos.MC_eff)
	p_snr= f_snr(MC_P, MC_R, grid=False)
	assert p_snr.ndim == 1

	if True or Isotropic:
		idet= p_snr >= np.random.uniform(0,1,MC_P.size)
	else:
		# draw same random number for S/N calc
		# correlated noise
		# makes a big difference for smaller planets?
		IDsys, toplanet= np.unique(MC_ID, return_inverse=True) 
		idet= p_snr >= np.random.uniform(0,1,IDsys.size)[toplanet]
	
		det_ID= MC_ID[idet]
		
	det_P= MC_P[idet]
	det_R= MC_R[idet]
	#if len(alldP)>0: 
	if (not Parametric) and 'all_Pratio' in sg: det_dP= MC_dP[idet]
	#print type(idet), idet.shape
	#print np.min(p_snr), np.max(p_snr)
	if Verbose and not Parametric: 
		print '  {} transiting planets, {} detectable'.format(idet.size, idet.sum())
		multi.frequency(det_ID, Verbose=True)


	# Compare with obs
	# Note: KS doesn't have normalization (!)
	if KS:
		tstart=time.time()
		gof={}
		# make sure that x=P, y=R (where?)
		ix= (epos.xzoom[0]<=det_P) & (det_P<=epos.xzoom[1])
		iy= (epos.yzoom[0]<=det_R) & (det_R<=epos.yzoom[1])
	
		gof['D xvar'], gof['p xvar']=  ks_2samp(epos.obs_zoom['x'], det_P[ix&iy])
		gof['D yvar'], gof['p yvar']=  ks_2samp(epos.obs_zoom['y'], det_R[ix&iy])
		# chi^2: (np-nobs)/nobs**0.5 -> p: e^-0.5 x^2
		gof['p n']= np.exp(-0.5*(epos.obs_zoom['x'].size-np.sum(ix&iy))**2. /epos.obs_zoom['x'].size )
		
		if Verbose:
			print '  logp= {:.2f}'.format(np.log(gof['p xvar']* gof['p yvar']*gof['p n']))
			print '  - p(x)={:.2e}'.format(gof['p xvar'])
			print '  - p(y)={:.2e}'.format(gof['p yvar'])
			print '  - p(n={})={:.2e}'.format(np.sum(ix&iy), gof['p n'])
			#print gof['p n'], epos.obs_zoom['x'].size,  (ix&iy).size
			# add to list: probability that they are drawn from the same distribution (0=bad, 1=good)
		# if chain
		# if para
		tgof= time.time()
		if Verbose: print '  observation comparison in {:.3f} sec'.format(tgof-tstart)


	if Store:	
		ss= epos.synthetic_survey= {}
		ss['R']= det_R # MC_R[idet]
		ss['P']= det_P # MC_P[idet]
		if not Parametric:
			ss['M']= MC_M[idet]
			ss['i sg']= MC_SG[idet]
		
		#if len(alldP)>0
		if not Parametric:
			if 'all_Pratio' in sg: ss['dP']= det_dP
			ss['multi']={}
			ss['multi']['bin'], ss['multi']['count']= multi.frequency(det_ID[ix&iy])
			ss['multi']['Pratio']= multi.periodratio(det_ID[ix&iy], det_P[ix&iy])
		
		epos.gof=gof # not parallel proof
		ss['P zoom']= det_P[ix&iy]
		ss['R zoom']= det_R[ix&iy]
	else:
		# normalization
		prob= gof['p xvar']* gof['p yvar']*gof['p n']
		return np.log(prob) if LogProb else prob

		

	
