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
		prep_obs(epos) # make pdf, cdf
		epos.Prep=True
	
	''' set weights / parameters '''
	if epos.populationtype is 'parametric':
		fpara= epos.p0 
	elif epos.populationtype is 'model':
		fpara= [sg['weight'] for sg in epos.groups]
		#if hasattr(epos,'p0'): fpara=epos.p0
	else: assert False
	
	''' Time the first MC run'''
	print '\nStarting the first MC run'
	tstart=time.time()
	MC(epos, fpara, Store=True)
	tMC= time.time()
	print 'Finished one MC in {:.3f} sec'.format(tMC-tstart)
	epos.tMC= tMC-tstart
	
def mcmc(epos, nMC=500, nwalkers=100, dx=0.1, nburn=50, threads=1):
	assert epos.Prep
	import emcee
	# multi-threading: not much of a performance increase past 4 threads
	
	''' set starting parameters '''
	if epos.populationtype is 'parametric':
		fpara= epos.p0
	elif epos.populationtype is 'model':
		if len(epos.groups) == 1:
			#fpara= [sg['weight'], 0.5, 0.8, 1.0] # singles, Pratio, inc (+corr SNR?)
			fpara=epos.p0 # defined where?
		else:
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
		runtime= (epos.tMC/3600.)*nsims # single-threaded run time
		if runtime>1:
			print '\nPredicted runtime {:.3f} hours for {} runs at {:.3f} sec'.format(
					runtime, nsims, epos.tMC)
		else:
			print '\nPredicted runtime {:.1f} minutes for {} runs at {:.3f} sec'.format(
					runtime*60., nsims, epos.tMC)
	
		''' Wrap function '''
		lnmc= partial(MC, epos, Verbose=False, LogProb= True)
	
		''' Set up the MCMC walkers '''
		p0 = [np.array(fpara)*np.random.uniform(1.-dx,1+dx,len(fpara)) 
				for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, len(fpara), lnmc, threads=threads)
	
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
	if epos.populationtype is 'parametric':
		pps= [np.sum(epos.func(epos.X, epos.Y,*para)) for para in epos.samples]
		eta= np.percentile(pps, [16, 50, 84]) 
		print '  eta= {:.3g} +{:.3g} -{:.3g}'.format(eta[1], 
				eta[2]-eta[1], eta[1]-eta[0])

		# planets in box (where is the X,Y grid stored?)
		#pps= [np.sum(epos.func(epos.X, epos.Y,*para)) for para in epos.samples]
		#eta= np.percentile(pps, [16, 50, 84]) 
		#print '  eta_box= {:.3g} +{:.3g} -{:.3g}'.format(eta[1], 
		#		eta[2]-eta[1], eta[1]-eta[0])

	
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
	MC(epos, Store=True, fpara=epos.pfit)
	
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
	z['multi']['cdf']= multi.cdf(epos.obs_starID[ix&iy])

def MC(epos, fpara, Store=False, Verbose=True, KS=True, LogProb=False):
	''' 
	Do the Monte Carlo Simulations
	Note:
	variable x/X is P
	variable y/Y is R/M
	'''	
	
	tstart=time.time()
	
	''' construct 1D arrays for P, R
	dimension equal to sample size * planets_per_star
	also keeping track off subgroup (SG), inc (I), and period ratio (dP)
	'''
	if epos.populationtype == 'parametric':
		''' Check if parameters make sense'''
		# normalization can't be zero
		if fpara[0] <= 0:
			#print 'oops: {}'.format(fpara) # rare
			if Store: raise ValueError('normalization needs to be larger than 0')
			return -np.inf
		
		# hack for double power law
		if len(fpara) is 7:
			if (fpara[1]<= 0) or (fpara[4]<= 0):
				#print 'oops II: {}'.format(fpara) # rare
				if Store: raise ValueError('breaks needs to be larger than 0')
				return -np.inf
		
		''' create PDF, CDF'''
		# assumes a seperable function of mass and radius
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
		allX= np.interp(np.random.uniform(0,pps_x,ndraw), cum_X, epos.MC_xvar)
		allY= np.interp(np.random.uniform(0,pps_y,ndraw), cum_Y, epos.MC_yvar)
		
		''' convert to observable parameters '''
		allP= allX
		if epos.RV:	allM= allY
		else:		allR= allY
				
	else:
		''' Fit parameters '''
		f_SNR= 0.5 # doesn't seem to be well-constrained, light degenracy with f_inc
		if len(fpara) is len(epos.groups):
			''' Draw systems from subgroups, array with length survey size '''
			weight= fpara
			f_iso= 0.0 # 0.6 is best fit?
			f_dP= 1.0
			f_inc= 1.0
			#f_SNR= 0.5 
			if hasattr(epos,'p0'): 
				weight= [epos.p0[0]]
				#f_iso, f_dP, f_inc, f_SNR= epos.p0[1:]
				f_iso, f_dP, f_inc= epos.p0[1:]
		else:
			try:
				weight= [fpara[0]]
				f_iso= fpara[1]
				f_dP= fpara[2]
				f_inc= fpara[3]
				#f_SNR= fpara[4]
			except:
				raise ValueError('\Wrong number of parameters ({})'.format(len(fpara)))
			
			if f_iso<0 or f_dP<=0 or f_inc<0 or not (0 <= weight[0] <= 1) or \
				not (0 <= f_SNR <= 1):
				#print 'Oops IV'
				#print fpara
				if Store: raise ValueError('parameters out of bounds')
				return -np.inf
				
		
		L_P, L_M, L_R, L_I= [], [], [], []
		L_SG=[] # keep track of subgroup
		L_ID=[] # keep track of multis
		L_dP=[] # Period ratios (input)
		IDstart=0
		for i, sg in enumerate(epos.groups):
			ndraw= int(round(1.*epos.nstars*weight[i]/sg['n']))
			# note: do a proper calc with // and % for low weights
			if Verbose: 
				print '  {} planets in subgroup {} with {} simulations'.format(
					sg['all_P'].size, sg['name'], sg['n'])
				print '  {} stars in survey, {} draws from subgroup (w={:.3f})'.format(
					epos.nstars, ndraw, weight[i])
			L_P.append(np.tile(sg['all_P'], ndraw))
			L_M.append(np.tile(sg['all_mass'], ndraw))
			L_SG.append(np.full_like(L_M[-1],i))
			for j in range(ndraw):
				L_ID.append(sg['all_ID']+ IDstart)
				IDstart= L_ID[-1][-1]
			
			if 'all_Pratio' in sg: L_dP.append(np.tile(sg['all_Pratio'], ndraw))
			if 'all_R' in sg: L_R.append(np.tile(sg['all_R'], ndraw))

			if not epos.Isotropic: L_I.append(np.tile(sg['all_inc'], ndraw))


		allP=np.concatenate(L_P)
		allM=np.concatenate(L_M)
		allSG=np.concatenate(L_SG)
		allID=np.concatenate(L_ID)
		if not epos.Isotropic: allI=np.concatenate(L_I)
		if 'all_Pratio' in sg: alldP=np.concatenate(L_dP)
		if 'all_R' in sg: allR=np.concatenate(L_R)
		#print len(np.unique(L_allID))/ndraw # == sg['n']
		
	''' 
	remove planets according to geometric transit probability 
	'''
	if epos.RV:
		#itrans= np.full(allP.size, True, np.bool) # keep all
		MC_P= allP
		MC_M= allM
	elif epos.Isotropic:
		p_trans= epos.fgeo_prefac *allP**epos.Pindex
		itrans= p_trans >= np.random.uniform(0,1,allP.size)

		MC_P= allP[itrans] 
		MC_R= allR[itrans]	
	else:
		# draw same numbers for multi-planet systems
		IDsys, toplanet= np.unique(allID, return_inverse=True) # return_counts=True
		if Verbose: print '  {}/{} systems'.format(IDsys.size, allID.size)
		
		# draw system viewing angle proportionate to sin theta
		inc_sys= np.arcsin(np.random.uniform(0,1,IDsys.size))
		inc_pl= inc_sys[toplanet]
		assert inc_pl.size == allP.size
		
		R_a= epos.fgeo_prefac *allP**epos.Pindex # == p_trans
		mutual_inc= allI * f_inc
		#mutual_inc= 0.0 # planar distribution
		#mutual_inc= 1.0 # fit 
		if Verbose: print '  Average mutual inc={:.1f} degrees'.format(np.median(allI))
		delta_inc= mutual_inc *np.cos(np.random.uniform(0,np.pi,allP.size)) * np.pi/180.
		itrans= np.abs(inc_pl+delta_inc) < np.arcsin(R_a)

		# allow for a fraction of isotropic systems
		if f_iso > 0:
			p_trans= epos.fgeo_prefac *allP**epos.Pindex
			itrans_iso= p_trans >= np.random.uniform(0,1,allP.size)
			
			iso_sys= (np.random.uniform(0,1,IDsys.size) < f_iso)
			iso_pl= iso_sys[toplanet]
			itrans = np.where(iso_pl, itrans_iso, itrans)

		MC_M= allM[itrans]	
		MC_P= allP[itrans] * (1.+0.1*np.random.normal(size=MC_M.size) )
		MC_SG= allSG[itrans]
		MC_ID= allID[itrans]
#		if len(alldP)>0: MC_dP= alldP[itrans]	
		if 'all_Pratio' in sg: MC_dP= alldP[itrans]
		
	# Print multi statistics	
	if Verbose and not (epos.Isotropic or epos.RV): 
		print '\n  {} planets, {} transit their star'.format(itrans.size, itrans.sum())
		multi.frequency(allID[itrans], Verbose=True)

	'''
	Set the observable MC_Y (R or Msin i) 
	'''
	if epos.RV:
		''' M sin i '''
		MC_Y= MC_Msini= MC_M* np.random.uniform(0,1,MC_P.size) # sin of arcsin?
	else:
		''' Convert Mass to Radius '''
		if epos.RadiusMassConversion:
			mean, dispersion= epos.RM(MC_M)
			MC_R= mean+ dispersion*np.random.normal(size=MC_M.size)
		MC_Y=MC_R
	
	''' Store transiting planet sample for verification plot'''
	if Store:
		tr= epos.transit={}

		tr['P']= MC_P
		tr['Y']= MC_Y
		if epos.populationtype is 'model':
			tr['i sg']= MC_SG
	
	'''
	remove planets due to SNR
	'''
	f_snr= interpolate.RectBivariateSpline(epos.MC_xvar, epos.MC_yvar, epos.MC_eff)
	p_snr= f_snr(MC_P, MC_Y, grid=False)
	assert p_snr.ndim == 1

	idet= p_snr >= np.random.uniform(0,1,MC_P.size)
	
	# draw same random number for S/N calc, 1=correlated noise
	if (not epos.Isotropic) and f_SNR >0:
		IDsys, toplanet= np.unique(MC_ID, return_inverse=True) 
		idet_cor= p_snr >= np.random.uniform(0,1,IDsys.size)[toplanet]

		cor_sys= (np.random.uniform(0,1,IDsys.size) < f_SNR)
		cor_pl= cor_sys[toplanet]
		idet = np.where(cor_pl, idet_cor, idet)

	# arrays with detected planets
	if not epos.Isotropic:
		det_ID= MC_ID[idet]		
	det_P= MC_P[idet]
	det_Y= MC_Y[idet]


	#if len(alldP)>0: 
	if epos.populationtype is 'model' and 'all_Pratio' in sg: 
		det_dP= MC_dP[idet]
	if Verbose and not epos.Isotropic: 
		print '  {} transiting planets, {} detectable'.format(idet.size, idet.sum())
		multi.frequency(det_ID, Verbose=True)


	'''
	Probability that simulated data matches observables
	'''
	if KS:
		tstart=time.time()
		gof={}
		# make sure that x=P, y=R (where?)
		ix= (epos.xzoom[0]<=det_P) & (det_P<=epos.xzoom[1])
		iy= (epos.yzoom[0]<=det_Y) & (det_Y<=epos.yzoom[1])
	
		gof['D xvar'], gof['p xvar']=  ks_2samp(epos.obs_zoom['x'], det_P[ix&iy])
		gof['D yvar'], gof['p yvar']=  ks_2samp(epos.obs_zoom['y'], det_Y[ix&iy])
		# chi^2: (np-nobs)/nobs**0.5 -> p: e^-0.5 x^2
		gof['p n']= np.exp(-0.5*(epos.obs_zoom['x'].size-np.sum(ix&iy))**2. 
							/ epos.obs_zoom['x'].size )
		prob= 1.* gof['p n']
		if epos.Isotropic:
			prob*= gof['p xvar']* gof['p yvar'] # increase acceptance fraction for models
		
		if not epos.Isotropic:
			# Period ratio
			gof['D dP'], gof['p dP']= ks_2samp(epos.obs_zoom['multi']['Pratio'], 
						f_dP*multi.periodratio(det_ID[ix&iy], det_P[ix&iy]))
			prob*= gof['p dP']
			
			# Multi-planets
			multis= multi.cdf(det_ID[ix&iy], Verbose=False)
			gof['D multi'], gof['p multi']=  \
				ks_2samp(epos.obs_zoom['multi']['cdf'], multis)
			prob*= gof['p multi']
		
		if Verbose:
			print '\nGoodness-of-fit'
			print '  logp= {:.1f}'.format(np.log(prob))
			print '  - p(x)={:.2g}'.format(gof['p xvar'])
			print '  - p(y)={:.2g}'.format(gof['p yvar'])
			print '  - p(n={})={:.2g}'.format(np.sum(ix&iy), gof['p n'])
			if not epos.Isotropic:
				print '  - p(multi)={:.2g}'.format(gof['p multi'])
				print '  - p(P ratio)={:.2g}'.format(gof['p dP'])
			

		tgof= time.time()
		if Verbose: print '  observation comparison in {:.3f} sec'.format(tgof-tstart)


	if Store:	
		ss= epos.synthetic_survey= {}
		ss['P']= det_P # MC_P[idet]
		ss['Y']= det_Y
		if not epos.Isotropic:
			ss['i sg']= MC_SG[idet]
		if epos.RadiusMassConversion:
			# or if has mass and radius
			ss['M']= MC_M[idet]
			ss['R']= MC_R[idet]
		
		#if len(alldP)>0
		if not epos.Isotropic:
			if 'all_Pratio' in sg: ss['dP']= det_dP
			ss['multi']={}
			ss['multi']['bin'], ss['multi']['count']= multi.frequency(det_ID[ix&iy])
			ss['multi']['Pratio']= multi.periodratio(det_ID[ix&iy], det_P[ix&iy]) # *f_dP
			ss['multi']['cdf']= multi.cdf(det_ID[ix&iy])
		
		epos.gof=gof
		ss['P zoom']= det_P[ix&iy]
		ss['Y zoom']= det_Y[ix&iy]
	else:
		# return probability
		lnprob= np.log(prob)
		if np.isnan(lnprob):
			print 'oops III'
			print gof['p xvar'], gof['p yvar'], gof['p n']
			print allP.size, MC_P.size, det_P.size
			return -np.inf
		return lnprob if LogProb else prob

		

	
