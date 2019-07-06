import numpy as np
from scipy import interpolate
from scipy.stats import ks_2samp, anderson_ksamp, norm, chi2_contingency, kstest
import os, sys, logging, time
from functools import partial

import cgs
import multi
import EPOS.analytics
from EPOS.fitfunctions import brokenpowerlaw1D
from EPOS.population import periodradius

try:
	import emcee
except ImportError:
	print '#WARNING# emcee could not be imported'

def once(epos, fac=1.0, Extra=None, goftype='KS'):
	'''
	Run EPOS once
	
	Description:
		Sets up some defaults and runs the MC simulation once (steps 1-4 in paper)
		
	Args:
		Extra: store the planet population as an extra for plotting, default None
	'''
	epos.goftype=goftype
	
	if not epos.Prep:
		
		print '\nPreparing EPOS run...'
		# prep the input population (this should be after MC run?)
		if epos.Parametric:
			fpar2d= epos.pdfpars.get2d(Init=True)
			print '  {} fit parameters'.format(len(fpar2d))
		
			# Call function once to see if it works, P=1, R=1
			try: 	epos.func(np.asarray([1.]), np.asarray([1.]), *fpar2d)
			except:	raise
			
			''' Check if all parameters are set, set defaults '''
			
			if epos.Multi:
				epos.pdfpars.default('f_cor',0.5)
				epos.pdfpars.default('f_iso',0.5)
				epos.pdfpars.default('inc',2)
				#epos.pdfpars.default('',)
				if epos.RandomPairing:
					epos.pdfpars.default('npl',5)
				else:
					epos.pdfpars.default('dR',0.1)
					if epos.spacing == 'brokenpowerlaw':
						epos.pdfpars.default('dP break',1.7)
						epos.pdfpars.default('dP 1',10)
						epos.pdfpars.default('dP 2',-3)
					elif epos.spacing == 'dimensionless':
						epos.pdfpars.default('log D',-0.3)
						epos.pdfpars.default('sigma',0.2)
					else:
						raise ValueError('no spacing defined')

		else:
			# set defaults for planet formation models here
			epos.fitpars.default('eta',0.5)
			epos.fitpars.default('f_cor',0.5)
			epos.fitpars.default('f_iso',0.5)
			epos.fitpars.default('f_inc',1.0)
			epos.fitpars.default('f_dP',1.0)

			if not epos.MassRadius and not 'R' in epos.pfm:
				raise ValueError('Supply radii or a mass-radius relation')
		
		# prep the detection efficiency / observations
		if not epos.Range: epos.set_ranges()
		prep_obs(epos) # make pdf, cdf
		epos.Prep=True
	
	''' set weights / parameters '''
	fpara= epos.fitpars.getfit(Init=True)
	
	''' Time the first MC run'''
	runtype= 'MC' if epos.MonteCarlo else 'noMC'
	if Extra is None:
		print '\nStarting the first {} run'.format(runtype)
	else:
		print '\nStarting extra {} run {}'.format(runtype, Extra)
	tstart=time.time()
	runonce= MC if epos.MonteCarlo else noMC
	runonce(epos, fpara, Store=True, Extra=Extra)
	tMC= time.time()
	print 'Finished one {} in {:.3f} sec'.format(runtype, tMC-tstart)
	epos.tMC= tMC-tstart
	
def mcmc(epos, nMC=500, nwalkers=100, dx=0.1, nburn=50, threads=1, npos=30, Saved=True):
	if not 'emcee' in sys.modules:
		raise ImportError('You need to install emcee')
	assert epos.Prep
	
	runonce= MC if epos.MonteCarlo else noMC
	
	''' set starting parameters '''
	fpara= epos.fitpars.getfit(Init=True)
		
	if not len(fpara)>0: raise ValueError('no fit paramaters defined')
	
	''' Load previous chain?'''
	ndim= len(fpara)
	shape= (nwalkers, nMC, ndim)
	
	# store, npy is uncompressed, savez returns as dict instead of array
	dir= 'chain/{}'.format(epos.name)
	fname= '{}/{}x{}x{}.npz'.format(dir, nwalkers, nMC, ndim)
	if not os.path.exists(dir): os.makedirs(dir)
	
	if os.path.isfile(fname) and Saved:
		print '\nLoading saved status from {}'.format(fname)
		npz= np.load(fname)
		
		epos.chain=npz['chain']
		assert epos.chain.shape == (nwalkers, nMC, ndim)
		
		if epos.seed!=npz['seed']: 
			print '\nNOTE: Random seed changed: {} to {}'.format(npz['seed'],epos.seed)
		epos.seed= npz['seed']
		
		# check if same keys (not very flexible)
		if 'fitkeys' in npz:
			for loadkey,key in zip(npz['fitkeys'],epos.fitpars.keysfit):
				if loadkey != key:
					raise ValueError('Stored key {} doesnt match {}'.format(loadkey,key))
	else:
		
		''' start the timer '''
		tstart=time.time()
		nsims= nMC*nwalkers
		runtime= (epos.tMC/3600.)*nsims # single-threaded run time
		print '\nPredicted runtime:'
		if runtime>1:
			print '  {:.3f} hours for {} runs at {:.3f} sec'.format(
					runtime, nsims, epos.tMC)
		else:
			print '  {:.1f} minutes for {} runs at {:.3f} sec'.format(
					runtime*60., nsims, epos.tMC)
		if threads > 1:
			if runtime/threads>1:
				print '  {:.3f} hours at 100% scaling'.format(runtime/threads)
			else:
				print '  {:.1f} minutes at 100% scaling'.format(runtime/threads*60.)
	
		''' Set up a log file '''
		# Note: if logging.debug is called before this, the log file is never created
		fdir= 'log/{}'.format(epos.name)
		if not os.path.exists(fdir): os.makedirs(fdir)
		flog='{}/{}.log'.format(fdir,time.strftime('%Y-%m-%d_%H.%M.%S'))
		logging.basicConfig(filename=flog,level=logging.DEBUG,filemode='w')
		logging.info('MCMC with random seed {}'.format(epos.seed))
		
		''' Wrap function '''
		lnmc= partial(runonce, epos, Verbose=False)
	
		''' Set up the MCMC walkers '''
		#p0 = [np.array(fpara)*np.random.uniform(1.-dx,1+dx,len(fpara)) 
		#		for i in range(nwalkers)]
		dx=np.array(epos.fitpars.getfit(attr='dx'))
		p0 = [np.array(fpara)+dx*np.random.uniform(-1,1,len(fpara)) 
				for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, len(fpara), lnmc, threads=threads)
	
		''' run the chain '''
		if True:
			# chop to pieces for progress bar?
			for i, result in enumerate(sampler.sample(p0, iterations=nMC)):
				amtDone= float(i)/nMC
				print '\r  [{:50s}] {:5.1f}%'.format('#' * int(amtDone * 50), amtDone * 100),
				os.sys.stdout.flush() 
		else:
			sampler.run_mcmc(p0, nMC)
		
		print '\nDone running\n'
		logging.info('Made it to the end')
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
		#np.save(fname, epos.chain)
		# compression slow on loading?
		np.savez_compressed(fname, chain=epos.chain, seed=epos.seed, keys=epos.fitpars.keysfit)
		
	''' the posterior samples after burn-in '''
	epos.samples= epos.chain[:, nburn:, :].reshape((-1, ndim))
	epos.burnin= nburn
	fitpars = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(epos.samples, [16, 50, 84],
                                                axis=0)))
	epos.fitpars.setfit([p[0] for p in fitpars])
	
	''' Generate posterior populations '''
	if npos is not None:
		epos.plotsample= epos.samples[np.random.randint(len(epos.samples), size=npos)]
		# run & store
		print '\nMC-ing the {} samples to plot'.format(npos)
		epos.ss_sample=[]
		for fpara in epos.plotsample:
			epos.ss_sample.append(\
				runonce(epos, fpara, Store=True, Sample=True, Verbose=False))
		# parallel
		# return & store in structure
	
	''' Estimate Solar System Analogs'''
	if epos.Parametric and epos.Multi:
		fMercury=[]
		fVenus=[]
		for sample in epos.samples:
			_, _, xpdf, _= periodradius(epos, fpara=sample)
			fMercury.append(np.sum(xpdf[epos.MC_xvar>88.])/np.sum(xpdf))
			fVenus.append(np.sum(xpdf[epos.MC_xvar>225.])/np.sum(xpdf))
		
		print
		for name, posterior in zip(['Mercury','Venus'],[fMercury, fVenus]):
			eta= np.percentile(posterior, [16, 50, 84]) 
			print '{} analogues < {:.1%} +{:.1%} -{:.1%}'.format(name, eta[1], 
					eta[2]-eta[1], eta[1]-eta[0])
			UL= np.percentile(posterior, [68.2, 95.4, 99.7])
			for i in range(3): print '  {} sigma UL {:.1%}'.format(i+1,UL[i])


	''' Best-fit values'''
	print '\nBest-fit values'
	for pname, fpar in zip(epos.fitpars.keysfit, fitpars): 
		print '  {}= {:.3g} +{:.3g} -{:.3g}'.format(pname,*fpar)

	print '\nStarting the best-fit MC run'	
	runonce(epos, np.array([p[0] for p in fitpars]), Store=True, BestFit=True)
	
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
	z['multi']['pl cnt']= z['multi']['bin']*z['multi']['count']
	
	z['multi']['Pratio'], z['multi']['Pinner']= \
		multi.periodratio(epos.obs_starID[ix&iy], epos.obs_xvar[ix&iy])
	z['multi']['cdf']= multi.cdf(epos.obs_starID[ix&iy])

def MC(epos, fpara, Store=False, Sample=False, StorePopulation=False, Extra=None, 
		BestFit=False, Verbose=True):
	''' 
	Do the Monte Carlo Simulations
	Note:
	variable x/X is P
	variable y/Y is R/M
	TODO: split into multiple functions
	'''	
	if Verbose: tstart=time.time()
	#if not Store: logging.debug(' '.join(['{:.3g}'.format(fpar) for fpar in fpara]))
	
	''' Seed the random number generator '''
	if epos.seed is not None: np.random.seed(epos.seed)
	
	''' construct 1D arrays allP, allR or allM
	dimension equal to sample size * planets_per_star
	also keeping track of:
		ID: star identifier
		I: inc
		N: Nth planet in system
		dP: period ratio
	'''
	if epos.Parametric:
		
		''' parameters within bounds? '''
		try:
			epos.pdfpars.checkbounds(fpara)
		except ValueError as message:
			if Store: raise
			else:
				logging.debug(message)
				return -np.inf
				
		''' Draw (inner) planet from distribution '''
		pps= epos.fitpars.getpps_fromlist(fpara)
		fpar2d= epos.fitpars.get2d_fromlist(fpara)
		npl= epos.fitpars.getmc('npl', fpara) if epos.RandomPairing else 1
		
		try: sysX, sysY= draw_from_2D_distribution(epos, pps, fpar2d, npl=npl)
		except ValueError:
			if Store: raise
			else: return -np.inf
		
		''' Multi-planet systems '''
		if not epos.Multi:
			allX= sysX
			allY= sysY
			
		elif epos.RandomPairing:			
			# set ID, nth planet in system
			isys= np.arange(sysX.size/npl)
			allID= np.repeat(isys,npl)
			order= np.lexsort((sysX,allID)) # sort by ID, then P
			assert np.all(allID[order]==allID)

			allX= sysX[order]
			allY= sysY #[order] # random order anyway
			allN= np.tile(np.arange(npl),isys.size) # ignores xzoom
			
			# isotropic or inc?
			dInc= epos.fitpars.getmc('inc', fpara)
			if dInc is not None:
				f_iso= epos.fitpars.getmc('f_iso', fpara)
				allI= np.random.rayleigh(dInc, allID.size)
			
			f_cor= epos.fitpars.getmc('f_cor', fpara)
			f_dP, f_inc= 1.0, 1.0 # no need to fudge these
			
		else:
			''' retrieve fit parameters for multi-planets'''
			npl= epos.fitpars.getmc('npl', fpara)
			dR= epos.fitpars.getmc('dR', fpara)
			dInc= epos.fitpars.getmc('inc', fpara)
			f_iso= epos.fitpars.getmc('f_iso', fpara)
			f_cor= epos.fitpars.getmc('f_cor', fpara)
			
			f_dP, f_inc= 1.0, 1.0 # no need to fudge these

			
			''' Parameter bounds '''
			if npl < 1:
				logging.debug('npl = {:.3g} < 1'.format(npl))
				if Store: raise ValueError('at least one planet per system')
				return -np.inf
			if (dInc <=0) or (dR <=0) or not (0 <= f_iso <= 1):
				if Store: raise ValueError('parameters out of bounds')
				return -np.inf
			
			''' Draw multiplanet distributions '''
			allX, allY, allI, allN, allID= \
				draw_multi(epos, sysX, sysY, npl, dInc, dR, fpara, Store)
					
		''' convert to observable parameters '''
		allP= allX
		if epos.RV or epos.MassRadius: allM= allY
		else:		allR= allY
				
	else:
		pfm= epos.pfm
			
		''' parameters within bounds? '''
		# move out of loop?
		try:
			epos.fitpars.checkbounds(fpara)
		except ValueError as message:
			if Store: raise
			else:
				logging.debug(message)
				return -np.inf
				
		''' Fit parameters'''
		pps= epos.fitpars.getpps_fromlist(fpara)

		f_cor= epos.fitpars.getmc('f_cor', fpara)
		f_iso= epos.fitpars.getmc('f_iso', fpara)
		f_inc= epos.fitpars.getmc('f_inc', fpara)
		f_dP= epos.fitpars.getmc('f_dP', fpara)
					
		# need this here?
		if not (0<=f_iso<=1) or not (0 < pps) or not (0 <= f_cor <= 1):
			#\or not (0<=f_dP<=10) or not (0 <= f_inc < 10):
			if Store: raise ValueError('parameters out of bounds')
			return -np.inf
		
		''' 
		Draw from all
		'''
		if not 'draw prob' in pfm:
			ndraw= int(round(1.*epos.nstars*pps/pfm['ns']))
			if Verbose: 
				print '  {} planets in {} simulations'.format(pfm['np'],pfm['ns'])
				print '  {} stars in survey, {} draws, eta={:.2g}'.format(epos.nstars, ndraw, pps)
		
			allP= np.tile(pfm['P'], ndraw)
			allM= np.tile(pfm['M'], ndraw)
			if 'R' in pfm:
				allR= np.tile(pfm['R'], ndraw)
		
			if epos.Multi:
				allI= np.tile(pfm['inc'], ndraw)
				allN= np.tile(pfm['kth'], ndraw)
				
				# ID or system index?
				allID= np.tile(pfm['ID'], ndraw) \
						+ np.repeat(np.arange(ndraw)*pfm['ns'], pfm['np'])
		else:
			'''
			Draw from some distributions according to 'tag' parameter
			TODO: functions to calculate draw probability from tag
			'''
			#draw planetary systems from simulations
			ndraw= int(round(1.*epos.nstars*pps))
			if Verbose: print '\nDraw {} systems'.format(ndraw) 
			system_index= np.random.choice(pfm['system index'], size=ndraw, 
							p=pfm['draw prob'])
			
			#create a list of planets
			indexlist= [pfm['planet index'][pfm['ID']==i] for i in system_index]
			starIDlist= [np.full((pfm['ID']==i).sum(),ID) 
							for ID,i in enumerate(system_index)]
			
			planets= np.concatenate(indexlist)
			allID= np.concatenate(starIDlist)
			allP= pfm['P'][planets]
			allM= pfm['M'][planets]
			allR= pfm['R'][planets]
			allI= pfm['inc'][planets]
			allN= pfm['kth'][planets]
			#allID= pfm['ID'][planets]
			
			if Verbose: print '  {} planets'.format(allP.size)
		
		allY= allM
			
		dInc=False # Isotropic inclinations not implemented
			
	''' 
	Identify transiting planets (itrans is a T/F array)
	'''
	if epos.RV:
		# RV keep all
		itrans= np.full(allP.size, True, np.bool)
	elif (not epos.Multi) or (epos.Multi and dInc==None):
		# geometric transit probability
		p_trans= epos.fgeo_prefac *allP**epos.Pindex
		itrans= p_trans >= np.random.uniform(0,1,allP.size)
	else:
		#multi-transit probability
		itrans= istransit(epos, allID, allI, allP, f_iso, f_inc, Verbose=Verbose)
		
	# Print multi statistics	
	if Verbose and epos.Multi and not epos.RV: 
		print '\n  {} planets, {} transit their star'.format(itrans.size, itrans.sum())
		multi.frequency(allID[itrans], Verbose=True)
	
	'''
	remove planets according to transit probability
	'''	
	MC_P= allP[itrans]
	if epos.MassRadius or epos.RV:	MC_M= allM[itrans]
	else:							MC_R= allR[itrans]	
	if epos.Multi:
		MC_ID= allID[itrans]	
		if epos.Parametric:
			MC_N= allN[itrans] # also for PFM?
		else:
			#MC_P*= (1.+0.1*np.random.normal(size=MC_ID.size) )
			pass

	'''
	Set the observable MC_Y (R or Msin i) 
	'''
	if epos.RV:
		''' M sin i.'''
		# Note different conventions for i in Msini (i=0 is pole-on)
		# sin(arccos(chi)) == cos(arcsin(chi)) == sqrt(1-chi^2)
		MC_Y= MC_Msini= MC_M* np.sqrt(1.-np.random.uniform(0,1,MC_P.size)**2.)
	else:
		''' Convert Mass to Radius '''
		if epos.MassRadius:
			mean, dispersion= epos.MR(MC_M)
			MC_R= mean+ dispersion*np.random.normal(size=MC_M.size)
		
		''' uncertainty in stellar radius? '''
		MC_Y=MC_R * (1.+epos.radiusError*np.random.normal(size=MC_R.size) )

	
	'''
	Store (transiting) planet sample for verification plot
	'''
	if Store:
		tr= epos.transit={}
		tr['P']= MC_P
		tr['Y']= MC_Y
	
	'''
	Identify detectable planets based on SNR (idet is a T/F array)
	'''
	f_snr= interpolate.RectBivariateSpline(epos.MC_xvar, epos.MC_yvar, epos.MC_eff)
	p_snr= f_snr(MC_P, MC_Y, grid=False)
	assert p_snr.ndim == 1

	idet= p_snr >= np.random.uniform(0,1,MC_P.size)
	
	# draw same random number for S/N calc, 1=correlated noise
	if epos.Multi and f_cor >0:
		IDsys, toplanet= np.unique(MC_ID, return_inverse=True) 
		idet_cor= p_snr >= np.random.uniform(0,1,IDsys.size)[toplanet]

		cor_sys= (np.random.uniform(0,1,IDsys.size) < f_cor)
		cor_pl= cor_sys[toplanet]
		idet = np.where(cor_pl, idet_cor, idet)

	''' 
	Remove undetectable planets
	'''
	# arrays with detected planets
	if epos.Multi:
		det_ID= MC_ID[idet]
		if not epos.RV and epos.Parametric: det_N= MC_N[idet]
	det_P= MC_P[idet]
	det_Y= MC_Y[idet]

	#if len(alldP)>0: 
	# 	if not epos.Parametric and 'all_Pratio' in sg: 
	# 		det_dP= MC_dP[idet]

	if Verbose and epos.Multi: 		
		print '  {} transiting planets, {} detectable'.format(idet.size, idet.sum())
		multi.frequency(det_ID, Verbose=True)

	'''
	Probability that simulated data matches observables
	TODO:
	- store ln prob, D
	- switches for in/excluding parameters
	- likelyhood from D instead of prob?
	'''
	# bugfix 
	#np.seterr(divide='raise')

	tstart=time.time()
	prob={}
	lnp={}
	
	# make sure that x=P, y=R (where?)
	ix= (epos.xzoom[0]<=det_P) & (det_P<=epos.xzoom[1])
	iy= (epos.yzoom[0]<=det_Y) & (det_Y<=epos.yzoom[1])
	if (ix&iy).sum()<1:
		if Store: raise ValueError('no planets detectable')
		return -np.inf			
	
	if epos.goftype=='KS':
		prob_2samp= _prob_ks
	elif epos.goftype=='AD': 
		prob_2samp= _prob_ad
	else:
		raise ValueError('{} not a goodness-of-fit type (KS, AD)'.format(epos.goftype))

	if 'xvar' in epos.summarystatistic:
		prob['xvar'], lnp['xvar']=  prob_2samp(epos.obs_zoom['x'], det_P[ix&iy])
	if 'yvar' in epos.summarystatistic:
		prob['yvar'], lnp['yvar']=  prob_2samp(epos.obs_zoom['y'], det_Y[ix&iy])

	if 'N' in epos.summarystatistic:
		# chi^2: (np-nobs)/nobs**0.5 -> p: e^-0.5 x^2
		chi2= (epos.obs_zoom['x'].size-np.sum(ix&iy))**2. / epos.obs_zoom['x'].size
		lnp['N']= -0.5* chi2
		prob['N']= np.exp(-0.5* chi2)
	
	if epos.Multi:
		''' Multi-planet frequency, pearson chi_squared '''
		k, Nk= multi.frequency(det_ID[ix&iy])
		Nk_obs= epos.obs_zoom['multi']['count']
		ncont= max(len(Nk),len(Nk_obs))
		
		# pad with zeros
		obs= np.zeros((2,ncont), dtype=int)
		obs[0,:len(Nk)]= Nk
		obs[1,:len(Nk_obs)]= Nk_obs
		# remove double zero frequencies
		obs=obs[:, ~((obs[0,:] == 0) & (obs[1,:]==0))]
		
		try:	
			_, prob['Nk'], _, _ = chi2_contingency(obs)
		except ValueError:
			#print Nk
			#print Nk_obs
			#print obs
			raise
			
		with np.errstate(divide='ignore'): lnp['Nk']= np.log(prob['Nk'])			
		
		''' Period ratio, innermost planet '''
		sim_dP, sim_Pinner= multi.periodratio(det_ID[ix&iy], det_P[ix&iy])

		if (len(sim_dP)>0) & (len(sim_Pinner)>0): 
			prob['dP'],lnp['dP']= prob_2samp(epos.obs_zoom['multi']['Pratio'],f_dP*sim_dP)					
			prob['Pin'],lnp['Pin']= prob_2samp(epos.obs_zoom['multi']['Pinner'],
											sim_Pinner)
		else:
			logging.debug('no multi-planet statistics, {}'.format(len(sim_dP)))
			prob['dP'], prob['Pin']= 0, 0 
			lnp['dP'], lnp['Pin']= -np.inf, -np.inf

	''' combine with Fischer's rule: '''
	lnprob= np.sum([lnp[key] for key in epos.summarystatistic])
	#dof= len(prob_keys)
	chi_fischer= -2. * lnprob

	if Verbose:
		print '\nGoodness-of-fit'
		print '  logp= {:.1f}'.format(lnprob)
		print '  - p(n={})={:.2g}'.format(np.sum(ix&iy), prob['N'])
		if 'xvar' in prob:	print '  - p(x)={:.2g}'.format(prob['xvar'])
		if 'yvar' in prob:	print '  - p(y)={:.2g}'.format(prob['yvar'])
		if 'Nk' in prob:	print '  - p(N_k)={:.2g}'.format(prob['Nk'])
		if 'dP' in prob:	print '  - p(P ratio)={:.2g}'.format(prob['dP'])
		if 'Pin' in prob:	print '  - p(P inner)={:.2g}'.format(prob['Pin'])

		if BestFit: 
			''' Akaike/Bayesian information criterion'''
			k_free= epos.fitpars.get_kfree()
			n_data= epos.obs_zoom['x'].size
			bic= EPOS.analytics.bic_loglike(lnprob, k_free, n_data)
			aic= EPOS.analytics.aic_loglike(lnprob, k_free)
			aic_c= EPOS.analytics.aic_c_loglike(lnprob, k_free, n_data)

			print '\n  Akaike/Bayesian Information Criterion'
			print '  - k={}, n={}'.format(k_free,n_data)
			print '  - BIC= {:.1f}'.format(bic)
			print '  - AIC= {:.1f}, AICc= {:.1f}'.format(aic,aic_c)
		
	tgof= time.time()
	if Verbose: print '  observation comparison in {:.3f} sec'.format(tgof-tstart)

	''' Store _systems_ with at least one detected planet '''
	# StorePopulation
	if Store and epos.Multi and (not epos.RV):
		# include/exclude singles w/ iso_pl?
		itransdet= np.copy(itrans)
		itransdet[itrans]= idet
		isysdet,isingle,imulti,order= storepopulation(allID, allP, det_ID, itransdet)
				
		pop= epos.population={}
		pop['order']= order
		pop['P']= allP
		pop['k']= allN
		pop['inc']= allI
		for key, subset in zip(['system', 'single', 'multi'],[isysdet, isingle, imulti]):
			pop[key]={}
			pop[key]['Y']= allY[subset] # R? M?
			pop[key]['P']= allP[subset]
			pop[key]['ID']= allID[subset]
			pop[key]['inc']= allI[subset]
			pop[key]['order']= order[subset]
			pop[key]['detectable']= itransdet[subset]
	
	''' Store detectable planet population '''	
	if Store:
		ss={}
		ss['P']= det_P # MC_P[idet]
		ss['Y']= det_Y
		if epos.MassRadius:
			# or if has mass and radius
			ss['M']= MC_M[idet]
			ss['R']= MC_R[idet]
		
		#if len(alldP)>0
		if epos.Multi:
			ss['ID']= det_ID
			ss['multi']={}
			ss['multi']['bin'], ss['multi']['count']= multi.frequency(det_ID[ix&iy])
			ss['multi']['pl cnt']=ss['multi']['bin']* ss['multi']['count']
			ss['multi']['Pratio'], ss['multi']['Pinner']= \
				multi.periodratio(det_ID[ix&iy], det_P[ix&iy]) # *f_dP
			ss['multi']['cdf']= multi.cdf(det_ID[ix&iy])
			if not epos.RV and epos.Parametric:
				ss['multi']['PN'],ss['multi']['dPN']= multi.periodratio(
						det_ID[ix&iy], det_P[ix&iy], N=det_N[ix&iy]) 
		
		epos.prob=prob
		epos.lnprob=lnprob
		ss['P zoom']= det_P[ix&iy]
		ss['Y zoom']= det_Y[ix&iy]
		ss['nobs']= ss['P zoom'].size
		
		# Store as an extra model 
		if Sample:
			return ss
		elif Extra is not None:
			ss['name']=Extra
			if not hasattr(epos,'ss_extra'):
				epos.ss_extra=[]
			epos.ss_extra.append(ss)
			#print 'saving extra {}'.format(Extra)
		else:
			epos.synthetic_survey= ss 
	else:
		# return probability
		if np.isnan(lnprob):
			return -np.inf
		return lnprob

def noMC(epos, fpara, Store=False, Sample=False, StorePopulation=False, Extra=None, 
		BestFit=False, Verbose=True):
	''' 
	Do the Simulations without Monte Carlo
	'''	
	if Verbose: tstart=time.time()
	#if not Store: logging.debug(' '.join(['{:.3g}'.format(fpar) for fpar in fpara]))

	if epos.Multi: raise ValueError('Multi-planets need Monte Carlo (?)')
	if not epos.Parametric: 
		raise ValueError('Planet Formation models need Monte Carlo (?)')
		
	''' parameters within bounds? '''
	try:
		epos.pdfpars.checkbounds(fpara)
	except ValueError as message:
		if Store: raise
		else:
			logging.debug(message)
			return -np.inf

	''' Generate observable period-radius distribution, in counts'''
	if epos.RV:
		pps, pdf, pdf_X, pdf_Y= periodradius(epos, fpara=fpara, fdet=epos.MC_eff*epos.nstars, 
			Convert=epos.Msini)
	elif epos.MassRadius:
		raise ValueError('Generate pdf on radius grid here')
	else:
		pps, pdf, pdf_X, pdf_Y= periodradius(epos, fpara=fpara, fdet=epos.f_det*epos.nstars)

	'''
	Probability that simulated data matches observables
	'''
	tstart=time.time()
	prob={}
	lnp={}
	
	pdf_x= np.interp(epos.noMC_zoom_x, epos.MC_xvar, pdf_X)
	pdf_y= np.interp(epos.noMC_zoom_y, epos.MC_yvar, pdf_Y)
	cdf_x= np.cumsum(pdf_x)
	cdf_y= np.cumsum(pdf_y)
	cdf_x/= cdf_x[-1]
	cdf_y/= cdf_y[-1] 
	
	''' Wrap functions '''
	func_cdf_X= partial(np.interp, xp=epos.noMC_zoom_x, fp=cdf_x, left=0, right=0)
	func_cdf_Y= partial(np.interp, xp=epos.noMC_zoom_y, fp=cdf_y, left=0, right=0)
		 	
	prob['xvar'], lnp['xvar']= _prob_ks_func(epos.obs_zoom['x'], func_cdf_X)
	prob['yvar'], lnp['yvar']= _prob_ks_func(epos.obs_zoom['y'], func_cdf_Y)

	#nobs= int(pdf[epos.trim_to_zoom].sum())
	nobs_x=int(np.sum(pdf_x)/epos.noMC_scale_x)
	nobs_y=int(np.sum(pdf_y)/epos.noMC_scale_y)
	nobs= int(np.sqrt(nobs_x*nobs_y))
	if Verbose: print 'nobs={} (x:{},y:{})'.format(nobs,nobs_x,nobs_y) # ??

	chi2= (epos.obs_zoom['x'].size-nobs)**2. / epos.obs_zoom['x'].size
	lnp['N']= -0.5* chi2
	prob['N']= np.exp(-0.5* chi2)
		
	prob_keys= ['N','xvar','yvar']

	# combine with Fischer's rule:
	lnprob= np.sum([lnp[key] for key in prob_keys])
	#dof= len(prob_keys)
	chi_fischer= -2. * lnprob
	
	
	if Verbose:
		print '\nGoodness-of-fit'
		print '  logp= {:.1f}'.format(lnprob)
		print '  - p(n={})={:.2g}'.format(nobs, prob['N'])
		print '  - p(x)={:.2g}'.format(prob['xvar'])
		print '  - p(y)={:.2g}'.format(prob['yvar'])

		if BestFit:
			''' Bayesian information criterion'''
			k_free= epos.fitpars.get_kfree()
			n_data= epos.obs_zoom['x'].size
			bic= EPOS.analytics.bic_loglike(lnprob, k_free, n_data)
			aic= EPOS.analytics.aic_loglike(lnprob, k_free)
			aic_c= EPOS.analytics.aic_c_loglike(lnprob, k_free, n_data)

			print '\n  Akaike/Bayesian Information Criterion'
			print '  - k={}, n={}'.format(k_free,n_data)
			print '  - BIC= {:.1f}'.format(bic)
			print '  - AIC= {:.1f}, AICc= {:.1f}'.format(aic,aic_c)
		
	tgof= time.time()
	if Verbose: print '\nObservation comparison in {:.3f} sec'.format(tgof-tstart)

	''' Store detectable planet population '''	
	if Store:
		ss={}
		ss['nobs']= nobs
		
		ss['pdf']= pdf

		ss['P']= epos.MC_xvar
		ss['P pdf']= pdf_X
		
		ss['Y']= epos.MC_yvar
		ss['Y pdf']= pdf_Y
		
		ss['P zoom']= epos.noMC_zoom_x
		ss['P zoom pdf']= pdf_x # /scale?
		ss['P zoom cdf']= cdf_x

		ss['Y zoom']= epos.noMC_zoom_y
		ss['Y zoom pdf']= pdf_y # /scale?
		ss['Y zoom cdf']= cdf_y
		
		epos.prob=prob
		epos.lnprob=lnprob

		if Sample:
			return ss
		elif Extra is not None:
			ss['name']=Extra
			if not hasattr(epos,'ss_extra'):
				epos.ss_extra=[]
			epos.ss_extra.append(ss)
			#print 'saving extra {}'.format(Extra)
		else:
			epos.synthetic_survey= ss
	else:
		''' return probability '''
		if np.isnan(lnprob):
			return -np.inf
		return lnprob
	
def draw_from_2D_distribution(epos, pps, fpara, npl=1):
	
	''' create PDF, CDF'''
	# assumes a separable function of mass and radius
	pdf= epos.func(epos.X_in, epos.Y_in, *fpara)

	pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
	cum_X, cum_Y= np.cumsum(pdf_X), np.cumsum(pdf_Y)
	#pps_x, pps_y=  cum_X[-1], cum_Y[-1]
	#planets_per_star= 0.5*(pps_x+pps_y) # should be equal
	
	try:
		ndraw= npl*int(round(pps*epos.nstars))
	except OverflowError:
		logging.debug('infinity encountered, pps= {},{}'.format(pps_x, pps_y))
		for pname, fpar in zip(epos.pname, fpara):
			logging.debug('  {}= {:.3g}'.format(pname,fpar))
		raise ValueError('Infinity encountered')
	
	if ndraw < 1: 
		logging.debug('no draws ({}, {}*{})'.format(ndraw, epos.nstars, planets_per_star))
		raise ValueError('no planets')
	elif ndraw > 1e8:
		logging.debug('>1e8 planets ({})'.format(ndraw))
		raise ValueError('too many planets')
	# 	elif planets_per_star > 100:
	# 		logging.debug('>100 planets per star ({})'.format(planets_per_star))
	# 		raise ValueError('too many planets per star')
	try:
		allX= np.interp(np.random.uniform(cum_X[0],cum_X[-1],ndraw), cum_X, epos.MC_xvar)
		allY= np.interp(np.random.uniform(cum_Y[0],cum_Y[-1],ndraw), cum_Y, epos.in_yvar)
	except MemoryError:
		raise ValueError('Memory error for n={}'.format(ndraw))
	except OverflowError:
		print ndraw
		print cum_X
		print cum_Y
		raise
	
	return allX, allY

def draw_multi(epos, sysX, sysY, npl, dInc, dR, fpara, Store):
	''' assign ID to each system '''
	sysID= np.arange(sysX.size)
	#allID= np.repeat(sysID, npl) # array with star ID for each planet
	npl_arr= np.random.uniform(npl, npl+1, sysID.size) # rounds down
	allID= np.repeat(sysID, npl_arr.astype(int))
	#print allID.size, sysID.size
	if allID.size > 1e7:
		logging.debug('Too many planets: {} > 1e7'.format(allID.size))
		for pname, fpar in zip(epos.pname, fpara):
			logging.debug('  {}= {:.3g}'.format(pname,fpar))
		if Store: raise ValueError('Too many planets: {}'.format(allID.size))
		return -np.inf				

	''' initialize planet parameters'''
	_, toplanet, sysnpl= np.unique(allID, return_inverse=True,return_counts=True)
	allX= sysX[toplanet]
	allY= sysY[toplanet]
	allI= np.random.rayleigh(dInc, allID.size)
	#allN= np.ones_like(allID) # index to planet in system
	allN= np.where(allX>=epos.xzoom[0],1,0) # also yzoom?
	#print allX[:3]

	# get index of first planet in each system
	if len(sysnpl) < 1:
		logger.debug('no planets')
		if Store: raise ValueError('no planets')
		return -np.inf
	di= np.roll(sysnpl,1)
	di[0]=0
	i1= np.cumsum(di) # index to first planet
	assert np.all(sysID == allID[i1])
	
	''' Draw period of 2nd, 3rd planet etc.'''
	Pgrid= np.logspace(0,1)
	# use population.periodratio here
	if epos.spacing == 'powerlaw':
		dPbreak= epos.fitpars.getmc('dP break', fpara)
		dP1= epos.fitpars.getmc('dP 1', fpara)
		dP2= epos.fitpars.getmc('dP 2', fpara)
		if (dPbreak<=0):
			if Store: raise ValueError('parameters out of bounds')
			return -np.inf
		cdf= np.cumsum(brokenpowerlaw1D(Pgrid, dPbreak, dP1, dP2))
	elif epos.spacing=='dimensionless':
		logD=  epos.fitpars.getmc('log D', fpara)
		sigma= epos.fitpars.getmc('sigma', fpara)
		if (sigma<=0):
			if Store: raise ValueError('parameters out of bounds')
			return -np.inf			

		with np.errstate(divide='ignore'): 
			Dgrid= np.log10(2.*(Pgrid**(2./3.)-1.)/(Pgrid**(2./3.)+1.))
		Dgrid[0]= -2
		#print Dgrid
		cdf= norm(logD,sigma).cdf(Dgrid)
	
	# loop over planet 2,3... n
	for i in range(2,len(np.bincount(sysnpl)) ):
		im= i1[sysnpl>=i] # index to ith planet in each system
		#print '  planet {} in system, {} systems'.format(i,im.size)
		#dP= draw_from_function(brokenpowerlaw1D,xgrid,im.size, dPbreak, dP1, dP2)
		dP= np.interp(np.random.uniform(cdf[0],cdf[-1],im.size), cdf, Pgrid)
		allX[im+(i-1)]= allX[im+(i-2)]* dP
		#np.random.uniform(dP-0.7, dP+0.7, im.size)
		#print allX[:3]
		#np.random.norm()
		allY[im+(i-1)]= allY[im+(i-2)]* 10.**np.random.normal(0, dR, im.size)
		
		# nth planet in system (zoom range)
		#allN[im+(i-1)]= allN[im+(i-2)]+1
		#allN[im+(i-1)]= allN[im+(i-2)]+np.where(allX[im+(i-1)]>=epos.xzoom[0],1,0)
		allN[im+(i-1)]= np.where(allX[im+(i-1)]>=epos.xzoom[0],allN[im+(i-2)]+1,0)

		
	#print '+{}={} story checks out?'.format(allID[i1].size, allID.size)
	#print allX[:3]
	
	''' Toss out planets (reduces memory footprint) '''
	Xinrange= allX<=epos.MC_xvar[-1]
	#if Xinrange.sum() < 0.5* allX.size: 
	#	print 'yeah, {}/{}'.format(Xinrange.sum(),allX.size)
	allX= allX[Xinrange]
	allY= allY[Xinrange]
	allI= allI[Xinrange]
	allN= allN[Xinrange]
	allID= allID[Xinrange]
	
	return allX, allY, allI, allN, allID
	
def istransit(epos, allID, allI, allP, f_iso, f_inc, Verbose=False):
	# draw same numbers for multi-planet systems
	IDsys, toplanet= np.unique(allID, return_inverse=True) # return_counts=True
	if Verbose: print '  {}/{} systems'.format(IDsys.size, allID.size)
	
	# draw system viewing angle proportionate to sin theta (i=0: edge-on)
	inc_sys= np.arcsin(np.random.uniform(0,1,IDsys.size))
	inc_pl= inc_sys[toplanet]
	assert inc_pl.size == allP.size
	
	R_a= epos.fgeo_prefac *allP**epos.Pindex # == p_trans
	mutual_inc= allI * f_inc
	#mutual_inc= 0.0 # planar distribution
	#mutual_inc= 1.0 # fit 
	if Verbose: 
		print '  Average mutual inc={:.1f} degrees'.format(np.median(allI))
		if f_inc != 1.0:
			print 'f_inc= {:.2g}, inc= {:.1f} deg'.format(f_inc, np.median(mutual_inc))
	delta_inc= mutual_inc *np.cos(np.random.uniform(0,np.pi,allP.size)) *np.pi/180. #rad
	itrans= np.abs(inc_pl+delta_inc) < np.arcsin(R_a)

	# allow for a fraction of isotropic systems
	if f_iso > 0:
		p_trans= epos.fgeo_prefac *allP**epos.Pindex
		itrans_iso= p_trans >= np.random.uniform(0,1,allP.size)
		
		iso_sys= (np.random.uniform(0,1,IDsys.size) < f_iso)
		iso_pl= iso_sys[toplanet]
		itrans = np.where(iso_pl, itrans_iso, itrans)
		
	return itrans

def storepopulation(allID, allP, det_ID, idetected):
	# add f_iso?

	# systems with at least one detectable planet
	try:
		isysdet= np.isin(allID, det_ID)
	except AttributeError:
		isysdet= np.in1d(allID, det_ID)
	
	# detected multis/singles, index to allID
	uniqueID, detcounts= np.unique(allID[idetected],return_counts=True)
	IDmulti= uniqueID[(detcounts>1)] # systems that are multi
	IDsingle= uniqueID[(detcounts==1)]
	# returns all planets in system (not just detectable ones)
	try:
		imulti=  np.isin(allID, IDmulti)
		isingle= np.isin(allID, IDsingle)
	except AttributeError:
		imulti=  np.in1d(allID, IDmulti)
		isingle= np.in1d(allID, IDsingle)

	# Sort order [0-1] based on inner period (also of undetected)
	IDsys, iunique, toplanet= np.unique(allID, return_index=True, return_inverse=True)
	assert np.all(IDsys[toplanet] == allID)
	Pinner= allP[iunique] # inner if lexsorted
	
	idx= np.argsort(Pinner) #[::-1]
	rank= np.argsort(idx)
	
	#Porder= rank[toplanet]
	Porder= 1.0 * rank[toplanet] / IDsys.size
 
	return isysdet, isingle&idetected, imulti&idetected, Porder

def _prob_ks_func(a,func):
	_, prob= kstest(a,func)
	with np.errstate(divide='ignore'):
		lnprob= np.log(prob)
	return prob, lnprob

def _prob_ks(a,b):
	_, prob= ks_2samp(a,b)
	with np.errstate(divide='ignore'):
		lnprob= np.log(prob)
	return prob, lnprob

def _prob_ad(a,b):
	_, _, prob= anderson_ksamp([a,b])
	with np.errstate(divide='ignore'):
		lnprob= np.log(prob)
	if prob > 1:
		print
		print anderson_ksamp([a,b])
		print ks_2samp(a,b)
		print
		print a
		print b
		print prob, lnprob
	return prob, lnprob

''' Old code '''
# def draw_from_function(f, grid, ndraw, *args):
# 	cdf= np.cumsum(f(grid, *args))
# 	return np.interp(np.random.uniform(0,cdf[-1],ndraw), cdf, grid)

#def make_pdf(epos, norm=None, Init=False):