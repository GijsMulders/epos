import numpy as np
from scipy import interpolate
from scipy.stats import ks_2samp, norm
import os, sys, logging, time
from functools import partial

import cgs
import multi
from EPOS.fitfunctions import brokenpowerlaw1D

try:
	import emcee
except ImportError:
	print '#WARNING# emcee could not be imported'

def once(epos, fac=1.0, Extra=False):
	'''
	A test run with equal weights
	TODO: throw in a bunch of assertions
	'''
	if not epos.Prep:
		
		print '\nPreparing EPOS run...'
		# prep the input population (this should be after MC run?)
		if epos.populationtype is 'parametric':
			fpar2d= epos.fitpars.get2d(Init=True)
			print '  {} fit parameters'.format(len(fpar2d))
		
			# Call function once to see if it works, P=1, R=1
			try: 	epos.func(np.asarray([1.]), np.asarray([1.]), *fpar2d)
			except:	raise
			
			''' Check if all parameters are set, set defaults '''
			
			if epos.Multi:
				epos.fitpars.default('f_cor',0.5)
				epos.fitpars.default('f_iso',0.5)
				epos.fitpars.default('inc',2)
				#epos.fitpars.default('',)
				if epos.RandomPairing:
					epos.fitpars.default('npl',5)
				else:
					epos.fitpars.default('dR',0.1)
					if epos.spacing == 'brokenpowerlaw':
						epos.fitpars.default('dP break',1.7)
						epos.fitpars.default('dP 1',10)
						epos.fitpars.default('dP 2',-3)
					elif epos.spacing == 'dimensionless':
						epos.fitpars.default('log D',-0.3)
						epos.fitpars.default('sigma',0.2)
					else:
						raise ValueError('no spacing defined')

		elif epos.populationtype is 'model':
			summedweight= np.sum([sg['weight'] for sg in epos.groups])
			for sg in epos.groups:
				sg['weight']*= fac/summedweight
				print 'set weight {} to {}'.format(sg['name'],sg['weight']) 
			#epos.weights= [fac/len(epos.groups)]* len(epos.groups) # equal weights
			assert epos.MassRadius, 'set mass-to-radius function'
		else: assert False
		
		# prep the detection efficiency / observations
		if not epos.Range: epos.set_ranges()
		prep_obs(epos) # make pdf, cdf
		epos.Prep=True
	
	''' set weights / parameters '''
	if epos.populationtype is 'parametric':
		fpara= epos.fitpars.getfit(Init=True)
	elif epos.populationtype is 'model':
		fpara= [sg['weight'] for sg in epos.groups]
		if hasattr(epos,'p0'): fpara=epos.p0
	else: assert False
	
	''' Time the first MC run'''
	print '\nStarting the first MC run'
	tstart=time.time()
	MC(epos, fpara, Store=True, Extra=Extra)
	tMC= time.time()
	print 'Finished one MC in {:.3f} sec'.format(tMC-tstart)
	epos.tMC= tMC-tstart
	
def mcmc(epos, nMC=500, nwalkers=100, dx=0.1, nburn=50, threads=1):
	if not 'emcee' in sys.modules:
		raise ImportError('You need to install emcee')
	assert epos.Prep
	
	''' set starting parameters '''
	if epos.populationtype is 'parametric':
		fpara= epos.fitpars.getfit(Init=True)
	elif epos.populationtype is 'model':
		if len(epos.groups) == 1:
			fpara=epos.p0 # [epos.ipfit] ??
		else:
			fpara= [sg['weight'] for sg in epos.groups]
	else: assert False
	
	if not len(fpara)>0: raise ValueError('no fit paramaters defined')
	
	''' Load previous chain?'''
	ndim= len(fpara)
	shape= (nwalkers, nMC, ndim)
	
	# store, npy is uncompressed, savez returns as dict instead of array
	dir= 'chain/{}'.format(epos.name)
	#fname= '{}/{}x{}x{}.npy'.format(dir, nwalkers, nMC, ndim)
	fname= '{}/{}x{}x{}.npz'.format(dir, nwalkers, nMC, ndim)
	if not os.path.exists(dir): os.makedirs(dir)
	if os.path.isfile(fname):
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
		lnmc= partial(MC, epos, Verbose=False, LogProb= True)
	
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
	if epos.populationtype is 'parametric':
		epos.fitpars.setfit([p[0] for p in fitpars])
		#print 'before: {}'.format(epos.pfit)
		#print 'after: {}'.format(epos.pfit)
	else:
		# set weight instead?
		epos.pfit=[p[0] for p in fitpars]
	
	''' Estimate Solar System Analogs'''
	if epos.populationtype is 'parametric' and epos.Multi:
		fMercury=[]
		fVenus=[]
		for sample in epos.samples:
			_, _, xpdf, _= _pdf(epos, fpara=sample)
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
	MC(epos, np.array([p[0] for p in fitpars]), Store=True)
	
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

def MC(epos, fpara, Store=False, StorePopulation=False, Extra=False, 
		Verbose=True, KS=True, LogProb=False):
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
		SG: subgroup
		ID: star identifier
		I: inc
		N: Nth planet in system
		dP: period ratio
	'''
	if epos.populationtype == 'parametric':
		
		''' parameters within bounds? '''
		try:
			epos.fitpars.checkbounds(fpara)
		except ValueError as message:
			if Store: raise
			else:
				logging.debug(message)
				return -np.inf
				
		''' Draw (inner) planet from distribution '''
		pps= epos.fitpars.getmc('pps', fpara)
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
		''' Fit parameters (redo) '''
		f_cor= 0.5 # doesn't seem to be well-constrained, light degenracy with f_inc
		
		if len(fpara) is len(epos.groups):
			''' Draw systems from subgroups, array with length survey size '''
			weight= fpara
			f_iso= 0.0 # 0.6 is best fit?
			f_dP= 1.0
			f_inc= 1.0
			if hasattr(epos,'p0'): 
				weight= [epos.p0[0]]
				#f_iso, f_dP, f_inc, f_cor= epos.p0[1:]
				f_iso, f_dP, f_inc= epos.p0[1:]
		else:
			try:
				weight= [fpara[0]]
				f_iso= fpara[1]
				f_dP= fpara[2]
				f_inc= fpara[3]
				#f_cor= fpara[4]
			except:
				raise ValueError('\Wrong number of parameters ({})'.format(len(fpara)))
			
			if not (0<=f_iso<=1) or not (0 <= weight[0] <= 1) or not (0 <= f_cor <= 1) \
				or not (0<=f_dP<=10) or not (0 <= f_inc < 3):
				if Store: raise ValueError('parameters out of bounds')
				return -np.inf		
		
		'''
		Draw from subgroups 
		'''
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
			
			# set allN, Nth planet is system
			
			if 'all_Pratio' in sg: L_dP.append(np.tile(sg['all_Pratio'], ndraw))
			if 'all_R' in sg: L_R.append(np.tile(sg['all_R'], ndraw))

			if epos.Multi: L_I.append(np.tile(sg['all_inc'], ndraw))

		dInc=False # Isotropic inclinations not implemented

		allP=np.concatenate(L_P)
		allM=np.concatenate(L_M)
		allSG=np.concatenate(L_SG)
		allID=np.concatenate(L_ID)
		if epos.Multi: allI=np.concatenate(L_I)
		if 'all_Pratio' in sg: alldP=np.concatenate(L_dP)
		if 'all_R' in sg: allR=np.concatenate(L_R)
		#print len(np.unique(L_allID))/ndraw # == sg['n']
	#if Verbose: print '{:.3f} sec'.format(time.time()-tstart)
			
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
		if epos.populationtype == 'parametric':
			MC_N= allN[itrans] # also for PFM?
		else:
			MC_P*= (1.+0.1*np.random.normal(size=MC_M.size) )
			MC_SG= allSG[itrans]
			if 'all_Pratio' in sg: MC_dP= alldP[itrans]

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
		MC_Y=MC_R
	
	'''
	Store (transiting) planet sample for verification plot
	'''
	if Store:
		tr= epos.transit={}

		tr['P']= MC_P
		tr['Y']= MC_Y
		if epos.populationtype is 'model':
			tr['i sg']= MC_SG
	
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
		if not epos.RV and not epos.populationtype == 'model': det_N= MC_N[idet]
	det_P= MC_P[idet]
	det_Y= MC_Y[idet]

	#if len(alldP)>0: 
	if epos.populationtype is 'model' and 'all_Pratio' in sg: 
		det_dP= MC_dP[idet]

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
	if KS:
		tstart=time.time()
		gof={}
		# make sure that x=P, y=R (where?)
		ix= (epos.xzoom[0]<=det_P) & (det_P<=epos.xzoom[1])
		iy= (epos.yzoom[0]<=det_Y) & (det_Y<=epos.yzoom[1])
		if (ix&iy).sum()<1:
			if Store: raise ValueError('no planets detectable')
			return -np.inf			
	
		_, gof['xvar']=  ks_2samp(epos.obs_zoom['x'], det_P[ix&iy])
		_, gof['yvar']=  ks_2samp(epos.obs_zoom['y'], det_Y[ix&iy])
		# chi^2: (np-nobs)/nobs**0.5 -> p: e^-0.5 x^2
		gof['n']= np.exp(-0.5*(epos.obs_zoom['x'].size-np.sum(ix&iy))**2. 
							/ epos.obs_zoom['x'].size )
		prob= 1.* gof['n']
		if not epos.Multi:
			prob*= gof['xvar']* gof['yvar'] # increase acceptance fraction for models
		elif epos.populationtype == 'parametric' and epos.Multi:
			prob*= gof['xvar']
		
		if epos.Multi:
			# Period ratio
			sim_dP, sim_Pinner= multi.periodratio(det_ID[ix&iy], det_P[ix&iy])
			_, gof['dP']= ks_2samp(epos.obs_zoom['multi']['Pratio'], f_dP*sim_dP)
			prob*= gof['dP']

			# Inner planet
			_, gof['Pinner']= ks_2samp(epos.obs_zoom['multi']['Pinner'], sim_Pinner)
			prob*= gof['Pinner']

			# Multi-planets
			multis= multi.cdf(det_ID[ix&iy], Verbose=False)
			if epos.obs_zoom['multi']['Pratio'].size > 20:
				_, gof['multi']= ks_2samp(epos.obs_zoom['multi']['cdf'], multis)		
			else:
				# not enough statistics
				nobs= epos.obs_zoom['multi']['Pratio'].size
				nsim= sim_dP.size
				gof['multi']= np.exp(-0.5*(nobs-nsim)**2. / nobs )
			
			prob*= gof['multi']
		
		if Verbose:
			print '\nGoodness-of-fit'
			print '  logp= {:.1f}'.format(np.log(prob))
			print '  - p(x)={:.2g}'.format(gof['xvar'])
			print '  - p(y)={:.2g}'.format(gof['yvar'])
			print '  - p(n={})={:.2g}'.format(np.sum(ix&iy), gof['n'])
			if epos.Multi:
				print '  - p(multi)={:.2g}'.format(gof['multi'])
				print '  - p(P ratio)={:.2g}'.format(gof['dP'])
				print '  - p(P inner)={:.2g}'.format(gof['Pinner'])
			

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
		for key, subset in zip(['system', 'single', 'multi'],[isysdet, isingle, imulti]):
			pop[key]={}
			pop[key]['Y']= allY[subset] # R? M?
			pop[key]['P']= allP[subset]
			pop[key]['order']= order[subset]
		
	if Store:
		ss={}
		ss['P']= det_P # MC_P[idet]
		ss['Y']= det_Y
		if epos.populationtype is 'model':
			ss['i sg']= MC_SG[idet]
			if 'all_Pratio' in sg:
				ss['dP']= det_dP
		if epos.MassRadius:
			# or if has mass and radius
			ss['M']= MC_M[idet]
			ss['R']= MC_R[idet]
		
		#if len(alldP)>0
		if epos.Multi:
			ss['multi']={}
			ss['multi']['bin'], ss['multi']['count']= multi.frequency(det_ID[ix&iy])
			ss['multi']['pl cnt']=ss['multi']['bin']* ss['multi']['count']
			ss['multi']['Pratio'], ss['multi']['Pinner']= \
				multi.periodratio(det_ID[ix&iy], det_P[ix&iy]) # *f_dP
			ss['multi']['cdf']= multi.cdf(det_ID[ix&iy])
			if not epos.RV and not epos.populationtype is 'model':
				ss['multi']['PN'],ss['multi']['dPN']= multi.periodratio(
						det_ID[ix&iy], det_P[ix&iy], N=det_N[ix&iy]) 
		
		epos.gof=gof
		ss['P zoom']= det_P[ix&iy]
		ss['Y zoom']= det_Y[ix&iy]
		
		# Store as an extra model 
		if Extra:
			if ~hasattr(epos,'ss_extra'):
				epos.ss_extra=[]
			epos.ss_extra.append(ss)
		else:
			epos.synthetic_survey= ss 
	else:
		# return probability
		with np.errstate(divide='ignore'): 
			lnprob= np.log(prob)
		if np.isnan(lnprob):
			# properly catch errors here?
			#print 'oops III'
			#print gof['p xvar'], gof['p yvar'], gof['p n']
			#print allP.size, MC_P.size, det_P.size
			return -np.inf
		return lnprob if LogProb else prob

# def draw_from_function(f, grid, ndraw, *args):
# 	cdf= np.cumsum(f(grid, *args))
# 	return np.interp(np.random.uniform(0,cdf[-1],ndraw), cdf, grid)

#def make_pdf(epos, norm=None, Init=False):
	
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
		allY[im+(i-1)]= allY[im+(i-2)]* np.random.normal(1, dR, im.size)
		
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

def _pdf(epos, Init=False, fpara=None, xbin=None, ybin=None):
	if fpara is None:
		pps= epos.fitpars.get('pps',Init=Init)
		fpar2d= epos.fitpars.get2d(Init=Init)
	else:
		pps= epos.fitpars.getmc('pps',fpara)
		fpar2d= epos.fitpars.get2d_fromlist(fpara)
		#print fpara
	
	_pdf= epos.func(epos.X_in, epos.Y_in, *fpar2d)
	_pdf_X, _pdf_Y= np.sum(_pdf, axis=1), np.sum(_pdf, axis=0)
	
	# calculate pdf on different grid?
	if (xbin is None) or (ybin is None):
		pdf=_pdf
		pdf_X=_pdf_X
		pdf_Y=_pdf_Y
	else:
		xgrid= np.logspace(np.log10(xbin[0]),np.log10(xbin[-1]),5)
		ygrid= np.logspace(np.log10(ybin[0]),np.log10(ybin[-1]),5)
		X,Y=np.meshgrid(xgrid, ygrid)
		pdf= epos.func(X,Y, *fpar2d)
		pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)

	# normalized to total area
	pdf_X= pps* pdf_X/np.sum(_pdf_X)
	pdf_Y= pps* pdf_Y/np.sum(_pdf_Y)
	pdf= pps* pdf/np.sum(_pdf)
	
	return pps, pdf, pdf_X, pdf_Y
	
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
	
	idx= np.argsort(Pinner)
	rank= np.argsort(idx)

	#Porder= rank[toplanet]
	Porder= 1.0 * rank[toplanet] / IDsys.size
 
	return isysdet, isingle&idetected, imulti&idetected, Porder
