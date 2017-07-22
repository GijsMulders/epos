import numpy as np
import time
from scipy import interpolate
from scipy.stats import ks_2samp
import os, logging
import cgs
import multi
from functools import partial
from EPOS.fitfunctions import brokenpowerlaw1D

#reload(logging)

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
		fpara= epos.p0[epos.ip_fit]
	elif epos.populationtype is 'model':
		fpara= [sg['weight'] for sg in epos.groups]
		if hasattr(epos,'p0'): fpara=epos.p0
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
		fpara= epos.p0[epos.ip_fit]
	elif epos.populationtype is 'model':
		if len(epos.groups) == 1:
			fpara=epos.p0 # [epos.ipfit] ??
		else:
			fpara= [sg['weight'] for sg in epos.groups]
	else: assert False
	
	''' Load previous chain?'''
	ndim= len(fpara)
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
		p0 = [np.array(fpara)*np.random.uniform(1.-dx,1+dx,len(fpara)) 
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
		np.save(fname, epos.chain) 
		
	''' the posterior samples after burn-in '''
	epos.samples= epos.chain[:, nburn:, :].reshape((-1, ndim))
	epos.burnin= nburn
	fitpars = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(epos.samples, [16, 50, 84],
                                                axis=0)))
	if epos.populationtype is 'parametric':
		epos.pfit= np.copy(epos.p0)
		#print 'before: {}'.format(epos.pfit)
		epos.pfit[epos.ip_fit]=[p[0] for p in fitpars]
		#print 'after: {}'.format(epos.pfit)
	else:
		epos.pfit=[p[0] for p in fitpars]
	
	''' estimate #planets/star '''
	#print epos.samples.shape # 30000, 7
	if epos.populationtype is 'parametric' and epos.Isotropic:
		# reconstruct posterior (from p0 and samples?)
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
	
	z['multi']['Pratio'], z['multi']['Pinner']= \
		multi.periodratio(epos.obs_starID[ix&iy], epos.obs_xvar[ix&iy])
	z['multi']['cdf']= multi.cdf(epos.obs_starID[ix&iy])

def MC(epos, fpara, Store=False, Verbose=True, KS=True, LogProb=False):
	''' 
	Do the Monte Carlo Simulations
	Note:
	variable x/X is P
	variable y/Y is R/M
	'''	
	
	if Verbose: tstart=time.time()
	
	#if not Store: logging.debug(' '.join(['{:.3g}'.format(fpar) for fpar in fpara]))
	
	''' Seed the random number generator '''
	if epos.seed is not None: np.random.seed(epos.seed)
	
	''' construct 1D arrays for P, R
	dimension equal to sample size * planets_per_star
	also keeping track off subgroup (SG), inc (I), and period ratio (dP)
	'''
	if epos.populationtype == 'parametric':
		
		''' Draw (inner) planet from distribution '''
		fpar2d= np.append(fpara[epos.ip2d_fit], epos.p0[epos.ip2d_fixed])
		
		try: sysX, sysY= draw_from_2D_distribution(epos, fpar2d)
		except ValueError:
			if Store: raise
			else: return -np.inf
		
		''' Add multi-planet systems? '''
		if epos.Isotropic:
			allX=sysX
			allY=sysY
		else:
			''' retrieve fit parameters for multi-planets (no dR)'''
			try: npl, dPbreak, dP1, dP2, dInc, f_iso= fpara[-6:]
			except: raise ValueError('Incorrect number of parameters?')
			f_dP, f_inc= 1.0, 1.0 # no need to fudge these
			f_SNR=0.5
			dR=0.01 # epos.p0[epos.ip_fixed[-1:]]
			
			''' Parameter bounds '''
			if npl < 1:
				logging.debug('npl = {:.3g} < 1'.format(npl))
				if Store: raise ValueError('at least one planet per system')
				return -np.inf
			if npl > 10:
				logging.debug('npl = {:.1f} > 10'.format(npl))
				if not Store: return -np.inf # out of memory
			if (dPbreak<1) or (dPbreak > 10):
				logging.debug('dP break = {:.3g} not between 1 and 10'.format(dPbreak))
				if Store: print ' 1 < Pbreak < 10'
				else: return -np.inf
			if (dInc <=0) or (dR <=0) or not (0 <= f_iso <= 1):
				if Store: raise ValueError('parameters out of bounds')
				return -np.inf
			
			# 
			if len(fpara[:-6]) is 7:
				# P2 > 1
				if (fpara[3]> 1):
					logging.debug('{} = {} > 1'.format(epos.pname[3], fpara[3]))
					if Store: raise ValueError('P2 < 1 for multis')
					return -np.inf
					
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
			
			''' Adjust parameter of 2nd, 3rd planet etc.'''
			Pgrid= np.logspace(0,1)
			cdf= np.cumsum(brokenpowerlaw1D(Pgrid, dPbreak, dP1, dP2))
			# loop over planet 2,3... n
			for i in range(2,len(np.bincount(sysnpl)) ):
				im= i1[sysnpl>=i] # index to ith planet in each system
				#print '  planet {} in system, {} systems'.format(i,im.size)
				#dP= draw_from_function(brokenpowerlaw1D,xgrid,im.size, dPbreak, dP1, dP2)
				dP= np.interp(np.random.uniform(0,cdf[-1],im.size), cdf, Pgrid)
				allX[im+(i-1)]= allX[im+(i-2)]* dP
				#np.random.uniform(dP-0.7, dP+0.7, im.size)
				#print allX[:3]
				#np.random.norm()
				allY[im+(i-1)]= allY[im+(i-2)]* np.random.normal(1, dR, im.size)
			#print '+{}={} story checks out?'.format(allID[i1].size, allID.size)
			#print allX[:3]
			
			''' Toss out planets (reduces memory footprint) '''
			Xinrange= allX<=epos.MC_xvar[-1]
			#if Xinrange.sum() < 0.5* allX.size: 
			#	print 'yeah, {}/{}'.format(Xinrange.sum(),allX.size)
			allX= allX[Xinrange]
			allY= allY[Xinrange]
			allI= allI[Xinrange]
			allID= allID[Xinrange]
		
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
			
			if not (0<=f_iso<=1) or not (0 <= weight[0] <= 1) or not (0 <= f_SNR <= 1) \
				or not (0<=f_dP<=10) or not (0 <= f_inc < 3):
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
	#if Verbose: print '{:.3f} sec'.format(time.time()-tstart)
		
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

		if epos.RadiusMassConversion:
			MC_M= allM[itrans]
		else:
			# parametric + not isotropic
			MC_R= allR[itrans]
		
		MC_ID= allID[itrans]
		if epos.populationtype == 'parametric':
			MC_P= allP[itrans]
		else:
			MC_P= allP[itrans] * (1.+0.1*np.random.normal(size=MC_M.size) )
			MC_SG= allSG[itrans]
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
		if (ix&iy).sum()<1:
			if Store: raise ValueError('no planets detectable')
			return -np.inf			
	
		_, gof['xvar']=  ks_2samp(epos.obs_zoom['x'], det_P[ix&iy])
		_, gof['yvar']=  ks_2samp(epos.obs_zoom['y'], det_Y[ix&iy])
		# chi^2: (np-nobs)/nobs**0.5 -> p: e^-0.5 x^2
		gof['n']= np.exp(-0.5*(epos.obs_zoom['x'].size-np.sum(ix&iy))**2. 
							/ epos.obs_zoom['x'].size )
		prob= 1.* gof['n']
		if epos.Isotropic:
			prob*= gof['xvar']* gof['yvar'] # increase acceptance fraction for models
		elif epos.populationtype == 'parametric' and not epos.Isotropic:
			prob*= gof['xvar']
		
		if not epos.Isotropic:
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
			if not epos.Isotropic:
				print '  - p(multi)={:.2g}'.format(gof['multi'])
				print '  - p(P ratio)={:.2g}'.format(gof['dP'])
				print '  - p(P inner)={:.2g}'.format(gof['Pinner'])
			

		tgof= time.time()
		if Verbose: print '  observation comparison in {:.3f} sec'.format(tgof-tstart)


	if Store:	
		ss= epos.synthetic_survey= {}
		ss['P']= det_P # MC_P[idet]
		ss['Y']= det_Y
		if epos.populationtype is 'model':
			ss['i sg']= MC_SG[idet]
			if 'all_Pratio' in sg:
				ss['dP']= det_dP
		if epos.RadiusMassConversion:
			# or if has mass and radius
			ss['M']= MC_M[idet]
			ss['R']= MC_R[idet]
		
		#if len(alldP)>0
		if not epos.Isotropic:
			ss['multi']={}
			ss['multi']['bin'], ss['multi']['count']= multi.frequency(det_ID[ix&iy])
			ss['multi']['Pratio'], ss['multi']['Pinner']= \
				multi.periodratio(det_ID[ix&iy], det_P[ix&iy]) # *f_dP
			ss['multi']['cdf']= multi.cdf(det_ID[ix&iy])
		
		epos.gof=gof
		ss['P zoom']= det_P[ix&iy]
		ss['Y zoom']= det_Y[ix&iy]
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

def draw_from_2D_distribution(epos, fpara):

		''' Check if parameters make sense'''
		# normalization can't be zero
		if fpara[0] <= 0: 
			raise ValueError('normalization needs to be larger than 0')

		# hack for double power law
		if len(fpara) is 7:
			if (fpara[1]<= 0) or (fpara[4]<= 0):
				raise ValueError('breaks needs to be larger than 0')
		
		''' create PDF, CDF'''
		# assumes a seperable function of mass and radius
		pdf= epos.func(epos.X, epos.Y, *fpara)
		pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
		cum_X, cum_Y= np.cumsum(pdf_X), np.cumsum(pdf_Y)
		pps_x, pps_y=  cum_X[-1], cum_Y[-1]
		planets_per_star= 0.5*(pps_x+pps_y) # should be equal
		
		try:
			ndraw= int(round(planets_per_star*epos.nstars))
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
		elif planets_per_star > 100:
			logging.debug('>100 planets per star ({})'.format(planets_per_star))
			raise ValueError('too many planets per star')
		try:
			allX= np.interp(np.random.uniform(0,pps_x,ndraw), cum_X, epos.MC_xvar)
			allY= np.interp(np.random.uniform(0,pps_y,ndraw), cum_Y, epos.MC_yvar)
		except MemoryError:
			raise ValueError('Memory error for n={}'.format(ndraw))
		
		return allX, allY
