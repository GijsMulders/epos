"""
EPOS classes docstring

This module defines the EPOS class, that contains the observed exoplanets,
the survey detection efficiency, and the synthetic planet population.
see example.py for a simple demonstration of the class
"""

import numpy as np
import cgs
from EPOS import multi
from EPOS.plot.helpers import set_pyplot_defaults

class fitparameters:
	''' Holds the fit parameters. Usually initialized in epos.fitpars '''
	def __init__(self):
		self.fitpars={} # not in order
		self.keysall=[]
		self.keysfit=[]
		self.keys2d=[]
	
	def add(self, key, value, fixed=False, min=-np.inf, max=np.inf, 
				dx=None, text=None, is2D=False):
		'''Add a fit parameter
		
		Args:
			key(str): fit parameter dictionary key
			value(float): starting guess
			fixed(bool): keep this parameter fixed
			min(float): lower bound
			max(float): upper bound
			dx(float): initial dispersion for MCMC
			text(str): plot safe name?
			is2D(bool): use this parameter in the 2D parametric :meth:`EPOS.fitfunctions`
		'''
		fp=self.fitpars[key]= {}
		
		# list of keys
		self.keysall.append(key)
		if is2D: self.keys2d.append(key)
		if not fixed: self.keysfit.append(key)
		
		fp['key']= key
		fp['value_init']= value

		# T/F
		fp['fixed']=fixed
		
		# parameters for fitting
		if not fixed:
			fp['min']=min
			fp['max']=max
			# initial walker positions, can't be zero, default 10% or 0.1 dex
			dx=0.1*value if dx is None else dx
			fp['dx']=abs(dx) if (dx!=0) else 0.1

	def default(self, key, value, Verbose=True):
		if not key in self.keysall: 
			if Verbose: print '  Set {} to default {}'.format(key, value)
			self.add(key, value, fixed=True)
	
	def set(self, key, value):
		self.fitpars[key]['value_fit']=value

	def setfit(self, mclist):
		for i,key in enumerate(self.keysfit):
			#self.fitpars[key]['value_fit']=mclist[self.keysfit.index(key)]
			self.fitpars[key]['value_fit']=mclist[i]
	
	def get(self, key, Init=False, attr=None):
		if attr is None:
			# return initial or fit value
			if Init or not 'value_fit' in self.fitpars[key]:
				return self.fitpars[key]['value_init']
			else:
				return self.fitpars[key]['value_fit']
		else:
			# list of attribute
			return self.fitpars[key][attr]
	
	def get2d(self, Init=False):
		# returns the values for the 2D distribution
		return [self.get(key, Init=Init) for key in self.keys2d]
	
	def getfit(self, Init=True, attr=None): 
		return [self.get(key, Init=Init, attr=attr) for key in self.keysfit]

	def getmc(self, key, parlist):
		# returns value for an mc run
		if self.fitpars[key]['fixed']:
			return self.fitpars[key]['value_init']
		else:
			return parlist[self.keysfit.index(key)]
			#try:				
			#except ValueError:
			#	raise ValueError('Parameter {} not defined'.format(key))

	def get2d_fromlist(self, parlist):
		# returns 2d from fit/fixed, fit supplied in list
		l2d= []
		for key in self.keys2d:
			if self.fitpars[key]['fixed']:
				l2d.append(self.fitpars[key]['value_init'])
				#pps= self.fitpars['pps']['value_init']
			else:
				l2d.append(parlist[self.keysfit.index(key)])
				#pps= parlist[self.keysfit.index('pps')]
		return l2d
	
	def checkbounds(self, parlist):
		for i, key in enumerate(self.keysfit):
			if parlist[i]<self.fitpars[key]['min']:
				raise ValueError('{} out of bounds, {} < {}'.format(
					key,parlist[i],self.fitpars[key]['min']))
			if parlist[i]>self.fitpars[key]['max']:
				raise ValueError('{} out of bounds, {} > {}'.format(
					key,parlist[i],self.fitpars[key]['max']))

class epos:
	"""The epos class
    
    Description:
    	Initialize

    Args:
    	name (str): name to use for directories
        RV(bool): Compare to radial velocity instead of transits
        Debug(bool): Log more output for debugging
        seed(int): Same random number for each simulation? True, None, or int
        Norm(bool): normalize pdf (deprecated?)
    
    Attributes:
		name(str): name
		plotdir(str): plot directory
		RV(bool): Compare to Radial Velocity instead of transit data
		Multi(bool): Do multi-planet statistics
		RandomPairing(bool): multis are randomly paired
		Isotropic(bool): Assume isotropic mutual inclinations
		populationtype(str): type of planet population. 'parametric' or 'model'
		Debug(bool): Verbose logging
		seed(int): Random seed, int or None
    """
	def __init__(self, name, RV=False, Debug=False, seed=True, Norm=False):
		"""
		Initialize the class
		"""
		self.name=name
		self.plotdir='png/{}/'.format(name)
		self.RV= RV

		self.Multi=False
		self.RandomPairing= False
		self.Isotropic= False # phase out?
		
		self.populationtype=None # ['parametric','model']

		# Seed for the random number generator
		if seed is None: self.seed= None
		else:
			if type(seed) is int: self.seed= seed
			else: self.seed= np.random.randint(0, 4294967295)
			print '\nUsing random seed {}'.format(self.seed)
		
		self.Debug= False
		self.Parallel= True # speed up a few calculations 
		set_pyplot_defaults() # nicer plots

		# switches to be set later (undocumented)	
		self.Observation=False
		self.Range=False
		self.DetectionEfficiency=False
		self.Occurrence= False # inverse detection efficiency (?)
		self.Prep= False # ready to run? EPOS.run.once()
		self.MassRadius= False
		self.Radius= False # is this used?
		
		self.plotpars={} # dictionary to hold some customization keywords
		
	def set_observation(self, xvar, yvar, starID, nstars=1.6862e5):
		''' Observed planet population
		
		Args:
			xvar: planet orbital period [list]
			yvar: planet radius or M sin i [list]
			ID: planet ID [list]
			nstars: number of stars surveyed
		
		Note:
			Some pre-defined planet populations from Kepler can be generated from 
			:mod:`EPOS.kepler`
		'''
		order= np.lexsort((xvar,starID)) # sort by ID, then P
		self.obs_xvar=np.asarray(xvar)[order]
		self.obs_yvar=np.asarray(yvar)[order]
		self.obs_starID=np.asarray(starID)[order]
		self.nstars=nstars
		
		assert self.obs_xvar.ndim == self.obs_yvar.ndim == self.obs_starID.ndim == 1, 'only 1D arrays'
		assert self.obs_xvar.size == self.obs_yvar.size == self.obs_starID.size, 'arrays not same length'
	
		# set plot limits in observation 5% wider than data
		xmin,xmax= min(self.obs_xvar), max(self.obs_xvar)
		ymin,ymax= min(self.obs_yvar), max(self.obs_yvar)
		dx= (xmax/xmin)**0.05
		dy= (ymax/ymin)**0.05
		self.obs_xlim=[xmin/dx,xmax*dx]
		self.obs_ylim=[ymin/dy,ymax*dy]
		
		self.Observation=True
		
		# print some stuff
		print '\nObservations:\n  {} stars'.format(int(nstars))
		print '  {} planets'.format(self.obs_starID.size)
		multi.indices(self.obs_starID, Verbose=True)
		epos.multi={}
		epos.multi['bin'], epos.multi['count']= \
			multi.frequency(self.obs_starID, Verbose=True)
		epos.multi['pl cnt']= epos.multi['bin']* epos.multi['count']
		epos.multi['Pratio'], epos.multi['Pinner']= \
			multi.periodratio(self.obs_starID, self.obs_xvar, Verbose=True)
		epos.multi['cdf']= multi.cdf(self.obs_starID, Verbose=True)	
		
	def set_survey(self, xvar, yvar, eff_2D, Rstar=1.0, Mstar=1.0):
		'''Survey detection efficiency (completeness)
		Args:
			xvar: planet orbital period grid [list]'
			yvar: planet radius or M sin i grid [list]
			eff_2D: 2D matrix of detection efficiency
			Rstar: stellar radius, for calculating transit probability
			Mstar: stellar mass, for period-semimajor axis conversion
		
		Note:
			Some pre-defined detection efficiencies from Kepler can be generated from 
			:mod:`EPOS.kepler`
		'''
		self.eff_xvar=np.asarray(xvar)
		self.eff_yvar=np.asarray(yvar)
		self.eff_2D=np.asarray(eff_2D)
		
		assert self.eff_xvar.ndim == self.eff_yvar.ndim == 1, 'only 1D arrays'
		assert self.eff_2D.ndim == 2, 'Detection efficiency must by a 2dim array'
		if self.eff_2D.shape != (self.eff_xvar.size, self.eff_yvar.size):
			raise ValueError('Mismatching detection efficiency'
			'\n: nx={}, ny={}, (nx,ny)=({},{})'.format(self.eff_xvar.size, self.eff_yvar.size, *self.eff_2D.shape))
	
		self.eff_xlim= [min(self.eff_xvar),max(self.eff_xvar)]
		self.eff_ylim= [min(self.eff_yvar),max(self.eff_yvar)]
		
		self.Mstar= Mstar
		if self.RV:
			self.completeness= self.eff_2D
		else:
			self.Rstar=Rstar # Solar radii
			self.Pindex= -2./3.
			fourpi2_GM= 4.*np.pi**2. / (cgs.G*self.Mstar*cgs.Msun)
			self.fgeo_prefac= self.Rstar*cgs.Rsun * fourpi2_GM**(1./3.) / cgs.day**(2./3.)
			#print self.fgeo_prefac
			P, R= np.meshgrid(self.eff_xvar, self.eff_yvar, indexing='ij')
			self.completeness= self.eff_2D * self.fgeo_prefac*P**self.Pindex

		self.DetectionEfficiency=True
	
	def set_ranges(self, xtrim=None, ytrim=None, xzoom=None, yzoom=None):
		
		if self.Range: raise ValueError('Range already defined')
		if not self.Observation: raise ValueError('No observation defined')
		if not self.DetectionEfficiency: raise ValueError('No detection effifiency defined')
		
		''' Define the region where completeness is calculated'''
		if xtrim is None:
			print 'Trimming x-axis from detection efficiency'
			self.xtrim= self.eff_xlim
		else:
			self.xtrim= [max(xtrim[0], self.eff_xlim[0]), min(xtrim[1], self.eff_xlim[1])]

		if ytrim is None:
			print 'Trimming y-axis from detection efficiency'
			self.ytrim= self.eff_ylim
		else:
			self.ytrim= [max(ytrim[0], self.eff_ylim[0]), min(ytrim[1], self.eff_ylim[1])]
		
		''' Define a smaller region where observational comparison is performed'''	
		if xzoom is None:
			print 'Not zooming in on x-axis for model comparison'
			self.xzoom= self.xtrim
		else:
			self.xzoom= [max(xzoom[0], self.xtrim[0]), min(xzoom[1], self.xtrim[1])]

		if yzoom is None:
			print 'Not zooming in on y-axis for model comparison'
			self.yzoom= self.ytrim
		else:
			self.yzoom= [max(yzoom[0], self.ytrim[0]), min(yzoom[1], self.ytrim[1])]
		
		if (xzoom is None) and (yzoom is None):
			self.Zoom=False
		elif (self.xzoom==self.xtrim) and (self.yzoom==self.ytrim):
			self.Zoom=False
			print 'Not a zoom'
		else:
			self.Zoom=True
		
		''' Prep the grid for the observable'''
		# make sure range _encompasses_ trim
		ixmin,ixmax= _trimarray(self.eff_xvar, self.xtrim)
		iymin,iymax= _trimarray(self.eff_yvar, self.ytrim)
	
		self.MC_xvar= self.eff_xvar[ixmin:ixmax]
		self.MC_yvar= self.eff_yvar[iymin:iymax]
		self.MC_eff= self.eff_2D[ixmin:ixmax,iymin:iymax]
		
		# scale factor to multiply pdf such that occurrence in units of dlnR dlnP
		self.scale_x= self.MC_xvar.size/np.log(self.MC_xvar[-1]/self.MC_xvar[0])
		self.scale_y= self.MC_yvar.size/np.log(self.MC_yvar[-1]/self.MC_yvar[0])	
		self.scale= self.scale_x * self.scale_y
		
		self.X, self.Y= np.meshgrid(self.MC_xvar, self.MC_yvar, indexing='ij')
		
		''' Prep the grid for the PDF, if using a mass-radius conversion'''
		if self.populationtype is 'parametric':
			if self.MassRadius:
				self.in_ytrim= self.masslimits
				self.in_yvar= np.logspace(*np.log10(self.in_ytrim))
				self.scale_in_y= \
					self.in_yvar.size/np.log(self.in_yvar[-1]/self.in_yvar[0])	
				self.scale_in= self.scale_x * self.scale_in_y
				self.X_in,self.Y_in= np.meshgrid(self.MC_xvar,self.in_yvar,indexing='ij')
			else:
				self.in_ytrim= self.ytrim
				self.in_yvar= self.MC_yvar
				self.scale_in_y= self.scale_y 
				self.scale_in= self.scale
				self.X_in, self.Y_in= self.X, self.Y
			
		''' plot ticks '''
		if self.RV:
			self.xticks= [1,10,100]
			self.yticks= [1,10,100,1000]
		else:
			self.xticks= [1,10,100,1000]
			self.yticks= [0.5,1,4,10]
			
		self.Range=True
	
	def set_bins(self, xbins=[[1,10]], ybins=[[1,10]], Sparse=False):
		focc= self.occurrence={}
		assert Sparse==False
		if not (np.ndim(xbins) == np.ndim(ybins) == 2):
			raise ValueError('wrong bin dimensions')
		assert len(xbins) == len(ybins)
		assert np.shape(xbins) == np.shape(ybins) 
		
		focc['bin']={}
		focc['bin']['x']= np.array(xbins)
		focc['bin']['y in']= np.array(ybins)

	def set_parametric(self, func):
		'''Define a parametric function to generate the planet size-period distribution
		
		Description:
			Function should be callable as func(X, Y, \*fitpars2d) with
			X(np.array): period
			Y(np.array): size (radius or mass)
			The list of fit parameters fitpars2d will be constructed from parameters 
			added using :func:`EPOS.fitparameters.add` with is2D=True
			
		Note:
			Some pre-defined functions can be found in :mod:`EPOS.fitfunctions`
		
		Args:
			func (function): callable function
		
		'''
		if self.populationtype is None:
			self.populationtype='parametric'
		elif self.populationtype is not 'parametric':
			raise ValueError('You have already defined a planet population ({})'.format(self.populationtype))
		
		if not callable(func): raise ValueError('func is not a callable function')
		
		self.func=func
		self.fitpars= fitparameters()
		
	def set_multi(self, spacing=None):
		if self.populationtype is not 'parametric':
			raise ValueError('Define a parametric planet population first')
		self.Multi=True
		
		self.RandomPairing= (spacing==None)
		self.spacing= spacing # None, brokenpowerlaw, dimensionless
	
	def add_population(self, name, sma, mass, 
					inc=None, tag1=None, Verbose=False, weight=1.):
		# tag is fit parameter, i.e. metallicity or surface density
		if self.populationtype is None:
			self.populationtype='model'
			self.groups=[]
			self.mod_xlim=[np.inf,0]
			self.mod_ylim=[np.inf,0]
		elif self.populationtype is not 'model':
			raise ValueError('You have already defined a planet population ({})'.format(self.populationtype))
		
		print '\nLoading subgroup {} '.format(name)
		
		try:
			if len(sma) is not len(mass): raise ValueError('sma ({}) and mass ({}) not same length'.format(sma,mass))
		except: raise ValueError('sma ({}) and mass ({}) have to be iterable'.format(type(sma), type(mass)))
		
		# model has mutual inclinations?
		Inc= inc is not None
		if not Inc:
			self.Multi=False
		else:
			self.Multi=True # no mix of T/F
		
		sg={}
		sg['name']= name
		sg['n']= len(sma)
		sg['system']=[] # remains empty if list
		sg['weight']= weight

		# list of lists or 1-dim list?
		#print type(sma[0]) 
		if type(sma[0]) is list or type(sma[0]) is np.ndarray:
			# loop over systems
			for k, (L_sma, L_mass) in enumerate(zip(sma, mass)):
				#assert type(L_sma) is list
				sg['system'].append({})
				order= np.argsort(L_sma)
				_sma=	sg['system'][-1]['sma']= np.array(L_sma)[order]
				_P=		sg['system'][-1]['P']= _sma**1.5 * 365.25 # 
				_mass=	sg['system'][-1]['mass']= np.array(L_mass)[order]
				sg['system'][-1]['np']= len(L_sma)
				sg['system'][-1]['ID']= np.array([k+1]*len(L_sma))
				
				# quick ratio, nan for first planet
				dP= sg['system'][-1]['P ratio']= np.full_like(_P, np.nan)
				if _P.size>1:
					dP[1:]= _P[1:]/_P[:-1]
				
				# inc (can't zip if variable is None)
				if Inc: sg['system'][-1]['inc']= np.array(inc[k])[order]
				
			sg['all_sma']= 	np.concatenate([plsys['sma'] for plsys in sg['system']])
			sg['all_mass']=	np.concatenate([plsys['mass'] for plsys in sg['system']])
			sg['all_P']= 	np.concatenate([plsys['P'] for plsys in sg['system']])
			sg['all_ID']= 	np.concatenate([plsys['ID'] for plsys in sg['system']])
			sg['all_Pratio']= np.concatenate([plsys['P ratio'] for plsys in sg['system']])
			if Inc: sg['all_inc']= np.concatenate([plsys['inc'] for plsys in sg['system']])
				
		else:
			sg['all_sma']= 		np.asarray(sma)
			sg['all_mass']=		np.asarray(mass)
			#sg['all_ID']= 	
			if Inc: sg['inc']=	np.asarray(inc)
			if tag1 is not None: sg['all_tag1']= np.asarray(tag1)
	
		#sg['all_P']= sg['all_sma']**1.5 * 365.25
			
		print '  {} planetary systems'.format(sg['n'])
		if Verbose:
			# print '\nLoaded subgroup {} with {} planetary systems'.format(sg['name'], sg['n'])
			for system in sg['system']:
				print 'system has {} planets, {:.1f} Mearth:'.format(
					system['np'],np.sum(system['mass']))
				for a,m in zip(system['sma'], system['mass']):
					print '  {:.2f} au, {:.1f} Mearth'.format(a,m)
				#print sg['all_Pratio']
				#print
		
		# set plot limits in model 5% wider than data
		xmin,xmax= min(sg['all_sma']), max(sg['all_sma'])
		ymin,ymax= min(sg['all_mass']), max(sg['all_mass'])
		dx= (xmax/xmin)**0.05
		dy= (ymax/ymin)**0.05
		self.mod_xlim=[min(xmin/dx, self.mod_xlim[0]),max(xmax*dx, self.mod_xlim[1])]
		self.mod_ylim=[min(ymin/dy, self.mod_ylim[0]),max(ymax*dy, self.mod_ylim[1])]
		
		self.groups.append(sg)
	
	def set_massradius(self, MR, name, masslimits= [0.01,1e3]):
		print '\nMass-Radius relation from {}'.format(name)
		if self.MassRadius:
			raise ValueError('Already defined a Mass-Radius conversion function ')
		# actually radius as funtion of mass (mass-to-radius)
		if not callable(MR): raise ValueError('Mass-Radius function is not callable')

		self.MassRadius=True
		self.MR=MR
		self.MR_label=name
		
		self.masslimits=masslimits 
		meanradius= MR(masslimits)[0]
		print 'Mass and Radius limits:'
		print '  min M = {:.3f}-> <R> ={:.2f}'.format(masslimits[0], meanradius[0] )
		print '  max M = {:.0f}-> <R> ={:.1f}'.format(masslimits[-1], meanradius[-1] )
		

def _trimarray(array,trim):
	# trims array of points not needed for interpolation
	if trim[0] < array[1]:
		imin=0
	else:
		imin = np.searchsorted(array, trim[0], side='left')-1
	
	if trim[1] > array[-2]:
		imax=len(array)
	else:
		imax = np.searchsorted(array, trim[1], side='left')+1
	
	return imin, imax
