"""
EPOS classes docstring

This module defines the EPOS class, that contains the observed exoplanets,
the survey detection efficiency, and the synthetic planet population.
see example.py for a simple demonstration of the class
"""

import numpy as np

import cgs
import EPOS.multi
from EPOS.plot.helpers import set_pyplot_defaults
from EPOS import __version__

class fitparameters:
	''' Holds the fit parameters. Usually initialized in epos.fitpars '''
	def __init__(self):
		self.fitpars={} # not in order
		self.keysall=[]
		self.keysfit=[]
		self.keys2d=[]
		self.keypps='pps'
	
	def add(self, key, value, fixed=False, min=-np.inf, max=np.inf, 
				dx=None, text=None, is2D=False, isnorm=False):
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
			isnorm(bool): this parameter is the normalization factor for the number of planet per star :meth:`EPOS.fitfunctions`
		'''
		fp=self.fitpars[key]= {}
		
		# list of keys
		self.keysall.append(key)
		if is2D: self.keys2d.append(key)
		if isnorm: self.keypps=key
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
		self.fitpars[key]['value_init']=value

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

	def getpps(self, Init=False):
		# returns the normalization for the 2D distribution
		return self.get(self.keypps, Init=Init)

	def getpps_fromlist(self, parlist):
		return self.getmc(self.keypps, parlist)
	
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

	def get_kfree(self):
		return len(self.keysfit)

class epos:
	"""The epos class
	
	Description:
		Initialize

	Args:
		name (str): name to use for directories
		RV(bool): Compare to radial velocity instead of transits
		MC(bool): Generate planet population by random draws 
		Msini(bool): Convert the planet mass distribution into an Msini distribution
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
		Parametric(bool): parametric planet population?
		Debug(bool): Verbose logging
		seed(): Random seed, can be any of int, True, or None
	"""
	def __init__(self, name, Debug=False, seed=True, title=None, 
		RV=False, Norm=False, MC=True, Msini=False):
		"""
		Initialize the class
		"""
		self.name=name
		self.title=name if title is None else title

		print '\n\n |~| epos {} |~|\n'.format(__version__)

		''' Directories '''
		self.plotdir='png/{}/'.format(name)
		self.jsondir='json/{}/'.format(name)
		#self.path= os.path.dirname(EPOS.__file__)

		''' EPOS mode'''
		self.RV= RV
		self.Msini= Msini # do an M -> Msini conversion
		self.Multi=False
		self.RandomPairing= False
		self.Isotropic= False # phase out?
		self.MonteCarlo= MC
		
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
		self.PDF=False
		
		self.plotpars={} # dictionary to hold some customization keywords

	def set_observation(self, xvar, yvar, starID, nstars=1.6862e5, 
		radiusError=0.1, score=None):
		''' Observed planet population
		
		Args:
			xvar: planet orbital period [list]
			yvar: planet radius or M sin i [list]
			ID: planet ID [list]
			nstars: number of stars surveyed
		
		Note:
			Some pre-defined planet populations from Kepler can be generated from 
			:mod:`EPOS.kepler.dr25()`
		'''
		order= np.lexsort((xvar,starID)) # sort by ID, then P
		self.obs_xvar=np.asarray(xvar)[order]
		self.obs_yvar=np.asarray(yvar)[order]
		self.obs_starID=np.asarray(starID)[order]
		self.nstars=nstars

		if score is not None:
			self.obs_score= np.asarray(score)[order]
		
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
		self.radiusError= radiusError
		
		# print some stuff
		print '\nObservations:\n  {} stars'.format(int(nstars))
		print '  {} planets'.format(self.obs_starID.size)
		EPOS.multi.indices(self.obs_starID, Verbose=True)
		epos.multi={}
		epos.multi['bin'], epos.multi['count']= \
			EPOS.multi.frequency(self.obs_starID, Verbose=True)
		epos.multi['pl cnt']= epos.multi['bin']* epos.multi['count']
		epos.multi['Pratio'], epos.multi['Pinner'], epos.multi['Rpair']= \
			EPOS.multi.periodratio(self.obs_starID, self.obs_xvar, R=self.obs_yvar, Verbose=True)
		epos.multi['cdf']= EPOS.multi.cdf(self.obs_starID, Verbose=True)	
	
	def set_survey(self, xvar, yvar, eff_2D, Rstar=1.0, Mstar=1.0, vet_2D=None):
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
			raise ValueError('Mismatching detection efficiency'+
			'\n: nx={}, ny={}, (nx,ny)=({},{})'.format(self.eff_xvar.size, 
				self.eff_yvar.size, *self.eff_2D.shape))
	
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
		
		if vet_2D is not None:
			self.vetting= np.asarray(vet_2D)
			
			if self.eff_2D.shape != (self.eff_xvar.size, self.eff_yvar.size):
				raise ValueError('Mismatching vetting efficiency'+
					'\n: (nx,ny)={}, (nx,ny={})'.format(self.eff_2D.shape,
								self.vetting.shape))

			self.completeness_novet= self.completeness
			self.completeness*= self.vetting

		self.DetectionEfficiency=True
	
	def set_ranges(self, xtrim=None, ytrim=None, xzoom=None, yzoom=None, 
		LogArea=False, Occ=False, UnitTicks=True, plotxgrid= None, plotygrid=None):
		
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
			print 'Trim equal to zoom' # so??
		else:
			self.Zoom=True
		
		''' 
		Prep the Monte Carlo grid for the observable
		'''
		# make sure range _encompasses_ trim
		ixmin,ixmax= _trimarray(self.eff_xvar, self.xtrim)
		iymin,iymax= _trimarray(self.eff_yvar, self.ytrim)
		self.trim_to_zoom= (slice(ixmin,ixmax),slice(iymin,iymax))
	
		self.MC_xvar= self.eff_xvar[ixmin:ixmax]
		self.MC_yvar= self.eff_yvar[iymin:iymax]
		self.MC_eff= self.eff_2D[ixmin:ixmax,iymin:iymax]
		if hasattr(self,'vetting'):
			self.MC_eff*= self.vetting[ixmin:ixmax,iymin:iymax]
		
		# scale factor to multiply pdf such that occurrence in units of dlnR dlnP
		if LogArea:
			area= np.log10
			self.plotpars['area']= 'dlog'
		else:
			area= np.log
			self.plotpars['area']= 'dln'

		self.scale_x= self.MC_xvar.size/area(self.MC_xvar[-1]/self.MC_xvar[0])
		self.scale_y= self.MC_yvar.size/area(self.MC_yvar[-1]/self.MC_yvar[0])
		self.scale= self.scale_x * self.scale_y
		
		self.X, self.Y= np.meshgrid(self.MC_xvar, self.MC_yvar, indexing='ij')
		
		''' Transit probability '''
		if not self.RV:
			self.f_geo= self.fgeo_prefac *self.X**self.Pindex
			self.f_det= self.f_geo * self.MC_eff
		
		'''
		Prep the grid for the non-MC simulation
		'''
		if not self.MonteCarlo:
			self.noMC_zoom_x=np.geomspace(self.xzoom[0],self.xzoom[-1])
			self.noMC_zoom_y=np.geomspace(self.yzoom[0],self.yzoom[-1])
			self.noMC_scale_x= self.noMC_zoom_x.size/ np.log(self.xzoom[-1]/self.xzoom[0])
			self.noMC_scale_y= self.noMC_zoom_y.size/ np.log(self.yzoom[-1]/self.yzoom[0])
			self.noMC_scale= self.noMC_scale_x * self.noMC_scale_y
		
		''' Prep the grid for the PDF, if using a mass-radius conversion'''
		if self.PDF:
			if self.MassRadius:
				self.in_ytrim= self.masslimits
				self.in_yvar= np.logspace(*np.log10(self.in_ytrim))
				self.scale_in_y= \
					self.in_yvar.size/area(self.in_yvar[-1]/self.in_yvar[0])	
				self.scale_in= self.scale_x * self.scale_in_y
				self.X_in,self.Y_in= np.meshgrid(self.MC_xvar,self.in_yvar,indexing='ij')
			elif self.RV and self.Msini:
				# extend grid to higher masses
				nex=13
				ext= self.MC_yvar[-nex:]*self.MC_yvar[-1]/self.MC_yvar[-(nex+1)]
				self.in_yvar= np.concatenate([self.MC_yvar, ext]) 
				self.in_ytrim= self.in_yvar[([0,-1])]
				self.scale_in_y= self.scale_y 
				self.scale_in= self.scale
				self.X_in, self.Y_in= np.meshgrid(self.MC_xvar,self.in_yvar,indexing='ij')	
			else:
				self.in_ytrim= self.ytrim
				self.in_yvar= self.MC_yvar
				self.scale_in_y= self.scale_y 
				self.scale_in= self.scale
				self.X_in, self.Y_in= self.X, self.Y
			
		''' plot ticks '''
		yr= 365.24
		if UnitTicks:
			Pticks=np.array([1,10,100, yr, 10.*yr])
			Pticklabels= np.array(['1','10','100','1yr','10yr'])
		else:
			Pticks= np.logspace(-1,6,8)
			Pticklabels= Pticks

		xrange= (self.xtrim[0]<=Pticks) & (Pticks<=self.xtrim[1])
		self.xticks= Pticks[xrange]
		self.xticklabels= Pticklabels[xrange]
		
		if UnitTicks:
			Mjup= cgs.Mjup/cgs.Mearth
			Mticks= np.array([0.1,1,10,100, Mjup, 10.*Mjup, 100.*Mjup])
			Mticklabels= np.array(['0.1','1','10','100', '$M_J$', '$10 M_J$', '$100 M_J$'])
			Rticks= np.array([0.25,0.5,1,2, 4, 10])
			Rticklabels= np.array(['0.25','0.5','1','2', '4','10'])
		else:
			Mticks= np.logspace(-1,6,8)
			Mticklabels= Mticks
			Rticks= np.logspace(-2,4,7)
			Rticklabels= Rticks

		if self.MassRadius:
			yinrange= (self.in_ytrim[0]<=Mticks) & (Mticks<=self.in_ytrim[1])
			self.y_inticks= Mticks[yinrange]
			self.y_inticklabels= Mticklabels[yinrange]
			yrange= (self.ytrim[0]<=Rticks) & (Rticks<=self.ytrim[1])
			self.yticks= Rticks[yrange]
			self.yticklabels= Rticklabels[yrange]			
		elif self.RV:
			yrange= (self.ytrim[0]<=Mticks) & (Mticks<=self.ytrim[1])
			self.yticks= Mticks[yrange]
			self.yticklabels= Mticklabels[yrange]
			yrange= (self.in_ytrim[0]<=Mticks) & (Mticks<=self.in_ytrim[1])
			self.y_inticks= Mticks[yrange]
			self.y_inticklabels= Mticklabels[yrange]
		else:			
			yrange= (self.ytrim[0]<=Rticks) & (Rticks<=self.ytrim[1])
			self.yticks= Rticks[yrange]
			self.yticklabels= Rticklabels[yrange]
			
		self.Range=True
		
		if Occ:
			if self.MassRadius:
				raise ValueError('Plotting occurrence with mass-radius not yet supported')
				#pass
				
			if not hasattr(self,'occurrence'):
				self.occurrence={}
			focc= self.occurrence
			
			if self.RV:
				xgrid= np.geomspace(self.MC_xvar[0], self.MC_xvar[-1], num=20)
				ygrid= np.geomspace(self.MC_yvar[0], self.MC_yvar[-1], num=20)
			else:
				xgrid= self.MC_xvar
				ygrid= self.MC_yvar	
			
			# custom grid for *plotting* occurrence rates
			if plotxgrid is not None: xgrid= plotxgrid
			if plotygrid is not None: ygrid= plotygrid

			focc['xzoom']={}
			#focc['xzoom']['grid']= xgrid
			#ygrid= np.exp(np.arange(np.log(self.MC_yvar[0]),np.log(self.MC_yvar[-1])+0))
			focc['xzoom']['x']= [self.xzoom]* (ygrid.size-1)
			focc['xzoom']['y']= [[i,j] for i,j in zip(ygrid[:-1],ygrid[1:])]

			focc['yzoom']={}
			#focc['yzoom']['grid']= ygrid
			#xgrid= np.exp(np.arange(np.log(self.MC_xvar[0]),np.log(self.MC_xvar[-1])+0))
			focc['yzoom']['x']= [[i,j] for i,j in zip(xgrid[:-1],xgrid[1:])]
			focc['yzoom']['y']= [self.yzoom]* (xgrid.size-1)

		#if hasattr(self, pfm):
	
	def set_bins(self, xbins=[[1,10]], ybins=[[1,10]],xgrid=None, ygrid=None,Grid=False):
		'''
		Initialize period-radius (or mass) bins for occurrence rate calculations
		
		Description:
			Bins can be generated from a grid, f.e. xgrid=[1,10,100], 
			or from a list of bin edges, f.e. xbins= [[1,10],[10,100]]
		
		Args:
			xbins(list):	(list of) period bin edges
			ybins(list):	(list of) radius/mass bin edges
			xgrid(list):	period bin in interfaces
			ygrid(list):	radius/mas bin interfaces
			Grid(bool):
				If true, create a 2D grid from bins: nbins = nx ``*`` ny.
				If false, pair x and y bins: nbins == nx == ny
		'''
		if not hasattr(self,'occurrence'):
			self.occurrence={}
		focc= self.occurrence
		
		#generate list of bin inner and outer edges
		if xgrid is None:
			if np.ndim(xbins) ==1:
				assert len(xbins)==2
				_xbins= [xbins]
			elif np.ndim(xbins) ==2: 
				_xbins= xbins
			else:
				raise ValueError('wrong bin dimensions')
		else:
			_xbins= np.array([[i,j] for i,j in zip(xgrid[:-1],xgrid[1:])])

		if ygrid is None:
			if np.ndim(ybins) ==1:
				assert len(ybins)==2
				_ybins= [ybins]
			elif np.ndim(ybins) ==2: 
				_ybins= ybins
			else:
				raise ValueError('wrong bin dimensions')
		else:
			_ybins= np.array([[i,j] for i,j in zip(ygrid[:-1],ygrid[1:])])
		
		focc['bin']={}
		if Grid:
			focc['bin']['x']= np.tile(_xbins, (len(_ybins),1))
			#focc['bin']['y in']= np.tile(_ybins, (len(_xbins),1))
			focc['bin']['y in']= np.repeat(_ybins, len(_xbins), axis=0)
			# TODO: store 1d -> 2d mapping
			
		else:		
			if np.shape(_xbins) != np.shape(_ybins):
				raise ValueError('unequal amount of bins. Use Grid=True?')
			focc['bin']['x']= _xbins
			focc['bin']['y in']= _ybins

	def set_bins_poly(self, polys, labels=None):
		'''
		Initialize polygonic bins for occurrence rate calculations
		
		Description:
			each polygone is a Nx2 numpy array of (x,y) coordinates
		
		Args:
			polys(list):	(list of) polygones
		'''
		if not hasattr(self,'occurrence'):
			self.occurrence={}
		focc= self.occurrence
		
		focc['poly']={}
		focc['poly']['coords']=[]

		print '\nTrying {} polygons'.format(len(polys))
		for coords in polys:
			npc= np.asarray(coords)
			assert npc.ndim==2, 'coords needs to de 2dim list'
			assert npc.shape[-1]==2, 'coords need to be 2xN array'
			focc['poly']['coords'].append(npc)
	
		if labels is not None:
			assert len(labels) == len(focc['poly']['coords']), 'Names don\'t match up'
			focc['poly']['labels']= np.asarray(labels)

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
		if not callable(func): raise ValueError('func is not a callable function')		
		self.func=func
		
		self.pdfpars= fitparameters()
		
		self.Parametric= True
		self.PDF=True	
		self.fitpars=self.pdfpars
		self.summarystatistic= ['N','xvar','yvar']
	
	def set_multi(self, spacing=None):
		if not self.Parametric:
			raise ValueError('Define a parametric planet population first')
		self.Multi=True
		
		self.RandomPairing= (spacing==None)
		self.spacing= spacing # None, brokenpowerlaw, dimensionless

		# skipping yvar here
		self.summarystatistic= ['N','xvar','Nk','dP','Pin']
	
	def set_population(self, name, sma, mass, 
		radius=None, inc=None, starID=None, ecc=None, tag=None, 
		Verbose=False, **kwargs):
		# tag is fit parameter, i.e. metallicity, surface density, or model #

		if Verbose:
			print '\nArguments not used by set_population:'
			print '  ',kwargs.keys()

		if hasattr(self, 'pfm'):
			raise ValueError('expand: adding multiple populations?')
		
		self.Parametric=False
		self.modelpars= fitparameters()
		self.fitpars= self.modelpars
		
		if inc is None:
			self.summarystatistic= ['N'] # xvar, yvar?
		else:
			self.summarystatistic= ['N','Nk','dP'] #,'Pin']
		
		# length checks
		try:
			if len(sma) != len(mass): 
				raise ValueError('sma ({}) and mass ({}) not same length'.format(len(sma),len(mass)))
		except: 
			raise ValueError('sma ({}) and mass ({}) have to be iterable'.format(type(sma), type(mass)))
		
		# model has mutual inclinations?
		self.Multi= (inc != None)
		
		pfm= self.pfm= {}
		pfm['name']= name
		
		# lexsort?
		pfm['sma']= np.asarray(sma)
		pfm['M']= np.asarray(mass)
		
		# ID is array 0...ns with np elements
		if starID is None:
			pfm['ID']= np.arange(len(sma))
		else:
			_, pfm['ID']= np.unique(starID, return_inverse=True)
			
		if radius is not None:
			pfm['R']= np.asarray(radius)
		if ecc is not None:
			pfm['ecc']= np.asarray(ecc)			
		if tag is not None:
			pfm['tag']= np.asarray(tag)		
		if inc is None:
			self.Multi=False
		else:
			self.Multi= True
			pfm['inc']= np.asarray(inc)
				
		pfm['P']= pfm['sma']**1.5 * 365.25 # update
		
		pfm['np']= pfm['ID'].size
		pfm['ns']= np.unique(pfm['ID']).size
		
		''' If multiple planets per stars: Lexsort, period ratio'''
		if pfm['np'] > pfm['ns']:
			
			order= np.lexsort((pfm['sma'],pfm['ID'])) # sort by ID, then sma
			for key in ['ID','sma','M','P','inc','tag']:
				if key in pfm: pfm[key]=pfm[key][order]
			
			EPOS.multi.indices(pfm['ID'], Verbose=True)
			pfm['k'], pfm['Nk']= EPOS.multi.frequency(pfm['ID'], Verbose=True)
			
			# period ratio, multi-planet index
			single, multi, ksys, multis= EPOS.multi.nth_planet(pfm['ID'],pfm['P'])
			pfm['dP']=np.ones_like(pfm['P'])
			pfm['kth']=np.zeros_like(pfm['P'], dtype=int)

			pfm['dP'][single]= 0 #np.nan
			for k, km in enumerate(multis[1:]):
				# 2nd, 3rd, 4th??
				pfm['dP'][km]= pfm['P'][km]/pfm['P'][np.array(km)-1]
				pfm['kth'][km]= k+1
			#print pfm['ID'][1:6]
			#print pfm['P'][1:6]
			#print pfm['dP'] # not ok?	
			#for a,b in zip(pfm['ID'], pfm['kth']):
			#	print a,b

			# innermost planet in multi
			pfm['Pin']= np.squeeze(pfm['P'][np.array(multis[1])-1])

			pfm['system index']= np.arange(pfm['ns']) # 0...ns with ns elements
			pfm['planet index']= np.arange(pfm['np']) # 0...np with np elements

			if 'tag' in pfm:
				_, index= np.unique(starID, return_index=True)
				pfm['system tag']= pfm['tag'][index]
			#pfm['tag'] == pfm['system tag'][pfm['ID']]
		else:
			pfm['system tag']= pfm['tag']

		pfm['M limits']=[np.min(pfm['M']),np.max(pfm['M'])]
		pfm['P limits']=[np.min(pfm['P']),np.max(pfm['P'])]
		
		# set plot limits in model 5% wider than data
		xmin,xmax= min(pfm['P']), max(pfm['P'])
		ymin,ymax= min(pfm['M']), max(pfm['M'])
		dx= (xmax/xmin)**0.05
		dy= (ymax/ymin)**0.05
		self.mod_xlim=[xmin/dx, xmax*dx]
		self.mod_ylim=[ymin/dy, ymax*dy]
		# ??
		self.mod_xlim=[min(xmin/dx, self.mod_xlim[0]),max(xmax*dx, self.mod_xlim[1])]
		self.mod_ylim=[min(ymin/dy, self.mod_ylim[0]),max(ymax*dy, self.mod_ylim[1])]
	
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
