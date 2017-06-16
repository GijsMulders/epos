import numpy as np
import cgs
from EPOS import multi
from EPOS.plot.helpers import set_pyplot_defaults

def readme():
	print '\nThis module defines the EPOS class, that contains the observed exoplanets,'
	print 'the survey detection efficiency, and the synthetic planet population.'
	print 
	print 'see example.py for a simple demonstration of the class'
	print 
	print 'epos(name, RV=False)'
	print '  Initialize the class'
	print 'in:'
	print '  name: identifier, plots will appear in png/name'
	print '  RV: transit or RV? [T/F]'
	print 
	print 'set_observation(xvar, yvar, starID, nstars=1.6862e5)'
	print '  Observed planet population'
	print 'in:'
	print '  xvar: planet orbital period [list]'
	print '  yvar: planet radius or M sin i [list]'
	print '  ID: planet ID [list]'
	print '  nstars: number of stars in the survey'
	print
	print 'set_survey(xvar, yvar, eff_2D, Rstar=1.0)'
	print '  Survey detection efficiency (completeness)'
	print 'in:'
	print '  xvar: planet orbital period grid [list]'
	print '  yvar: planet radius or M sin i grid [list]'
	print '  eff_2D: 2D matrix of detection efficiency'
	print '  Rstar: stellar radius for calculating transit probability'

class epos:
	
	def __init__(self, name, RV=False, Debug=False):
		self.name=name
		self.plotdir='png/{}/'.format(name)
		self.RV= RV

		# switches to be set later		
		self.Observation=False
		self.Range=False
		self.DetectionEfficiency=False
		self.Occurrence= False # inverse detection efficiency (?)
		self.Prep= False # ready to run? EPOS.run.once()
		self.RadiusMassConversion= False
		self.Radius= False
		self.Isotropic= False
		
		self.populationtype=None # ['parametric','model']
				
		self.Debug= False
		set_pyplot_defaults() # nicer plots
		
	def set_observation(self, xvar, yvar, starID, nstars=1.6862e5):
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
		print '\nObservations:\n  {:.0f} stars'.format(nstars)
		print '  {} planets'.format(self.obs_starID.size)
		multi.indices(self.obs_starID, Verbose=True)
		epos.multi={}
		epos.multi['bin'], epos.multi['count']= \
			multi.frequency(self.obs_starID, Verbose=True)
		epos.multi['Pratio']= \
			multi.periodratio(self.obs_starID, self.obs_xvar, Verbose=True)
		epos.multi['cdf']= multi.cdf(self.obs_starID, Verbose=True)	
		
	def set_survey(self, xvar, yvar, eff_2D, Rstar=1.0):
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
		
		if self.RV:
			self.completeness= self.eff_2D
		else:
			self.Rstar=Rstar # Solar radii
			self.Pindex= -2./3.
			self.fgeo_prefac= self.Rstar*cgs.Rsun/ (cgs.au /365.24**(2./3.))
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
		
		''' Prep the grid '''
		# make sure range _encompasses_ trim
		ixmin,ixmax= _trimarray(self.eff_xvar, self.xtrim)
		iymin,iymax= _trimarray(self.eff_yvar, self.ytrim)
	
		self.MC_xvar= self.eff_xvar[ixmin:ixmax]
		self.MC_yvar= self.eff_yvar[iymin:iymax]
		self.MC_eff= self.eff_2D[ixmin:ixmax,iymin:iymax]
	
		self.X, self.Y= np.meshgrid(self.MC_xvar, self.MC_yvar, indexing='ij')
		
		#print '\nTrimming {} to {}'.format(self.eff_2D.shape,self.eff_trim.shape) 
		#print 'xlim: {}-{}'.format(*self.xtrim)
		#print 'ylim: {}-{}'.format(*self.ytrim)
		#print '\nTrim: {}'.format(self.eff_xvar)
		#print 'To  : {}'.format(self.MC_xvar)
		#print '\nTrim: {}'.format(self.eff_yvar)
		#print 'To  : {}'.format(self.MC_yvar)

		self.Range=True
	
	def set_parametric(self, func=None, p0=[], pname=None):
		if self.populationtype is None:
			self.populationtype='parametric'
		elif self.populationtype is not 'parametric':
			raise ValueError('You have already defined a planet population ({})'.format(self.populationtype))
		
		if not callable(func): raise ValueError('func is not a callable function')
		if not len(p0)>0: raise ValueError('nonzero list of starting params')
		
		# Call function once to see if it works, P=1, R=1
		try: 	func(np.asarray([1.]), np.asarray([1.]), *p0)
		except:	raise
		
		self.func=func
		self.p0= p0
		self.pname= ['c{}'.format(i) for i in range(len(p0))] if pname is None else pname
		
		self.Isotropic=True
		
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
			self.Isotropic=False # isotropic if any subgroup lacks inclinations
			self.Multi=False
		
		sg={}
		sg['name']= name
		sg['n']= len(sma)
		sg['system']=[] # remains empty if list
		sg['weight']= weight

		# list of lists or 1-dim list?
		if type(sma[0]) is list:
			# loop over systems
			for k, (L_sma, L_mass) in enumerate(zip(sma, mass)):
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
				print 'system has {} planets:'.format(system['np'])
				for a,m in zip(system['sma'], system['mass']):
					print '  {:.2f} au, {:.1f} Mearth'.format(a,m)
				print sg['all_Pratio']
				print
		
		# set plot limits in model 5% wider than data
		xmin,xmax= min(sg['all_sma']), max(sg['all_sma'])
		ymin,ymax= min(sg['all_mass']), max(sg['all_mass'])
		dx= (xmax/xmin)**0.05
		dy= (ymax/ymin)**0.05
		self.mod_xlim=[min(xmin/dx, self.mod_xlim[0]),max(xmax*dx, self.mod_xlim[1])]
		self.mod_ylim=[min(ymin/dy, self.mod_ylim[0]),max(ymax*dy, self.mod_ylim[1])]
		
		self.groups.append(sg)
	
	def set_massradius(self, RM, name):
		if self.RadiusMassConversion:
			raise ValueError('Already defined a Radius-Mass conversion function ')
		# actually radius as funtion of mass (mass-to-radius)
		if not callable(RM): raise ValueError('Radius-Mass function is not callable')

		self.RadiusMassConversion=True
		self.RM=RM
		self.RM_label=name

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