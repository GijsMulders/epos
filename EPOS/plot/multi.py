import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import scipy.stats # norm

import helpers
import EPOS.multi
from EPOS.population import periodradius as draw_PR

# backwards compatible colors 
import matplotlib
if matplotlib.__version__[0] != 2: 
	helpers.default_pyplot2_colors(matplotlib.colors)

def polar(epos):
	'''
	Plot plane populations as a half circle (not quite working) 
	'''
	# plot multiplicity
	f, axlist = plt.subplots(2,2, subplot_kw=dict(projection='polar'))
	
	pop=epos.population

	for ax, key in zip([axlist[0,0],axlist[0,1],axlist[1,0]],['system','single','multi']):
		ax.set_title(key)
		ax.set_xlabel('fraction')
		ax.set_ylabel('P [days]')

		ax.set_xlim(-0.05, 1.05)
		ax.set_ylim(0.1, 1000) # 7.
		ax.set_yscale('log')
		
		ax.set_rgrids([1, 10, 100])

		try:
			ax4.set_thetamin(0)
			ax4.set_thetamax(90)
		except AttributeError:
			print 'Update pyplot'
			raise
		
		ax.plot(pop[key]['order']*np.pi, pop[key]['P'],
				ls='', marker='.', mew=0, ms=3, color='k')

	ax4=axlist[1,1]
	ax4.set_title('all')
	ax4.set_xlabel('fraction')
	ax4.set_ylabel('P [days]')

	ax4.set_xlim(-0.05, 1.05)
	ax4.set_ylim(0.1, 1000) # 7.
	ax4.set_yscale('log')
	
	try:
		ax4.set_thetamin(0)
		ax4.set_thetamax(90)
	except AttributeError:
		print 'Update pyplot'
		raise
	
	ax4.plot(pop['order']*np.pi, pop['P'],
		ls='', marker='.', mew=0, ms=3, color='gray')

	helpers.save(plt, '{}/polar_test'.format(epos.plotdir))

def periodradius(epos, Nth=False, MC=True):

	f, (ax, axR, axP)= helpers.make_panels(plt)
	
	if MC:
		sim=epos.synthetic_survey
		ID= sim['ID']
		P= sim['P']
		Y= sim['Y']
		outdir='output'
		title= 'detectable planets'
	else:
		ID= epos.obs_starID
		P= epos.obs_xvar
		Y= epos.obs_yvar
		outdir='survey'
		title= 'detected planets'
	
	''' plot R(P), main panel'''
	ax.set_title(title)
	helpers.set_axes(ax, epos, Trim=True)
		
	''' Period side panel '''
	helpers.set_axis_distance(axP, epos, Trim=True)
	#axP.set_yscale('log')
	#axP.set_ylim([2e-3,5])	
	#axP.set_yticks([0.01,0.1,1])
	#axP.set_yticklabels(['1%','10%','100%'])
	axP.yaxis.tick_right()
	axP.yaxis.set_ticks_position('both')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')
	
	#axP.hist(sim['P'], bins=epos.MC_xvar, color='0.7')

	''' Radius side panel'''
	helpers.set_axis_size(axR, epos, Trim=True, In= epos.MassRadius)

	#axR.set_xscale('log')
	#axR.set_xlim([2e-3,5])
	#axR.set_xticks([1,10,100,1000])
	#axR.set_xticklabels(['1','10','100','1000'], rotation=70)
	for tick in axR.get_xticklabels():
		tick.set_rotation(70)
	#axR.tick_params(axis='x', which='minor',top='off',bottom='off')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')

	''' which multiplanets to color '''
# 	single, multi= EPOS.multi.indices(sim['ID'])
# 	for k, (label, subset) in enumerate(zip(['single','multi'],[single, multi])):
# 		ax.plot(sim['P'][subset], sim['Y'][subset], ls='', marker='.', mew=0, ms=5.0, \
# 			label=label)

	if Nth:
		single, multi, ksys, multis= EPOS.multi.nth_planet(ID, P)
		suffix= '.nth'
		label_single='single'
	else:
		single, multi, ksys, multis= EPOS.multi.indices(ID)
		suffix= ''
		label_single='1'
	
	ax.plot(P[single], Y[single], ls='', marker='.', \
			color='0.7', label=label_single)
			
	for k, subset in zip(ksys, multis):
		ht, = ax.plot(P[subset], Y[subset], ls='', marker='.', \
			label=k)
		
		if True:
			# pdf
			axP.hist(P[subset], bins=epos.MC_xvar, 
				color=ht.get_color(), histtype='step') #, cumulative=True, normed=1)
			if k==ksys[-1]:
				axP.hist(P[single], bins=epos.MC_xvar, color='0.7', histtype='step') 	
		else:
			# cumulative
			Plist= np.sort(P[subset])
			axP.step(Plist,np.arange(Plist.size,dtype=float)/Plist.size )

	#axR.hist(Y,orientation='horizontal', bins=epos.MC_yvar, color='0.7')
	axR.hist(Y,orientation='horizontal', bins=epos.MC_yvar, color='k',histtype='step')
	axR.hist(Y[single],orientation='horizontal', bins=epos.MC_yvar, color='0.7')

	
	#ax.legend(loc='lower left', shadow=False, prop={'size':14}, numpoints=1)
	ax.legend(bbox_to_anchor=(1.0, 1.0), markerscale=3)

	helpers.save(plt, '{}{}/PR.multi{}'.format(epos.plotdir, outdir, suffix))

def multiplicity(epos, MC=False, Planets=False, MCMC=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('planet multiplicity')
	ax.set_xlabel('planets per system')
	ax.set_ylabel('number of planets' if Planets else 'number of systems')

	ax.set_xlim(0.5, 7.5)
	if not Planets:
		ax.set_ylim(0.5, 1e4) # 7.	
		ax.set_yscale('log')
		#ax.get_yaxis().set_major_formatter(tck.ScalarFormatter())

	key = 'pl cnt' if Planets else 'count'

	# MC data
	if MC or MCMC:
		ss=epos.synthetic_survey
		if Planets:
			ax.bar(ss['multi']['bin'], ss['multi'][key], color='C0',label=epos.name, width=1)

			f_iso= epos.fitpars.get('f_iso')
			if f_iso > 0:
				fsingle= np.sum(ss['multi']['count'])*f_iso
				ax.bar(1, fsingle, bottom= ss['multi'][key][0]-fsingle, 
					color='',label='dichotomy', width=1, hatch='xx') #, ec='k')
				#ax.plot(1, ss['multi'][key][0]-fsingle, marker='+', ms=10, ls='', color='k', label='no dichotomy')
		elif MCMC:
			for ss in epos.ss_sample:
				ax.hlines(ss['multi'][key], ss['multi']['bin']-0.5,ss['multi']['bin']+0.5,
					color='b', alpha=0.1)
				
		else:
			#ax.step(ss['multi']['bin'], ss['multi'][key], color='C0',label=epos.name, where='mid')
			ax.plot(ss['multi']['bin'], ss['multi'][key], 
				ls='', marker='+', ms=10, color='C0',label=epos.name)
		
		if hasattr(epos, 'ss_extra'):
			# advance color cycle
			# ax._get_lines.get_next_color()
			ax.plot([], [])
			ax.plot([], [])

			# plot extra epos runs
			for ss in epos.ss_extra:
				ax.plot(ss['multi']['bin'], ss['multi'][key], 
					ls='', marker='+', ms=10, label=ss['name'])
	
		# observations in same region 
		ax.step(epos.obs_zoom['multi']['bin'], epos.obs_zoom['multi'][key], 
			where='mid', color='C1', label='Kepler subset')

	else:
		# observations
		ax.plot(epos.multi['bin'], epos.multi[key], drawstyle='steps-mid', 
			ls='--', marker='', color='gray', label='Kepler all')

	ax.legend(loc='upper right', shadow=False, prop={'size':14}, numpoints=1)
	
	prefix= 'output' if MC else 'survey'
	if MCMC: prefix= 'mcmc'
	suffix='.planets' if Planets else ''
	
	helpers.save(plt, '{}{}/multiplicity{}'.format(epos.plotdir,prefix,suffix))

def multiplicity_cdf(epos, MC=False):
	# plot multiplicity cdf
	f, ax = plt.subplots()
	ax.set_title('planet multiplicity')
	ax.set_xlabel('planets per system')
	ax.set_ylabel('cumulative number of planets')

	ax.set_xlim(0, 9)
	ax.set_ylim(-0.01,1.05)
	#ax.set_yscale('log')

	# MC data
	if MC:
		ss=epos.synthetic_survey
		bincounts= np.sort(ss['multi']['cdf'])
		cdf= np.arange(bincounts.size, dtype=float)/bincounts.size
		ax.plot(bincounts, cdf, 
			drawstyle='steps-mid', 
			ls='-', marker='', mew=2, ms=10, color='b',label=epos.name)
	
		# observations in same region 
		bincounts= np.sort(epos.obs_zoom['multi']['cdf'])
		cdf= np.arange(bincounts.size, dtype=float)/bincounts.size
		ax.plot(bincounts, cdf,  
			drawstyle='steps-mid', 
			ls='-', marker='', color='k', label='Kepler subset')

	# observations
	#ax.plot(epos.multi['bin'], epos.multi['count'], drawstyle='steps-mid', 
	#	ls='--', marker='', color='gray', label='Kepler all')
		
	ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
	
	prefix= 'output' if MC else 'survey'
	
	helpers.save(plt, '{}{}/cdf'.format(epos.plotdir,prefix))
	
def periodratio(epos, MC=False, N=False, Input=False, MCMC=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('period ratio adjacent planets')
	ax.set_xlabel('period outer/inner')
	ax.set_ylabel('PDF')

	ax.set_xlim(1, 10)
	ax.set_xscale('log')
	for s in [ax.set_xticks,ax.set_xticklabels]: s([1,2,3,4,5,7,10])
	ax.set_xticks([], minor=True) # minor ticks generate labels

	bins=np.logspace(0,1, 15) # 1,10
	
	# MC data
	if MC or MCMC:
		ss=epos.synthetic_survey
		if MCMC:
			for ss in epos.ss_sample:
				# bar?
				ax.hist(ss['multi']['Pratio'], bins=bins, histtype='step', color='b', alpha=0.1)
				#ax.hist(ss['multi']['Pratio'], bins=bins, color='b', alpha=1./len(epos.ss_sample))
		else:
			ax.hist(ss['multi']['Pratio'], bins=bins, 
					color='C0', histtype='step', label=epos.name)
	
		if N:
			# planets inbetween?
			#ax.hist(ss['multi']['dPN'][0], \
			#	bins=bins, ec='k', histtype='step', label='Adjacent planet')
			ax.hist(np.concatenate(ss['multi']['dPN'][1:]), 
				bins=bins, hatch='xx',
				histtype='stepfilled', label='Planet inbetween') # ec, color
		else:
			#ax.axvline(np.median(ss['multi']['Pratio']), color='C0', ls='--')
			pass
		
		if hasattr(epos, 'ss_extra'):
			# advance color cycle
			# ax._get_lines.get_next_color()
			ax.plot([], [])
			ax.plot([], [])
			ax.plot([], [])
			
			# plot extra epos runs
			for ss in epos.ss_extra:
				ax.hist(ss['multi']['Pratio'], bins=bins, 
					histtype='step', label=ss['name'])
		
	else:
		# observed all
		ax.hist(epos.multi['Pratio'], bins=bins, color='0.7', label='Kepler all')
		#ax.axvline(np.median(epos.multi['Pratio']), color='0.7', ls='--')
	
	''' input distribution '''	
	if Input and epos.populationtype=='parametric':
		#Pgrid= np.logspace(0,1)
		Pgrid=bins
		if epos.spacing == 'powerlaw':
			dPbreak= epos.fitpars.get('dP break')
			dP1= epos.fitpars.get('dP 1')
			dP2= epos.fitpars.get('dP 2')
			pdf= EPOS.fitfunctions.brokenpowerlaw1D(Pgrid, dPbreak, dP1, dP2)
		elif epos.spacing=='dimensionless':
			logD=  epos.fitpars.get('log D')
			sigma= epos.fitpars.get('sigma')	

			with np.errstate(divide='ignore'): 
				Dgrid= np.log10(2.*(Pgrid**(2./3.)-1.)/(Pgrid**(2./3.)+1.))
			Dgrid[0]= -2
			#print Dgrid
			pdf= scipy.stats.norm(logD,sigma).pdf(Dgrid)
		pdf*= 0.95* ax.get_ylim()[1] / max(pdf)
		
		ax.plot(Pgrid, pdf, marker='', ls='-', color='r',label='input')

	elif epos.Zoom:
		# Observed zoom
		ax.hist(epos.obs_zoom['multi']['Pratio'], \
			bins=bins, ec='C1', histtype='step', label='Kepler subset')
		if not MCMC:
			ax.axvline(np.median(epos.obs_zoom['multi']['Pratio']), color='C1', ls='--')
	
	prefix= 'output' if MC else 'survey'
	suffix= '.index' if N else ''
	if Input: suffix+= '.input' 
	if MCMC: prefix= 'mcmc'

	if MC: ax.legend(loc='upper right', shadow=False, prop={'size':14}, numpoints=1)
		
	helpers.save(plt, '{}{}/periodratio{}'.format(epos.plotdir, prefix, suffix))	

def periodratio_cdf(epos, MC=False):

	# plot multiplicity CDF
	f, ax = plt.subplots()
	ax.set_title('period ratio adjacent planets')
	ax.set_xlabel('period outer/inner')
	ax.set_ylabel('CDF')

	ax.set_xlim(1, 10)
	ax.set_xscale('log')
	for s in [ax.set_xticks,ax.set_xticklabels]: s([1,2,3,4,5,7,10])

	# MC data
	if MC:
		ss=epos.synthetic_survey
		Psort= np.sort(ss['multi']['Pratio'])
		cdf= np.arange(Psort.size, dtype=float)/Psort.size
		ax.plot(Psort, cdf, color='b', label=epos.name)
	
		# Observed zoom
		Psort= np.sort(epos.obs_zoom['multi']['Pratio'])
		cdf= np.arange(Psort.size, dtype=float)/Psort.size
		ax.plot(Psort, cdf, color='k', label='Kepler subset')

	else:
		for resonance in [2./1., 3./2.]: ax.axvline(resonance, ls=':', color='g')

	# all observed
	Psort= np.sort(epos.multi['Pratio'])
	cdf= cdf= np.arange(Psort.size, dtype=float)/Psort.size
	ax.plot(Psort, cdf ,color='gray' if MC else 'k', label='Kepler all')	
	
	prefix= 'output' if MC else 'survey'
	
	if MC: ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
		
	helpers.save(plt, '{}{}/periodratio.cdf'.format(epos.plotdir,prefix))	

''' these are practically identical to periodratio -> merge?'''

def periodinner(epos, MC=False, N=False, Input=False, MCMC=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('period innermost planet')
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('PDF')

	ax.set_xscale('log')
	ax.set_xlim(epos.xtrim)
		
	#bins= np.geomspace(*epos.xtrim)
	bins= epos.MC_xvar

	#pdf= regression.sliding_window_log(P, None, xgrid) #, width=2. )
	#ax3.plot(xgrid, pdf, ls='-', marker='', color=clrs[k % 4],label='{} x{:.3f}'.format(sg['name'], sg['weight']))
				
	# MC data
	if MC or MCMC:
		ss=epos.synthetic_survey
		if MCMC:
			for ss in epos.ss_sample:
				ax.hist(ss['multi']['Pinner'], bins=bins, 
						color='b', alpha=0.1, histtype='step')			
		else:
			ax.hist(ss['multi']['Pinner'], bins=bins, 
					color='C0', histtype='stepfilled', label=epos.name)
	
		if N:
			# Innermost is nth planet
			#ax.hist(ss['multi']['PN'][0], bins=bins, 
			#		ec='k', histtype='step', label='actual inner planet')
			# Not actual inner planet
			ax.hist(np.concatenate(ss['multi']['PN'][1:]), bins=bins, 
					color='r', histtype='stepfilled', label='not inner planet')

	else:
		# observed all
		ax.hist(epos.multi['Pinner'], bins=bins, color='0.7', label='Kepler all')

	''' Initial distribution or zoomed observations '''
	if Input and epos.populationtype=='parametric':
		_, _, pdf0_X, _= draw_PR(epos, Init=True, ybin=epos.yzoom)
		norm= 0.95* ax.get_ylim()[1] / max(pdf0_X)
		ax.plot(epos.MC_xvar, pdf0_X*norm, marker='',ls='-',color='r',label='input')
				
	elif epos.Zoom and not N:
		ax.hist(epos.obs_zoom['multi']['Pinner'], bins=bins, 
			ec='C1', histtype='step', label='Kepler subset')
	
	prefix= 'output' if MC else 'survey'
	suffix= '.index' if N else '' 
	suffix= '.input' if Input else suffix
	if MCMC: prefix= 'mcmc'

	if MC: ax.legend(loc='upper right', shadow=False, prop={'size':14}, numpoints=1)
		
	helpers.save(plt, '{}{}/innerperiod{}'.format(epos.plotdir, prefix,suffix))	

def periodinner_cdf(epos, MC=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('period innermost planet')
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('CDF')

	ax.set_xscale('log')
	ax.set_xlim(epos.xtrim)

	ax.set_ylim(-0.05, 1.05)
		
	bins= np.logspace(*np.log10(epos.xtrim))

	#pdf= regression.sliding_window_log(P, None, xgrid) #, width=2. )
	#ax3.plot(xgrid, pdf, ls='-', marker='', color=clrs[k % 4],label='{} x{:.3f}'.format(sg['name'], sg['weight']))
				
	# MC data
	if MC:
		ss=epos.synthetic_survey
		Psort= np.sort(ss['multi']['Pinner'])
		cdf= np.arange(Psort.size, dtype=float)/Psort.size
		ax.plot(Psort, cdf, color='b', label=epos.name)
	
		# Observed zoom
		Psort= np.sort(epos.obs_zoom['multi']['Pinner'])
		cdf= np.arange(Psort.size, dtype=float)/Psort.size
		ax.plot(Psort, cdf, color='k', label='Kepler subset')
	else:
		# observed all
		Psort= np.sort(epos.multi['Pinner'])
		cdf= cdf= np.arange(Psort.size, dtype=float)/Psort.size
		ax.plot(Psort, cdf ,color='gray' if MC else 'k', label='Kepler all')	
	
	prefix= 'output' if MC else 'survey'

	if MC: ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
		
	helpers.save(plt, '{}{}/innerperiod.cdf'.format(epos.plotdir, prefix))	
