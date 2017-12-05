import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import helpers


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

def multiplicity(epos, MC=False, Planets=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('planet multiplicity')
	ax.set_xlabel('planets per system')
	ax.set_ylabel('number of planets' if Planets else 'number of systems')

	ax.set_xlim(0, 9)
	if not Planets:
		ax.set_ylim(0.5, 1e4) # 7.	
		ax.set_yscale('log')
		#ax.get_yaxis().set_major_formatter(tck.ScalarFormatter())

	key = 'pl cnt' if Planets else 'count'

	# MC data
	if MC:
		ss=epos.synthetic_survey
		ax.plot(ss['multi']['bin'], ss['multi'][key], 
			ls='', marker='+', mew=2, ms=10, color='k',label=epos.name)
		
		if hasattr(epos, 'ss_extra'):
			for ss in epos.ss_extra:
				ax.plot(ss['multi']['bin'], ss['multi'][key], 
					ls='', marker='+', mew=2, ms=10, color='g',label='no dichotomy')
	
		# observations in same region 
		ax.plot(epos.obs_zoom['multi']['bin'], epos.obs_zoom['multi'][key], 
			drawstyle='steps-mid', 
			ls='-', marker='', color='k', label='Kepler subset')

	# observations
	ax.plot(epos.multi['bin'], epos.multi[key], drawstyle='steps-mid', 
		ls='--', marker='', color='gray', label='Kepler all')

	ax.legend(loc='upper right', shadow=False, prop={'size':14}, numpoints=1)
	
	prefix= 'output' if MC else 'survey'
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
	
def periodratio(epos, MC=False, N=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('period ratio adjacent planets')
	ax.set_xlabel('period outer/inner')
	ax.set_ylabel('PDF')

	ax.set_xlim(1, 10)
	ax.set_xscale('log')
	for s in [ax.set_xticks,ax.set_xticklabels]: s([1,2,3,4,5,7,10])
		
	#bins=np.linspace(1,10, 9*5+1)
	bins=np.logspace(0,1, 15)
	
	# MC data
	if MC:
		ss=epos.synthetic_survey
		ax.hist(ss['multi']['Pratio'], bins=bins, 
				ec='b', histtype='step', label=epos.name)
	
		if not N:
			ax.axvline(np.median(ss['multi']['Pratio']), color='b', ls='--')

			# Observed zoom
			ax.hist(epos.obs_zoom['multi']['Pratio'], \
				bins=bins, ec='k', histtype='step', label='Kepler subset')
			ax.axvline(np.median(epos.obs_zoom['multi']['Pratio']), color='k', ls='--')
		else:
			# planets inbetween?
			#ax.hist(ss['multi']['dPN'][0], \
			#	bins=bins, ec='k', histtype='step', label='Adjacent planet')
			ax.hist(np.concatenate(ss['multi']['dPN'][1:]), \
				bins=bins, ec='k', histtype='stepfilled', label='Planet inbetween')
			
			if epos.populationtype=='parametric' and epos.Multi and \
				'dP break' in epos.fitpars.keysall:
				
				from EPOS.fitfunctions import brokenpowerlaw1D
				pbreak= epos.fitpars.get('dP break')
				p1= epos.fitpars.get('dP 1')
				p2= epos.fitpars.get('dP 2')
			
				xx= np.logspace(0,1)
				yy= 150.*brokenpowerlaw1D(xx, pbreak, p1, p2)

				ax.plot(xx, yy, marker='', ls=':', color='r',label='input')
		
	else:
		# observed all
		ax.hist(epos.multi['Pratio'], bins=bins, color='b', label='Kepler all')
		ax.axvline(np.median(epos.multi['Pratio']), color='b', ls='--')

		# quick fit to Pratio distributio: what form is this??		
		from scipy import stats
		xx= np.logspace(0,1)
		#yy= 40.* stats.norm.pdf(np.log10(xx-1), np.log10(0.5), 0.1)
		#yy= 50.* stats.norm.pdf(np.log10(xx-0.3), np.log10(1.5), 0.12)
		#yy= 50.* stats.rayleigh.pdf(np.log10(xx), scale=np.log10(2.0)) # too wide
		#yy= 150.* (xx/2)**(-2.5) * stats.norm.cdf(xx, 1.6, 0.2)  #np.exp(1.-1./(xx/2.))
		
		# skewed norm is reasonable but can't draw from it.
		#frozen= stats.norm(loc=np.log10(1.7), scale=0.12)
		#skewnorm= 2.*frozen.pdf(np.log10(xx))*frozen.cdf(2.0*np.log10(xx))
		#yy=30. *skewnorm
		#yy= 150. * stats.skewnorm.pdf(np.log10(xx), 0.0, np.log10(1.5), 0.12)
		#ax.plot(xx, yy, marker='', ls=':', color='r')
		
		if epos.populationtype=='parametric' and epos.Multi and \
				'dP break' in epos.fitpars.keysall:
			from EPOS.fitfunctions import brokenpowerlaw1D
			pbreak= epos.fitpars.get('dP break',Init=True)
			p1= epos.fitpars.get('dP 1',Init=True)
			p2= epos.fitpars.get('dP 2',Init=True)
			yy= 170.*brokenpowerlaw1D(xx, pbreak, p1, p2)
			ax.plot(xx, yy, marker='', ls=':', color='r')
	
	prefix= 'output' if MC else 'survey'
	suffix= '.index' if N else '' 

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

def periodinner(epos, MC=False, N=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('period innermost planet')
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('PDF')

	ax.set_xscale('log')
	ax.set_xlim(epos.xtrim)
		
	bins= np.logspace(*np.log10(epos.xtrim))

	#pdf= regression.sliding_window_log(P, None, xgrid) #, width=2. )
	#ax3.plot(xgrid, pdf, ls='-', marker='', color=clrs[k % 4],label='{} x{:.3f}'.format(sg['name'], sg['weight']))
				
	# MC data
	if MC:
		ss=epos.synthetic_survey
		ax.hist(ss['multi']['Pinner'], bins=bins, 
				ec='b', histtype='step', label=epos.name)
	
		if not N:
			# Observed zoom
			ax.hist(epos.obs_zoom['multi']['Pinner'], bins=bins, 
				ec='k', histtype='step', label='Kepler subset')
		else:
			# Innermost is nth planet
			ax.hist(ss['multi']['PN'][0], bins=bins, 
					ec='k', histtype='step', label='actual inner planet')
			# Not actual inner planet
			ax.hist(np.concatenate(ss['multi']['PN'][1:]), bins=bins, 
					ec='k', histtype='stepfilled', label='not inner planet')
			
	else:
		# observed all
		ax.hist(epos.multi['Pinner'], bins=bins, color='b', label='Kepler all')
	
	# Input distribution (not very useful)
# 	if epos.populationtype=='parametric' and epos.Multi:
# 		try:
# 			pdf= epos.func(epos.X, epos.Y, *epos.pfit[epos.ip2d])
# 		except:
# 			raise
# 			pdf= epos.func(epos.X, epos.Y, *epos.p0[epos.ip2d])
# 		pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
# 		ax.plot(epos.MC_xvar, 500.*pdf_X, marker='', ls=':', color='r',label='input')
	
	prefix= 'output' if MC else 'survey'
	suffix= '.index' if N else '' 

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
