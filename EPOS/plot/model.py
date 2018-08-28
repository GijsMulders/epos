import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from matplotlib import gridspec

import helpers
from EPOS.population import periodradius
from EPOS.fitfunctions import brokenpowerlaw1D

# Kepler obs
# Kepler pdf mass
# kepler pdf radius
# Kepler pdf inner (M/R ??)

fmt_symbol= {'ls':'', 'marker':'o', 'mew':1, 'mec':'k', 'ms':8,'alpha':0.6}

def panels_mass(epos, Population=False, color='C1'):
	f, (ax, axM, axP)= helpers.make_panels(plt)
	pfm=epos.pfm
	eta= epos.modelpars.get('eta',Init=True)

	''' Bins '''
	dw= 0.5 # bin width in ln space
	xbins= np.exp(np.arange(np.log(pfm['P limits'][0]),np.log(pfm['P limits'][-1])+dw,dw))
	ybins= np.exp(np.arange(np.log(pfm['M limits'][0]),np.log(pfm['M limits'][-1])+dw,dw))

	''' Posterior '''
	if Population:
		assert hasattr(epos, 'func')
		fname='.pop'
		
		# side panels marginalized over M and P limits
		pps, pdf, pdf_X, pdf_Y= periodradius(epos, Init=True)
		_, _, pdf_X, _= periodradius(epos, Init=True, ybin=ybins)
		_, _, _, pdf_Y= periodradius(epos, Init=True, xbin=xbins)
		pps, _ , _, _ = periodradius(epos, Init=True, xbin=xbins, ybin=ybins)
		#pdf/= np.max(pdf)
		#pdflog= np.log10(pdf) # in %
		levels= np.linspace(0,np.max(pdf))
		lines= np.array([0.1, 0.5]) * np.max(pdf)
		
		ax.contourf(epos.X_in, epos.Y_in, pdf, cmap='Purples', levels=levels)
		#ax.contour(epos.X_in, epos.Y_in, pdf, levels=lines)
		
		# Side panels
		#print 'pps model= {}'.format(eta)
		scale=dw
		axP.plot(epos.MC_xvar, pdf_X*scale, marker='',ls='-',color='purple')
		axM.plot(pdf_Y*scale, epos.in_yvar, marker='',ls='-',color='purple')
	else:
		fname=''

	''' plot main panel'''
	ax.set_title(epos.name)
	#helpers.set_axes(ax, epos, Trim=True)

	ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim(epos.mod_xlim)
	#ax.set_ylim(epos.mod_ylim)
		
	ax.plot(pfm['P'],pfm['M'], color=color, **fmt_symbol)
	xlim= ax.get_xlim()
	ylim= ax.get_ylim()
	
	''' Period side panel '''
	#axP.yaxis.tick_right()
	#axP.yaxis.set_ticks_position('both')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')

	axP.set_xscale('log')
	axP.set_xlim(xlim)
	#ax.set_xlabel('Semi-Major Axis [au]')
	axP.set_xlabel('Orbital Period [days]')

	axP.hist(pfm['P'], color=color, bins=xbins, weights=np.full(pfm['np'], eta/pfm['np']) )

# 	try:
# 		xbins= np.geomspace(*xlim, num=20)
# 		ybins= np.geomspace(*ylim, num=10)
# 	except:
# 		xbins= np.logspace(*np.log10(xlim), num=20)
# 		ybins= np.logspace(*np.log10(ylim), num=10)	

	''' Mass side panel'''
	#helpers.set_axis_size(axR, epos, Trim=True) #, In= epos.MassRadius)
	axM.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
	axM.set_yscale('log')
	axM.set_ylim(ylim)

	axM.hist(pfm['M'], bins=ybins, orientation='horizontal', weights=np.full(pfm['np'], eta/pfm['np']), color=color )
	
	helpers.save(plt, '{}model/input.mass{}'.format(epos.plotdir, fname))
	
def panels_radius(epos, Population=False, Occurrence=False, Observation=False, Tag=False, color='C1'):
	f, (ax, axR, axP)= helpers.make_panels(plt)
	pfm=epos.pfm
	eta= epos.modelpars.get('eta',Init=True)
	
	if not 'R' in pfm:
		pfm['R'], _= epos.MR(pfm['M'])

	if Tag:
		# function that return a simulation subset based on the tag
		subset={'Fe/H<=0': lambda tag: tag<=0,
			'Fe/H>0': lambda tag: tag>0}

	''' Bins '''
	dwR=0.2 # bin width in ln space
	dwP=0.5
	xbins= np.exp(np.arange(np.log(epos.xzoom[0]),np.log(epos.xzoom[-1])+dwP,dwP))
	ybins= np.exp(np.arange(np.log(epos.yzoom[0]),np.log(epos.yzoom[-1])+dwR,dwR))

	''' Plot model occurrence or observed counts'''
	if Observation:
		# plot model planets * completeness
		weights= eta *epos.occurrence['model']['completeness'] / pfm['ns']
	else:
		weights=np.full(pfm['np'], eta/pfm['ns'])
	
	if 'draw prob' in pfm and not Tag:
		prob= pfm['draw prob'][pfm['ID']] 
		weights*= prob*pfm['ns'] # system weights sum up to 1
		#nonzero= np.where(prob>0, 1., 0.)
		#weights*= nonzero*(pfm['np']/nonzero.sum())
	
	# histograms
	if Tag:
		for key, f in subset.iteritems():
			toplot= f(pfm['tag'])
			#weights= eta*epos.occurrence['model']['completeness'] \
			#		*np.where(toplot,1.,0.)/f(pfm['system tag']).sum()
			weights= np.where(toplot,eta,0.)/f(pfm['system tag']).sum()
			axP.hist(pfm['P'], bins=xbins, weights=weights, histtype='step', label=key)
			axR.hist(pfm['R'], bins=ybins, orientation='horizontal',
				weights=weights, histtype='step')
	else:
		axP.hist(pfm['P'], bins=xbins, weights=weights, color=color)
		axR.hist(pfm['R'], bins=ybins, orientation='horizontal', weights=weights, 
					color=color)
				 
	''' Overplot observations? '''
	if Population:
		assert hasattr(epos, 'func')
		fname='.pop'
		
		ax.set_title(epos.name)
		
		pps, pdf, pdf_X, pdf_Y= periodradius(epos, Init=True)
		_, _, pdf_X, _= periodradius(epos, Init=True, ybin=ybins)
		_, _, _, pdf_Y= periodradius(epos, Init=True, xbin=xbins)
		pps, _ , _, _ = periodradius(epos, Init=True, xbin=xbins, ybin=ybins)
		#pdf/= np.max(pdf)
		#pdflog= np.log10(pdf) # in %
		levels= np.linspace(0,np.max(pdf))
		lines= np.array([0.1, 0.5]) * np.max(pdf)
		
		ax.contourf(epos.X_in, epos.Y_in, pdf, cmap='Purples', levels=levels)
		#ax.contour(epos.X_in, epos.Y_in, pdf, levels=lines)
		
		# Side panels
		#print 'pps model= {}'.format(eta)
		scale=dwR
		axP.plot(epos.MC_xvar, pdf_X*scale, marker='',ls='-',color='purple')
		axR.plot(pdf_Y*scale, epos.in_yvar, marker='',ls='-',color='purple')
	elif Observation:
		fname='.obs'
		ax.set_title(epos.name+': Counts')
		
		ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', ms=5.0, color='0.5')
		
		weights= np.full(epos.obs_xvar.size, 1./epos.nstars)
		axP.hist(epos.obs_xvar,bins=xbins,weights= weights, histtype='step', color='0.5')
		axR.hist(epos.obs_yvar,bins=ybins,weights= weights, 
			orientation='horizontal', histtype='step', color='0.5')
					
	elif Occurrence:
		fname='.occ'
		ax.set_title(epos.name+': Occurrence')

		ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', ms=5.0, color='0.5')
		
		cut= epos.obs_yvar > 0.45
		
		weights= 1. / (epos.occurrence['planet']['completeness'][cut]*epos.nstars)
		axP.hist(epos.obs_xvar[cut],bins=xbins,weights=weights,histtype='step',color='k')
		axR.hist(epos.obs_yvar[cut],bins=ybins,weights= weights, 
			orientation='horizontal', histtype='step', color='k')
			
	elif Tag:
		fname='.tag'
		ax.set_title(epos.name+': Tag')
		
		axP.legend(frameon=False, fontsize='small')
# 		for k, tag in enumerate(subset):
# 			axP.text(0.98,0.95-0.05*k,tag,ha='right',va='top',color='C1', 
# 				transform=axP.transAxes)
		
	else:
		fname=''	

	''' plot main panel'''
	#helpers.set_axes(ax, epos, Trim=True)

	helpers.set_axes(ax, epos, Trim=True)
	if Tag:
		for key, f in subset.iteritems():
			todraw= f(pfm['tag'])
			ax.plot(pfm['P'][todraw],pfm['R'][todraw], **fmt_symbol)	
	elif 'draw prob' in pfm:
		#fmt_symbol['alpha']= 0.6*pfm['draw prob'][pfm['ID']] # alpha can't be array
		todraw= pfm['draw prob'][pfm['ID']]>0
		ax.plot(pfm['P'][todraw],pfm['R'][todraw], color=color, **fmt_symbol)
	else:
		ax.plot(pfm['P'],pfm['R'], color=color, **fmt_symbol)

	
	''' Period side panel '''
	#axP.yaxis.tick_right()
	#axP.yaxis.set_ticks_position('both')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')
	helpers.set_axis_distance(axP, epos, Trim=True)
	
	''' Mass side panel'''
	helpers.set_axis_size(axR, epos, Trim=True) #, In= epos.MassRadius)
	
	helpers.save(plt, '{}model/input.radius{}'.format(epos.plotdir, fname))

def inclination(epos, color='C1', imin=1e-2):
	pfm=epos.pfm
	
	gs = gridspec.GridSpec(1, 2,
                       width_ratios=[20,6],
                       )
	f= plt.figure()
	f.subplots_adjust(wspace=0)
	
	ax = plt.subplot(gs[0, 0])	
	axh = plt.subplot(gs[0, 1])
	axh.tick_params(direction='in', which='both', left=False, right=True, labelleft=False)
	axh.yaxis.set_label_position('right')
	axh.axis('off')

	''' Inc-sma'''	
	ax.set_title('Input Inclination {}'.format(epos.name))

	ax.set_xlabel('Semi-Major Axis [au]')
	#ax.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
	ax.set_ylabel('Inclination [degree]')

	#ax.set_xlim(epos.mod_xlim)
	ax.set_ylim(imin,90)

	ax.set_xscale('log')
	ax.set_yscale('log')

	ax.axhline(np.median(pfm['inc']), ls='--', color=color)
	ax.plot(pfm['sma'], pfm['inc'], color=color, **fmt_symbol)

	''' Histogram'''
	axh.set_yscale('log')
	axh.set_ylim(imin,90)

	inc= np.logspace(np.log10(imin),2)
	axh.hist(pfm['inc'], bins=inc, orientation='horizontal', color=color) 
	
	#Model best-fit
	xmax= axh.get_xlim()[-1]
	for scale, ls in zip([1,2,2.7],[':','--',':']):
		pdf= scipy.stats.rayleigh(scale=scale).pdf(inc)
		#pdf/= np.log(inc) #log scale?
		pdf*= xmax/max(pdf)
		axh.plot(pdf, inc, ls=ls, color='purple')
		
		ax.axhline(scale, color='purple', ls=ls)	
	
# 		if ls=='--':
# 			cdf= np.cumsum(pdf*inc)
# 			h1, h2= np.interp([0.05*cdf[-1],0.95*cdf[-1]], cdf, inc)
# 			ax.axhspan(h1,h2, facecolor='purple', alpha=0.5,lw=0.1)
	

	helpers.save(plt, epos.plotdir+'model/inc-sma')
	
def periodratio(epos, color='C1'):
	pfm=epos.pfm
	
	f, (ax, axR, axP)= helpers.make_panels_right(plt)

	''' Inc-sma'''	
	ax.set_title('Input Multi-planets {}'.format(epos.name))

	axP.set_xlabel('Semi-Major Axis [au]')
	ax.set_ylabel('Period ratio')

	#ax.set_xlim(epos.mod_xlim)
	ax.set_ylim(0.9,10)

	ax.set_xscale('log')
	ax.set_yscale('log')

	#ax.axhline(np.median(pfm['inc']), ls='--')
	single= pfm['dP'] == np.nan
	inner= pfm['dP'] == 1
	nth= pfm['dP']> 1
	
	# print pfm['dP'][single]
	# print pfm['dP'][nth]
	# print pfm['dP'][inner]
	# print pfm['dP'][nth].size, pfm['np'] # ok
	
	ax.plot(pfm['sma'][single], pfm['dP'][single], color='0.7', **fmt_symbol)
	ax.plot(pfm['sma'][nth], pfm['dP'][nth], color=color, **fmt_symbol)
	ax.plot(pfm['sma'][inner], pfm['dP'][inner], color='C1', **fmt_symbol)

	''' Histogram Period Ratio'''
	axR.set_yscale('log')
	axR.set_ylim(0.9,10)

	dP= np.logspace(0,1,15)
	axR.hist(pfm['dP'][nth], bins=dP, orientation='horizontal', color=color)
	ax.axhline(np.median(pfm['dP'][nth]), ls='--', color=color)
	
	#Model best-fit
	dP= np.logspace(0,1)
	xmax= axR.get_xlim()[-1]
	with np.errstate(divide='ignore'): 
		Dgrid= np.log10(2.*(dP**(2./3.)-1.)/(dP**(2./3.)+1.))
	Dgrid[0]= -2

	for scale, ls in zip([-0.30,-0.37,-0.41],[':','--',':']):
		pdf= scipy.stats.norm(scale,0.19).pdf(Dgrid)
		pdf*= xmax/max(pdf)
		axR.plot(pdf, dP, ls=ls, color='purple')
		
		pscale= np.interp(scale,Dgrid,dP)
		ax.axhline(pscale, color='purple', ls=ls)	

	''' Histogram Inner Planet'''
	axP.set_xscale('log')
	axP.set_xlim(ax.get_xlim())
	sma= np.geomspace(*ax.get_xlim(), num=15)
	axP.hist(pfm['sma'][inner], bins=sma, color='C1', label='Inner Planet')
	
	ymax= axP.get_ylim()[-1]
	P= np.geomspace(0.5,730)
	smaP= (P/365.25)**(1./1.5)
	pdf= brokenpowerlaw1D(P, 10,1.5,-0.8)
	pdf*= ymax/max(pdf)
	
	#axP.legend(frameon=False)
	axP.text(0.98,0.95,'Inner Planet',ha='right',va='top',color='C1', 
		transform=axP.transAxes)
	
	axP.plot(smaP, pdf, ls='-', color='purple')

	
	#helpers.set_axes(ax, epos, Trim=True)
	#helpers.set_axis_distance(axP, epos, Trim=True)
	#helpers.set_axis_size(axR, epos, Trim=True) #, In= epos.MassRadius)

	helpers.save(plt, epos.plotdir+'model/Pratio-sma')

def periodratio_size(epos, color='C1'):
	pfm=epos.pfm
	
	f, (ax, axR, axP)= helpers.make_panels_right(plt)

	''' Inc-sma'''	
	ax.set_title('Input Multi-planets {}'.format(epos.name))

	axP.set_xlabel('Period ratio')
	#ax.set_ylabel(r'Size [$R_\bigoplus$]')

	ax.set_xlim(0.9,10)
	#ax.set_ylim(0.3,20)

	ax.set_xscale('log')
	#ax.set_yscale('log')

	helpers.set_axis_size(ax, epos, Trim=True)

	''' grids '''
	dP= np.logspace(0,1)
	dP_bins= np.logspace(0,1,15)

	# exoplanet data + hist
	ax.plot(epos.multi['Pratio'], epos.multi['Rpair'], ls='', marker='.', ms=5.0, color='0.5')

	#ax.axhline(np.median(pfm['inc']), ls='--')
	single= pfm['dP'] == np.nan
	inner= pfm['dP'] == 1
	nth= pfm['dP']> 1
	
	# print pfm['dP'][single]
	# print pfm['dP'][nth]
	# print pfm['dP'][inner]
	# print pfm['dP'][nth].size, pfm['np'] # ok
	
	ax.plot(pfm['dP'][single], pfm['R'][single], color='0.7', **fmt_symbol)
	ax.plot(pfm['dP'][nth], pfm['R'][nth], color=color, **fmt_symbol)
	ax.plot(pfm['dP'][inner], pfm['R'][inner], color='C1', **fmt_symbol)

	''' Histogram Period Ratio'''
	axP.set_xscale('log')
	axP.set_xlim(0.9,10)

	axP.hist(pfm['dP'][nth], bins=dP_bins, color=color)
	ax.axvline(np.median(pfm['dP'][nth]), ls='--', color=color)
	
	''' Model best-fit '''
	xmax= axP.get_ylim()[-1]
	with np.errstate(divide='ignore'): 
		Dgrid= np.log10(2.*(dP**(2./3.)-1.)/(dP**(2./3.)+1.))
	Dgrid[0]= -2

	for scale, ls in zip([-0.30,-0.37,-0.41],[':','--',':']):
		pdf= scipy.stats.norm(scale,0.19).pdf(Dgrid)
		pdf*= xmax/max(pdf)
		axP.plot(dP, pdf, ls=ls, color='purple')
		
		pscale= np.interp(scale,Dgrid,dP)
		ax.axvline(pscale, color='purple', ls=ls)	

	''' Raw data'''
	scale= 1.* pfm['dP'][nth].size/ epos.multi['Pratio'].size
	weights= np.full(epos.multi['Pratio'].size, scale)
	axP.hist(epos.multi['Pratio'], bins=dP_bins, weights=weights, histtype='step', color='0.5', zorder=1)

	''' Histogram Planet Size'''
	axR.set_yscale('log')
	axR.set_ylim(ax.get_ylim())
	radius= np.geomspace(*ax.get_ylim(), num=15)
	axR.hist(pfm['R'][inner], bins=radius, orientation='horizontal', color='C1', label='Inner Planet')
		
	#helpers.set_axes(ax, epos, Trim=True)
	#helpers.set_axis_distance(axP, epos, Trim=True)

	''' Linear regression '''
	try:
		# data
		slope, intercept, r_value, p_value, std_err= \
			scipy.stats.linregress(np.log(epos.multi['Pratio']), np.log(epos.multi['Rpair']))
		ax.plot(dP, np.exp(intercept+slope*np.log(dP)), label='r={:.2f}'.format(r_value), 
				marker='', ls='-', color='0.5')
		#print slope, intercept, r_value, p_value

		slope, intercept, r_value, p_value, std_err= \
			scipy.stats.linregress(np.log(pfm['dP'][nth]), np.log(pfm['R'][nth]))
		ax.plot(dP, np.exp(intercept+slope*np.log(dP)), label='r={:.2f}'.format(r_value), 
				marker='', ls='-', color=color)

		ax.legend(loc='lower right')
	except Exception as e: print(e)

	helpers.save(plt, epos.plotdir+'model/Pratio-size')


def multiplicity(epos, color='C1', Planets=False, Kepler=False):
	# plot multiplicity
	f, ax = plt.subplots()
	ax.set_title('Input Multi-planets {}'.format(epos.name))
	ax.set_xlabel('planets per system')
	ax.set_ylabel('number of planets' if Planets else 'number of systems')

	if Kepler:
		ax.set_xlim(0.5, 7.5)	
	#else:
		#ax.set_xlim(0.5, 10.5)

	''' Model planets '''
	_ , counts= np.unique(epos.pfm['ID'],return_counts=True)
	#bins= np.arange(10)
	#ax.hist(counts, bins=bins, align='left', color=color)

	if Kepler:
		bins= np.arange(1,9)
		bins[-1]=1000
	else:
		bins= np.arange(20)
	hist, bin_edges = np.histogram(counts, bins=bins)
	ax.bar(bins[:-1], hist, width=1, color=color)

	''' Kepler '''
# 	intrinsic= np.zeros_like(hist)
# 	intrinsic[0]= 0.5*epos.pfm['ns']
# 	intrinsic[-1]= 0.5*epos.pfm['ns']
# 	#intrinsic=np.arange(len(hist))
# 	ax.bar(bins[:-1], intrinsic, width=1, color='', ec='k')
	
	#ax.plot(, drawstyle='steps-mid', 
	#	ls='--', marker='', color='gray', label='Kepler all')

	ax.legend(loc='upper right', shadow=False, prop={'size':14}, numpoints=1)

	helpers.save(plt, epos.plotdir+'model/multi')

def period(epos, Population=False, Occurrence=False, Observation=False, Tag=False, color='C1'):
	''' Model occurrence as function of orbital period'''
	f, ax = plt.subplots()
	helpers.set_axis_distance(ax, epos, Trim=True)
	ax.set_xlim(0.5,200)
	
	ax.set_ylabel(r'Planet Occurrence 1-4 $R_\bigoplus$')
	ax.set_yscale('log')
	
	pfm=epos.pfm
	eta= epos.modelpars.get('eta',Init=True)
	
	if not 'R' in pfm:
		pfm['R'], _= epos.MR(pfm['M'])

	if Tag:
		# function that return a simulation subset based on the tag
		subset={'Fe/H<=0': lambda tag: tag<=0,
			'Fe/H>0': lambda tag: tag>0}

	''' Bins '''
	#xbins= np.geomspace(1,1000,10)/(10.**(1./3.))
	xbins= np.geomspace(0.5,200,15)


	''' Plot model occurrence or observed counts'''
	weights=np.full(pfm['np'], eta/pfm['ns'])
	
	if 'draw prob' in pfm and not Tag:
		prob= pfm['draw prob'][pfm['ID']] 
		weights*= prob*pfm['ns'] # system weights sum up to 1
	
	# histograms
	if Tag:
		for key, f in subset.iteritems():
			toplot= f(pfm['tag']) #& (pfm['R']>1.)
			weights= np.where(toplot,eta,0.)/f(pfm['system tag']).sum()			

			ax.hist(pfm['P'], bins=xbins, weights=weights, histtype='step', label=key,
				color='#88CCEE' if key=='Fe/H<=0' else '#332288')
	else:
		ax.hist(pfm['P'], bins=xbins, weights=weights, color=color)
				 
			
	if Tag:
		fname='.tag'
		ax.set_title(epos.name+': Disk Fe/H')
		ax.legend(fontsize='small', loc='lower right')		
	else:
		fname=''	

	
	helpers.save(plt, '{}model/input.period{}'.format(epos.plotdir, fname))
