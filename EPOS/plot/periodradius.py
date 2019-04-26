import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.patches as patches
import helpers
from EPOS import regression

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

# backwards compatible colors 
import matplotlib
if matplotlib.__version__[0] != 2: 
	helpers.default_pyplot2_colors(matplotlib.colors)

def periodradius(epos, SNR=True, Model=False, color='C1'):

	f, (ax, axR, axP)= helpers.make_panels(plt, Fancy=True)
	
	sim=epos.synthetic_survey
	title='Simulated Detections: {}'.format(epos.name)
	fsuffix='detect'
	xbins= epos.MC_xvar
	ybins= epos.MC_yvar

	if SNR:
		transit=epos.transit
		title= 'Simulated Transiting Planets'
		fsuffix='transit'
	elif Model:
		pfm= epos.pfm
		fsuffix='population'

		dwR=0.2 # bin width in ln space
		dwP=0.3
		
		xbins= np.exp(np.arange(np.log(epos.mod_xlim[0]),np.log(epos.mod_xlim[-1])+dwP,dwP))
		ybins= np.exp(np.arange(np.log(epos.ytrim[0]),np.log(epos.ytrim[-1])+dwR,dwR))	

	
	''' plot R(P), main panel'''
	f.suptitle(title)
	helpers.set_axes(ax, epos, Trim=True)
	if SNR: 
		ax.plot(transit['P'], transit['Y'], ls='', marker='.', color='C6')
	elif Model:
		ax.plot(pfm['P'],pfm['R'], ls='', marker='.', ms=5.0, color='0.5')
	ax.plot(sim['P'], sim['Y'], ls='', marker='.', color=color)

	if Model:
		xlim= ax.get_xlim()
		ax.set_xlim(xlim[0], epos.mod_xlim[-1])
		#ax.set_ylim(epos.mod_ylim)

	''' Period side panel '''
	#helpers.set_axis_distance(axP, epos, Trim=True)
	axP.set_xlabel(ax.get_xlabel())
	#axP.set_yscale('log')
	#axP.set_ylim([2e-3,5])	
	#axP.set_yticks([0.01,0.1,1])
	#axP.set_yticklabels(['1%','10%','100%'])
	axP.yaxis.tick_right()
	axP.yaxis.set_ticks_position('both')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')
	
	if SNR: 
		axP.hist(transit['P'], bins=xbins, color='C6', density=True)
	elif Model:
		axP.hist(pfm['P'], bins=xbins, color='0.5', density=True, log=True)
		axP.hist(pfm['P'], bins=xbins, color='k', density=True, log=True, 
			histtype='step', zorder=10, ls='-')
	axP.hist(sim['P'], bins=xbins, density=Model, log=Model, color=color)


	''' Radius side panel'''
	#helpers.set_axis_size(axR, epos, Trim=True, In= epos.MassRadius)
	axR.set_ylabel(ax.get_ylabel())

	#axR.set_xscale('log')
	#axR.set_xlim([2e-3,5])
	#axR.set_xticks([1,10,100,1000])
	#axR.set_xticklabels(['1','10','100','1000'], rotation=70)
	for tick in axR.get_xticklabels():
		tick.set_rotation(70)
	#axR.tick_params(axis='x', which='minor',top='off',bottom='off')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')

	if SNR: 
		axR.hist(transit['Y'],orientation='horizontal', bins=ybins, color='C6', 
			density=True)
	elif Model:
		axR.hist(pfm['R'],orientation='horizontal', bins=ybins, color='0.5', 
			density=True, log=True)
		axR.hist(pfm['R'],orientation='horizontal', bins=ybins, color='k', 
			density=True, log=True, histtype='step', zorder=10, ls='-')
	axR.hist(sim['Y'],orientation='horizontal', bins=ybins, 
		density=Model, log=Model, color=color)
	
	# labels
	if SNR:
		xpos= epos.MC_xvar[-1]
		ypos= axP.get_ylim()[-1]/1.05
		axP.text(xpos, ypos, 'Detected ', color=color, ha='right', va='top')
		axP.text(xpos, ypos/1.3, 'Undetected ', color='C6', ha='right', va='top')
	elif Model:
		xpos= ax.get_xlim()[-1]
		ypos= ax.get_ylim()[-1]/1.05
		ax.text(xpos, ypos, 'Detected ', color=color, ha='right', va='top')
		ax.text(xpos, ypos/1.3, 'Undetected ', color='0.5', ha='right', va='top')

	
	#ax.legend(loc='lower left', shadow=False, prop={'size':14}, numpoints=1)
	helpers.save(plt, '{}output/periodradius.{}'.format(epos.plotdir,fsuffix))

def panels(epos, MCMC=False, color='C1'):

	f, (ax, axR, axP)= helpers.make_panels(plt)
	
	sim=epos.synthetic_survey
	
	clr_bf= 'g'
		
	''' plot R(P), main panel'''
	ax.set_title('Simulated Detections')
	helpers.set_axes(ax, epos, Trim=True)
	if epos.MonteCarlo:
		ax.plot(sim['P'], sim['Y'], ls='', marker='.', color=clr_bf if MCMC else 'C0')
	else:
		levels= np.linspace(0,np.max(sim['pdf']))		
		ax.contourf(epos.MC_xvar, epos.MC_yvar, sim['pdf'].T, cmap='Blues', levels=levels)

	''' Period side panel '''
	helpers.set_axis_distance(axP, epos, Trim=True)
	#axP.set_yscale('log')
	#axP.set_ylim([2e-3,5])	
	#axP.set_yticks([0.01,0.1,1])
	#axP.set_yticklabels(['1%','10%','100%'])
	axP.yaxis.tick_right()
	axP.yaxis.set_ticks_position('both')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')
			
	''' Radius side panel'''
	helpers.set_axis_size(axR, epos, Trim=True) #, In= epos.MassRadius)

	#axR.set_xscale('log')
	#axR.set_xlim([2e-3,5])
	#axR.set_xticks([1,10,100,1000])
	#axR.set_xticklabels(['1','10','100','1000'], rotation=70)
	for tick in axR.get_xticklabels():
		tick.set_rotation(70)
	#axR.tick_params(axis='x', which='minor',top='off',bottom='off')
	#axP.tick_params(axis='y', which='minor',left='off',right='off')

	''' Histograms / posterior samples '''
	try:
		xbins= np.geomspace(*epos.xzoom, num=20)
		ybins= np.geomspace(*epos.yzoom, num=10)
	except:
		xbins= np.logspace(*np.log10(epos.xzoom), num=20)
		ybins= np.logspace(*np.log10(epos.yzoom), num=10)
	xscale= np.log(xbins[1]/xbins[0])
	yscale= np.log(ybins[1]/ybins[0])

	if MCMC:
		histkeys= {'color':'b', 'alpha':0.1}
		for ss in epos.ss_sample:
			if epos.MonteCarlo:
				axP.hist(ss['P zoom'], bins=xbins, histtype='step', **histkeys)
				axR.hist(ss['Y zoom'], bins=ybins, orientation='horizontal', \
					histtype='step', **histkeys)
			else:
				axP.plot(ss['P zoom'], ss['P zoom pdf']*xscale, 
					marker='', ls='-', **histkeys)
				axR.plot(ss['Y zoom pdf']*yscale, ss['Y zoom'], 
					marker='', ls='-', **histkeys)
		histdict= {'histtype':'step', 'color':clr_bf}
	else:
		histdict={}
	
	if epos.MonteCarlo:
		axP.hist(sim['P zoom'], bins=xbins, **histdict)
		axR.hist(sim['Y zoom'], bins=ybins, orientation='horizontal', **histdict)
	else:
		axP.plot(sim['P zoom'], sim['P zoom pdf']*xscale, marker='', ls='-')
		axR.plot(sim['Y zoom pdf']*yscale, sim['Y zoom'], marker='', ls='-')

	''' Observations'''
	axP.hist(epos.obs_zoom['x'], bins=xbins,histtype='step', color='C1')
	axR.hist(epos.obs_zoom['y'], bins=ybins, orientation='horizontal',histtype='step', color='C1')
	
	''' Box/ lines'''
	if epos.Zoom:
		for zoom in epos.xzoom: axP.axvline(zoom, ls='--', color='k')
		for zoom in epos.yzoom: axR.axhline(zoom, ls='--', color='k')
		ax.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
			epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
	
	#ax.legend(loc='lower left', shadow=False, prop={'size':14}, numpoints=1)
	
	fdir= 'mcmc' if MCMC else 'output'
	helpers.save(plt, '{}{}/pdf.zoom'.format(epos.plotdir, fdir))

def pdf_3d(epos):
	#plot in 3D
	# doesn't support log scale on axes, need a workaround
	
	fig = plt.figure()
	ax = fig.gca(projection='3d')

	ax.set_title('Synthetic Model Populations')
	ax.set_xlabel('Orbital Period [days]')
	ax.set_zlabel('Planets/bin')
	
	if epos.RV:	ax.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
	else:		ax.set_ylabel(r'Planet Radius [R$_\bigoplus$]')
	
	ax.set_xlim(np.log10(epos.xtrim))
	ax.set_ylim(np.log10(epos.ytrim))
	xplane, yplane= np.log10(epos.xtrim[0]), np.log10(epos.ytrim[-1])
	ax.set_zlim(0,2000)

	#ax.set_xscale('log')
	#ax.set_yscale('log')

	sim=epos.synthetic_survey
	
	# PDF, individual contributions
	# 	if epos.populationtype is 'model':
	# 		for k, sg in enumerate(epos.groups):
	# 			subset= sim['i sg']==k
	# 			P= sim['P'][subset]
	# 			R= sim['Y'][subset]
	# 			ax.plot(np.log10(P),np.log10(R), zs=0,zdir='z',
	# 				ls='',marker='.',mew=0,ms=5.0,color=clrs[k % 4])
	# 				
	# 			xgrid= np.logspace(*np.log10(epos.xtrim))
	# 			pdf= regression.sliding_window_log(P, None, xgrid) #, width=2. )
	# 			ax.plot(np.log10(xgrid), pdf, zs=yplane,zdir='y', 
	# 				ls='-', marker='', color=clrs[k % 4],
	# 				label='{} x{:.3f}'.format(sg['name'], sg['weight']))
	# 
	# 			ygrid= np.logspace(*np.log10(epos.ytrim))
	# 			pdf= regression.sliding_window_log(R, None, ygrid) #, width=2. )
	# 			ax.plot(np.log10(ygrid), pdf, zs=xplane, zdir='x', 
	# 				ls='-', marker='', color=clrs[k % 4])
	# 	else:

	# top left panel (P,R)
	ax.plot(np.log10(sim['P']), np.log10(sim['Y']),zs=0,zdir='z', 
			ls='', marker='.', mew=0, ms=5.0, color='k')
		
	# PDF, all combined, 2 panels
	xgrid= np.logspace(*np.log10(epos.xtrim))
	pdf= regression.sliding_window_log(sim['P'], None, xgrid) #, width=2. )
	weight= np.sum([sg['weight'] for sg in epos.groups]) if epos.populationtype is 'model' else epos.fitpars.get('pps',Init=True)
	ax.plot(np.log10(xgrid), pdf, zs=yplane,zdir='y', 
		ls='-', marker='', color='k',label='combined x {:.3f}'.format(weight))

	ygrid= np.logspace(*np.log10(epos.ytrim))
	pdf= regression.sliding_window_log(sim['Y'], None, ygrid) #, width=2. )
	ax.plot(np.log10(ygrid), pdf, zs=xplane,zdir='x', ls='-', marker='', color='k')

	# observations
	if epos.Observation and epos.DetectionEfficiency:
		ax.plot(np.log10(epos.obs_xvar), np.log10(epos.obs_yvar), zs=-10,zdir='z',
			ls='', marker='.', mew=0, ms=5.0, color='0.7',zorder=0)
		
		pdf= regression.sliding_window_log(epos.obs_xvar, None, xgrid)
		ax.plot(np.log10(xgrid), pdf, zs=yplane,zdir='y', ls=':', marker='', color='k', label='Kepler')

		pdf= regression.sliding_window_log(epos.obs_yvar, None, ygrid) 
		ax.plot(np.log10(ygrid), pdf, zs=xplane,zdir='x', ls=':', marker='', color='k')
		
	# 	xmax= ax.get_ylim()[-1]
	# 	ymax= ax.get_zlim()[-1]
	# 	ax.errorbar(xmax/1.5, ymax*0.9, yerr=(xmax/1.5*(1.-np.sqrt(1./2.)) ), 
	# 				zs=0,zdir='y', color='k')
	# 	
	# 	xmax= ax.get_xlim()[-1]
	# 	ymax= ax.get_zlim()[-1]
	# 	ax.errorbar(xmax/1.5, ymax*0.9, xerr=(xmax/1.5*(1.-np.sqrt(1./2.)) ), 
	# 				zs=0,zdir='x', color='k')
		
	#ax2.tick_params(labelbottom='on') # does not re-enable axis
	
	''' Legend instead of 4th plot'''
	#ax.legend()
	ax.legend(shadow=False, prop={'size':14},bbox_to_anchor=(0.7,0.0))
	
	ax.view_init(elev=20., azim=-35)
	
	helpers.save(plt, epos.plotdir+'output/pdf.3D.diag')

def pdf(epos):
	#assert epos.populationtype is 'model'
	
	f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
	f.subplots_adjust(hspace=0, wspace=0)
	# 1 2
	# 3 4
	
	ax1.set_title('Synthetic Model Populations',loc='left')

	ax3.set_xlabel('Orbital Period [days]')
	if epos.RV:	ax1.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
	else:		ax1.set_ylabel(r'Planet Radius [R$_\bigoplus$]')
	
	ax2.set_xlabel('Planets/bin')
	ax3.set_ylabel('Planets/bin')
	
	ax1.set_xlim(epos.xtrim)
	ax1.set_ylim(epos.ytrim)

	ax1.set_xscale('log')
	ax1.set_yscale('log')

	sim=epos.synthetic_survey
	
	# PDF, individual contributions
	if epos.populationtype is 'model':
		for k, sg in enumerate(epos.groups):
			subset= sim['i sg']==k
			P= sim['P'][subset]
			R= sim['Y'][subset]
			ax1.plot(P,R, ls='', marker='.', mew=0, ms=5.0, color=clrs[k % 4])
				
			xgrid= np.logspace(*np.log10(epos.xtrim))
			pdf= regression.sliding_window_log(P, None, xgrid) #, width=2. )
			ax3.plot(xgrid, pdf, ls='-', marker='', color=clrs[k % 4],label='{} x{:.3f}'.format(sg['name'], sg['weight']))

			ygrid= np.logspace(*np.log10(epos.ytrim))
			pdf= regression.sliding_window_log(R, None, ygrid) #, width=2. )
			ax2.plot(pdf, ygrid, ls='-', marker='', color=clrs[k % 4])
	else:
		# top left panel (P,R)
		ax1.plot(sim['P'], sim['Y'], ls='', marker='.', mew=0, ms=5.0, color='k')
		
	# PDF, all combined, 2 panels
	xgrid= np.logspace(*np.log10(epos.xtrim))
	pdf= regression.sliding_window_log(sim['P'], None, xgrid) #, width=2. )
	weight= np.sum([sg['weight'] for sg in epos.groups]) if epos.populationtype is 'model' else epos.fitpars.get('pps',Init=True)
	ax3.plot(xgrid, pdf, ls='-', marker='', color='k',label='combined x {:.3f}'.format(weight))

	ygrid= np.logspace(*np.log10(epos.ytrim))
	pdf= regression.sliding_window_log(sim['Y'], None, ygrid) #, width=2. )
	ax2.plot(pdf, ygrid, ls='-', marker='', color='k')

	# observations
	if epos.Observation and epos.DetectionEfficiency:
		pdf= regression.sliding_window_log(epos.obs_xvar, None, xgrid)
		ax3.plot(xgrid, pdf, ls=':', marker='', color='k', label='Kepler')

		pdf= regression.sliding_window_log(epos.obs_yvar, None, ygrid) 
		ax2.plot(pdf, ygrid, ls=':', marker='', color='k')
	
	xmax= ax3.get_xlim()[-1]
	ymax= ax3.get_ylim()[-1]
	ax3.errorbar(xmax/1.5, ymax*0.9, xerr=(xmax/1.5*(1.-np.sqrt(1./2.)) ), color='k')
	
	xmax= ax2.get_xlim()[-1]
	ymin= ax2.get_ylim()[0]
	ax2.errorbar(xmax*0.9, ymin*1.5, yerr=(ymin*1.5*(1.-np.sqrt(1./2.)) ), color='k')
	
	#ax2.tick_params(labelbottom='on') # does not re-enable axis
	
	''' Legend instead of 4th plot'''
	ax3.legend(loc='center left', shadow=False, prop={'size':14}, numpoints=1,bbox_to_anchor=(1,0.5))
	ax4.axis('off')
	
	helpers.save(plt, epos.plotdir+'output/pdf.diag')

def cdf(epos, color='C1'):
	
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	f.set_size_inches(9, 7) # default 7, 5

	f.suptitle(epos.name)

	''' 
	top left: synthetic obsservation
	'''
	sim= epos.synthetic_survey
	helpers.set_axes(ax1, epos, Trim=True)
	ax1.set_title('Synthetic ({})'.format(sim['nobs']))
	if epos.MonteCarlo:
		ax1.plot(sim['P'], sim['Y'], ls='', marker='.', mew=0, ms=5.0, color=color, alpha=0.5)
	else:
		levels= np.linspace(0,np.max(sim['pdf']))		
		ax1.contourf(epos.MC_xvar, epos.MC_yvar, sim['pdf'].T, cmap='Reds', levels=levels)
	
	if epos.Zoom:
		ax1.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
			epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
		
	''' 
	Top right: observed population 
	'''
	ax2.set_title('Observed ({})'.format(epos.obs_zoom['x'].size))
	helpers.set_axes(ax2, epos, Trim=True)

	ax2.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='C3')		
	if epos.Zoom:
		ax2.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
			epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )

	''' 
	cdf orbital period 
	'''
	if 'xvar' in epos.prob:
		ax3.set_title('Period, p={:.3g}'.format(epos.prob['xvar']))
	else:
		ax3.set_title('Orbital Period')
	ax3.set_xlabel('Orbital Period [days]')
	ax3.set_ylabel('CDF')
	ax3.set_xscale('log')
	ax3.set_ylim([-0.05,1.05])

	ax3.set_xlim(*epos.xzoom)

	if epos.MonteCarlo:
		#model histogram x
		P= sim['P zoom']
		ax3.plot(np.sort(P), np.arange(P.size, dtype=float)/P.size, ls='-', marker='', color=color)
	else:
		ax3.plot(sim['P zoom'], sim['P zoom cdf'], ls='-', marker='', color='r')
	
	#obs histogram x
	P= epos.obs_zoom['x']
	ax3.plot(np.sort(P), np.arange(P.size, dtype=float)/P.size, ls='-', marker='', color='C3')		

	''' 
	CDF planet radius
	'''
	if 'yvar' in epos.prob:
		ax4.set_title('{}, p={:.3g}'.format('M sin i' if epos.RV else 'Radius', epos.prob['yvar']))
	else:
		ax4.set_title('{}'.format('M sin i' if epos.RV else 'Radius'))
	if epos.RV:	ax4.set_xlabel(r'Planet M sin i [M$_\bigoplus$]')
	else:		ax4.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
	ax4.set_ylabel('CDF')
	ax4.set_xscale('log')
	ax4.set_ylim([-0.05,1.05])
	ax4.set_xlim(*epos.yzoom)
	
	ax4.set_xticks(epos.yticks)
	ax4.get_xaxis().set_major_formatter(tck.ScalarFormatter())

	if epos.MonteCarlo:
		#model histogram x
		R= sim['Y zoom']
		ax4.plot(np.sort(R), np.arange(R.size, dtype=float)/R.size, ls='-', marker='', color=color)
	else:
		ax4.plot(sim['Y zoom'], sim['Y zoom cdf'], ls='-', marker='', color=color)

	#obs histogram x
	R= epos.obs_zoom['y']
	ax4.plot(np.sort(R), np.arange(R.size, dtype=float)/R.size, ls='-', marker='', color='C3')		
		
	f.tight_layout()
	helpers.save(plt, epos.plotdir+'output/cdf.diag')