import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.patches as patches
import helpers
from EPOS import regression

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def periodradius(epos, SNR=True, Parametric=False):

	if SNR:
		sim=epos.synthetic_survey
		suffix=''
		fsuffix=''
	else:
		sim=epos.transit
		suffix=' (no SNR)'
		fsuffix='_noSNR'
	
	# plot R(P)
	f, ax = plt.subplots()
	ax.set_title('detectable planet sample'+suffix)
	helpers.set_axes(ax, epos, Trim=True)
	if Parametric or len(epos.groups)==1:
		ax.plot(sim['P'], sim['Y'], ls='', marker='.', mew=0, ms=5.0, color='k')
	else:
		for k, sg in enumerate(epos.groups):
			subset= sim['i sg']==k
			ax.plot(sim['P'][subset], sim['Y'][subset], ls='', marker='.', mew=0, ms=5.0, color=clrs[k % 4], label=sg['name'])

	
	#ax.legend(loc='lower left', shadow=False, prop={'size':14}, numpoints=1)
	helpers.save(plt, '{}output/periodradius{}'.format(epos.plotdir,fsuffix))

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
	if epos.populationtype is 'model':
		for k, sg in enumerate(epos.groups):
			subset= sim['i sg']==k
			P= sim['P'][subset]
			R= sim['Y'][subset]
			ax.plot(np.log10(P),np.log10(R), zs=0,zdir='z',
				ls='',marker='.',mew=0,ms=5.0,color=clrs[k % 4])
				
			xgrid= np.logspace(*np.log10(epos.xtrim))
			pdf= regression.sliding_window_log(P, None, xgrid) #, width=2. )
			ax.plot(np.log10(xgrid), pdf, zs=yplane,zdir='y', 
				ls='-', marker='', color=clrs[k % 4],
				label='{} x{:.3f}'.format(sg['name'], sg['weight']))

			ygrid= np.logspace(*np.log10(epos.ytrim))
			pdf= regression.sliding_window_log(R, None, ygrid) #, width=2. )
			ax.plot(np.log10(ygrid), pdf, zs=xplane, zdir='x', 
				ls='-', marker='', color=clrs[k % 4])
	else:
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

def cdf(epos):
	
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	f.set_size_inches(9, 7) # default 7, 5

	''' 
	top left: synthetic obsservation
	'''
	sim= epos.synthetic_survey	
	ax1.set_title('Synthetic ({})'.format(sim['P zoom'].size))
	helpers.set_axes(ax1, epos, Trim=True)
	ax1.plot(sim['P'], sim['Y'], ls='', marker='.', mew=0, ms=5.0, color='r', alpha=0.5)
		
	if epos.Zoom:
		ax1.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
			epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
		
	''' 
	Top right: observed population 
	'''
	ax2.set_title('Observed ({})'.format(epos.obs_zoom['x'].size))
	helpers.set_axes(ax2, epos, Trim=True)

	ax2.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='b')		
	if epos.Zoom:
		ax2.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
			epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )

	''' 
	cdf orbital period 
	'''
	ax3.set_title('Period, p={:.3g}'.format(epos.gof['xvar']))
	ax3.set_xlabel('Orbital Period [days]')
	ax3.set_ylabel('CDF')
	ax3.set_xscale('log')
	ax3.set_ylim([-0.05,1.05])

	ax3.set_xlim(*epos.xzoom)

	#model histogram x
	P= sim['P zoom']
	ax3.plot(np.sort(P), np.arange(P.size, dtype=float)/P.size, ls='-', marker='', color='r')
	
	#obs histogram x
	P= epos.obs_zoom['x']
	ax3.plot(np.sort(P), np.arange(P.size, dtype=float)/P.size, ls='-', marker='', color='b')		

	''' 
	CDF planet radius
	'''
	ax4.set_title('{}, p={:.3g}'.format('M sin i' if epos.RV else 'Radius', epos.gof['yvar']))
	if epos.RV:	ax4.set_xlabel(r'Planet M sin i [M$_\bigoplus$]')
	else:		ax4.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
	ax4.set_ylabel('CDF')
	ax4.set_xscale('log')
	ax4.set_ylim([-0.05,1.05])
	ax4.set_xlim(*epos.yzoom)
	
	ax4.set_xticks(epos.yticks)
	ax4.get_xaxis().set_major_formatter(tck.ScalarFormatter())

	#model histogram x
	R= sim['Y zoom']
	ax4.plot(np.sort(R), np.arange(R.size, dtype=float)/R.size, ls='-', marker='', color='r')

	#obs histogram x
	R= epos.obs_zoom['y']
	ax4.plot(np.sort(R), np.arange(R.size, dtype=float)/R.size, ls='-', marker='', color='b')		
		
	f.tight_layout()
	helpers.save(plt, epos.plotdir+'output/cdf.diag')