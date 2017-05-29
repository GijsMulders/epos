#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import aux
import numpy as np
from EPOS import regression

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def readme():
	print '\nMakes all the plots'
	print 
	print 'observed:'
	print 'model:'
	print 'result:'

def input(epos):
	print '\nPlotting input...'

	observed(epos, PlotBox=False)
	detection_efficiency(epos, PlotBox=False)
	
	observed(epos, PlotBox=True)
	detection_efficiency(epos, PlotBox=True)
	
	if epos.populationtype is 'parametric':
		parametric_1D(epos)
		parametric_2D(epos)		
	elif epos.populationtype is 'model':
		population_input(epos, PlotBox=False)
		population_input_diag(epos, PlotBox=False)
		if epos.RadiusMassConversion: 
			massradius(epos, Input=True, MC=False)
			massradius(epos, Input=True, MC=False, Log=True)
		if 'all_Pratio' in epos.groups[0]:
			Pratio(epos)
		if 'all_inc' in epos.groups[0]:
			inclination(epos)
			
	else:
		assert False

def observed(epos, PlotBox=True):
	assert epos.Observation
	f, ax = plt.subplots()
	ax.set_title('Input Observed Population')
	
	_set_axes(ax, epos, Trim=True)
	
	ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='k')
	# add multis? 
	
	if PlotBox:
		fname='.box'
		if not epos.Range: epos.set_ranges()
		ax.add_patch(patches.Rectangle( (epos.xtrim[0],epos.ytrim[0]), 
			epos.xtrim[1]-epos.xtrim[0], epos.ytrim[1]-epos.ytrim[0],fill=False, zorder=1, ls='--') )
		if epos.Zoom:
			ax.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
				epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
	else: fname=''

	aux.save(plt, epos.plotdir+'input/observed_planets'+fname)
	
def detection_efficiency(epos, PlotBox=False):
	assert epos.DetectionEfficiency
	f, ax = plt.subplots()
	ax.set_title('Detection efficiency [%]')
	_set_axes(ax, epos, Trim=False, Eff=True)

	# contour with labels (not straightforward)
	# missing: vmin, vmx, a log scale, contour linestyles
	labels=['0.001','0.01', '0.1','1.0','10','100']
	n=len(labels)-1
	levels_log= np.linspace(-n,0,10*n+1)
	ticks_log= np.linspace(-n,0,n+1)
	ticks= np.logspace(-n,0, n+1)
	
	cmap = plt.cm.get_cmap('PuRd')
	cs= ax.contourf(epos.eff_xvar, epos.eff_yvar, epos.eff_2D_log.T, cmap=cmap, levels=levels_log) #vmin=-3., vmax=0.)	# vmax keywords don't work on colorbar
	cbar= f.colorbar(cs, ax=ax, shrink=0.95, ticks=ticks_log)
	cbar.ax.set_yticklabels(labels)
		
	if PlotBox:
		fname='.box'
		if not epos.Range: epos.set_ranges()
		ax.add_patch(patches.Rectangle( (epos.xtrim[0],epos.ytrim[0]), 
			epos.xtrim[1]-epos.xtrim[0], epos.ytrim[1]-epos.ytrim[0],fill=False, zorder=1, ls='--') )
		if epos.Zoom:
			ax.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
				epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
	else: 
		#ax.contour(epos.eff_xvar, epos.eff_yvar, epos.eff_2D.T,levels= ticks, colors = ['k','k','k','k'], linewidths=3.)	
		for tick, label in zip(ticks,labels):
			cs= ax.contour(epos.eff_xvar, epos.eff_yvar, epos.eff_2D.T,levels= [tick], colors = ['k'], linewidths=3.)
			plt.clabel(cs, inline=1, fmt=label+'%%')

		fname=''

	aux.save(plt, epos.plotdir+'input/detection_efficiency'+fname)
	
	pass

def massradius(epos, Input=True, MC=False, Log=False):
	assert epos.RadiusMassConversion
	# plot R(M)
	f, ax = plt.subplots()
	ax.set_title('mass-radius relation + dispersion')
	ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	ax.set_ylabel(r'Planet Radius [R$_\bigoplus$]')

	if Log:
		ax.set_xlim(0.5, 1000.) # 20.
		ax.set_ylim(0.5, 20.) # 7.	
		ax.set_xscale('log')
		ax.set_yscale('log')
	else:
		ax.set_xlim(0, 20.) # 20.
		ax.set_ylim(0, 7.) # 7.
	
	# MC data
	if MC:
		tr=epos.transit
		ax.plot(tr['M'], tr['R'], ls='', marker='.', mew=0, ms=5.0, color='k', zorder=1)
		
	if Input:	
		xM= np.logspace(*np.log10(ax.get_xlim())) if Log else np.linspace(*ax.get_xlim()) #np.max(tr['M']))
		xR, dispersion= epos.RM(xM)
		#ax.plot(xM, epos.MR['fRock'](xM), ls='-', marker='', color='r', label='Rocky constraint', zorder=0)
		#xGas= epos.MR['fGas'](xM)
		#dispersion= epos.MR['dispersion']
		ax.plot(xM, xR, ls='-', marker='', color='b', label=epos.RM_label)
		ax.plot(xM, xR-dispersion, ls='--', marker='', color='b')
		ax.plot(xM, xR+dispersion, ls='--', marker='', color='b')

		ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
	
	prefix= 'output' if MC else 'input'
	
	suffix='_log' if Log else ''	
	if (not MC) and Input: suffix+= '_nodata'
	elif MC and (not Input): suffix+= '_nomodels'
	
	aux.save(plt, '{}{}/massradius{}'.format(epos.plotdir,prefix,suffix))		

def parametric_1D(epos, PlotZoom=False):
	# these are NOT integrated over other axes (TODO) -> use epos.pdf
	
	# where does this go? parametric.py?
	assert epos.populationtype is 'parametric'
	if not epos.Range: epos.set_ranges()
	ngrid= 100.
	xgrid= np.logspace(np.log10(epos.xtrim[0]),np.log10(epos.xtrim[1]),ngrid)
	ygrid= np.logspace(np.log10(epos.ytrim[0]),np.log10(epos.ytrim[1]),ngrid)
	parax= epos.norm0* epos.xfunc(xgrid, *epos.xp0)
	paray= epos.norm0* epos.yfunc(ygrid, *epos.yp0)
	
	f, ax = plt.subplots()
	ax.set_title('Starting guess planet population')
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('Occurrence []')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.xtrim)
	ax.set_ylim([1e-3,1e2])
	ax.plot(xgrid, 100.* parax, marker='',ls=':',color='k')
	aux.save(plt, epos.plotdir+'input/parametric_initial_x')

	f, ax = plt.subplots()
	ax.set_title('Starting guess planet population')
	if epos.Radius:	ax.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
	else:			ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	ax.set_ylabel('Occurrence []')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.ytrim)
	ax.set_ylim([1e-3,1e2])
	ax.plot(ygrid, 100.* paray, marker='',ls=':',color='k')
	aux.save(plt, epos.plotdir+'input/parametric_initial_y')

def parametric_2D(epos, PlotZoom=False):
	
	# where does this go -> run.py
	assert epos.populationtype is 'parametric'
	if not epos.Range: epos.set_ranges()
	ngrid= 100. # use native MCgrid instead
	xgrid= np.logspace(np.log10(epos.xtrim[0]),np.log10(epos.xtrim[1]),ngrid)
	ygrid= np.logspace(np.log10(epos.ytrim[0]),np.log10(epos.ytrim[1]),ngrid)
	X,Y= np.meshgrid(xgrid, ygrid)
	para2D= epos.norm0* epos.xfunc(X, *epos.xp0)* epos.yfunc(Y,*epos.yp0)
	paralog= np.log10(para2D)
	
	f, ax = plt.subplots()
	ax.set_title('Starting guess planet population')
	_set_axes(ax, epos, Trim=True)

	# contour with labels (not straightforward)
	# missing: vmin, vmx, a log scale, contour linestyles
	labels=['0.001','0.01', '0.1','1.0','10','100']
	n=len(labels)-1
	levels_log= np.linspace(-n,0,10*n+1)
	ticks_log= np.linspace(-n,0,n+1)
	ticks= np.logspace(-n,0, n+1)
	
	cmap = plt.cm.get_cmap('jet')
	cs= ax.contourf(X, Y, paralog, cmap=cmap, levels=levels_log) #vmin=-3., vmax=0.)	# vmax keywords don't work on colorbar
	cbar= f.colorbar(cs, ax=ax, shrink=0.95, ticks=ticks_log)
	cbar.ax.set_yticklabels(labels)

	aux.save(plt, epos.plotdir+'input/parametric_initial')
	
def population_input(epos, PlotBox=True):
	assert epos.populationtype is 'model'
	
	for k, sg in enumerate(epos.groups):
		f, ax = plt.subplots()
		ax.set_title('Input Model Population {}'.format(sg['name']))
	
		ax.set_xlabel('Semi-Major Axis [au]')
		ax.set_ylabel(r'Planet Mass [M$_\bigoplus$]')

		ax.set_xlim(epos.mod_xlim)
		ax.set_ylim(epos.mod_ylim)

		ax.set_xscale('log')
		ax.set_yscale('log')

		ax.plot(sg['all_sma'], sg['all_mass'], color=clrs[k % 4], **fmt_symbol)
		
		#xx= np.linspace(*ax.get_xlim())
		#ax.plot(xx, (5./22.5)**(1.5) * xx)

		aux.save(plt, epos.plotdir+'input/model_planets.{}'.format(sg['name']))

	f, ax = plt.subplots()
	ax.set_title('Input Model Populations')

	ax.set_xlabel('Semi-Major Axis [au]')
	ax.set_ylabel(r'Planet Mass [M$_\bigoplus$]')

	ax.set_xlim(epos.mod_xlim)
	ax.set_ylim(epos.mod_ylim)

	ax.set_xscale('log')
	ax.set_yscale('log')

	for k, sg in enumerate(epos.groups):
		ax.plot(sg['all_sma'], sg['all_mass'], label=sg['name'], color=clrs[k % 4], **fmt_symbol)

	ax.legend(loc='upper right', shadow=False, prop={'size':14}, numpoints=1)
	aux.save(plt, epos.plotdir+'input/model_planets.all')

def population_input_diag(epos, PlotBox=True):
	assert epos.populationtype is 'model'
	
	f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
	f.subplots_adjust(hspace=0, wspace=0)
	# 1 2
	# 3 4
	
	ax1.set_title('Input Model Populations',loc='left') #  (<1 au)

	ax3.set_xlabel('Semi-Major Axis [au]')
	ax1.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
	
	ax2.set_xlabel('Occurrence/bin')
	ax3.set_ylabel('Occurrence/bin')
	
	ax1.set_xlim(epos.mod_xlim)
	ax1.set_ylim(epos.mod_ylim)

	ax1.set_xscale('log')
	ax1.set_yscale('log')

	for k, sg in enumerate(epos.groups):
		ax1.plot(sg['all_sma'], sg['all_mass'], ls='', marker='.', mew=0, ms=5.0, color=clrs[k % 4])
		
		xgrid= np.logspace(*np.log10(epos.mod_xlim))
		pdf= regression.sliding_window_log(sg['all_sma'], None, xgrid) /sg['n'] #, width=2. ) [sg['all_sma']<1]
		ax3.plot(xgrid, pdf, ls='-', marker='', color=clrs[k % 4],label=sg['name'])

		ygrid= np.logspace(*np.log10(epos.mod_ylim))
		pdf= regression.sliding_window_log(sg['all_mass'], None, ygrid) /sg['n'] #, width=2. )
		ax2.plot(pdf, ygrid, ls='-', marker='', color=clrs[k % 4])
	
	xmax= ax3.get_xlim()[-1]
	ymax= ax3.get_ylim()[-1]
	ax3.errorbar(xmax/1.5, ymax*0.9, xerr=(xmax/1.5*(1.-np.sqrt(1./2.)) ), color='k')
	
	xmax= ax2.get_xlim()[-1]
	ymin= ax2.get_ylim()[0]
	ax2.errorbar(xmax*0.9, ymin*1.5, yerr=(ymin*1.5*(1.-np.sqrt(1./2.)) ), color='k')
	
	#ax2.tick_params(labelbottom='on') # does not re-enable axis

	''' add observations? '''
	if epos.Observation and epos.DetectionEfficiency:
		from scipy import interpolate
		#feff_log= interpolate.RectBivariateSpline(epos.eff_xvar, epos.eff_yvar, epos.eff_2D_log) # broken
		#occ_inv_log= feff_log(epos.obs_xvar, epos.obs_yvar, grid=False)
		feff= interpolate.RectBivariateSpline(epos.eff_xvar, epos.eff_yvar, epos.eff_2D)
		occ_inv= feff(epos.obs_xvar, epos.obs_yvar, grid=False)	

		
#		subset= (occ_inv_log>-3) & (epos.obs_xvar<300) & (epos.obs_xvar>1) # some cut
#		occ= 1./(10.**(occ_inv_log[subset]))/epos.obs_nstars
		subset= (occ_inv>1e-3) & (epos.obs_xvar<300) & (epos.obs_xvar>1) # some cut, interpolation gives <0 -> occ=nan	
		occ= 1./occ_inv[subset]/epos.nstars

		sma= (epos.obs_xvar[subset]/365.)**(2./3.)
		#print np.min(occ)
		#print np.max(occ)
		
		xgrid= np.logspace(*np.log10(epos.mod_xlim))
		pdf= regression.sliding_window_log(sma, occ, xgrid) #, width=2. )
		#ax3.plot(xgrid, 3.*pdf, ls=':', marker='', color='k', label='Kepler x3')
		ax3.plot(xgrid, pdf, ls=':', marker='', color='k', label='Kepler')
	
	''' Legen instead of 4th plot'''
	ax3.legend(loc='center left', shadow=False, prop={'size':14}, numpoints=1,bbox_to_anchor=(1,0.5))
	ax4.axis('off')
	
	aux.save(plt, epos.plotdir+'input/model_planets.diag')

def inclination(epos):
	assert epos.populationtype is 'model'
	
	for k, sg in enumerate(epos.groups):
		f, ax = plt.subplots()
		ax.set_title('Input Inclination {}'.format(sg['name']))
	
		ax.set_xlabel('Semi-Major Axis [au]')
		#ax.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
		ax.set_ylabel('Inclination [degree]')

		ax.set_xlim(epos.mod_xlim)
		ax.set_ylim(0.1,90)

		ax.set_xscale('log')
		ax.set_yscale('log')

		ax.axhspan(1.0, 2.2, facecolor='0.5', alpha=0.5,lw=0.1)
		ax.axhline(1.6, color='0.5', ls='-')
		
		ax.axhline(np.median(sg['all_inc']), color=clrs[k % 4], ls='--')
		ax.plot(sg['all_sma'], sg['all_inc'], color=clrs[k % 4], **fmt_symbol)

		aux.save(plt, epos.plotdir+'input/inc.sma.{}'.format(sg['name']))

def Pratio(epos):
	# these plots don't seem very useful due to unfinished planet formation...
	assert epos.populationtype is 'model'
	
	for k,sg in enumerate(epos.groups):
		f, ax = plt.subplots()
		ax.set_title('Input Period Ratio {}'.format(sg['name']))
	
		ax.set_xlabel('Semi-Major Axis [au]')
		#ax.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
		ax.set_ylabel('P2/P1')

		ax.set_xlim(epos.mod_xlim)
		ax.set_ylim(0.1,100)

		ax.set_xscale('log')
		ax.set_yscale('log')

		ax.axhspan(1.2, 4.0, facecolor='0.5', alpha=0.5)
		ax.axhline(2.2, color='0.5', ls='-')
		
		ax.axhline(np.median(sg['all_Pratio']), color=clrs[k % 4], ls='--')
		#Pratio= np.isfinite()
		ax.plot(sg['all_sma'], sg['all_Pratio'], color=clrs[k % 4], **fmt_symbol)

		aux.save(plt, epos.plotdir+'input/Pratio.sma.{}'.format(sg['name']))


		f, ax = plt.subplots()
		ax.set_title('Input Period Ratio {}'.format(sg['name']))
	
		ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
		ax.set_ylabel('P2/P1')

		ax.set_xlim(epos.mod_ylim)
		ax.set_ylim(0.1,100)

		ax.set_xscale('log')
		ax.set_yscale('log')

		ax.axhspan(1.2, 4.0, facecolor='0.5', alpha=0.5)
		ax.axhline(2.2, color='0.5', ls='-')
		
		#Pratio= np.isfinite()
		ax.plot(sg['all_mass'], sg['all_Pratio'], color=clrs[k % 4], **fmt_symbol)

		aux.save(plt, epos.plotdir+'input/Pratio.mass.{}'.format(sg['name']))
	
''' output '''
def output(epos):

	print '\nPlotting output...'
	if epos.populationtype is 'parametric':
		periodradius(epos, Parametric=True, SNR=False)
		periodradius(epos, Parametric=True, SNR=True)
	elif epos.populationtype is 'model':
		massradius(epos, Input=True, MC=True)
		massradius(epos, Input=True, MC=True, Log=True)
		periodradius(epos, SNR=False)
		periodradius(epos, SNR=True)
		population_output_pdf(epos)
		if 'all_Pratio' in epos.groups[0]: 
			out_Pratio(epos)
			hist_Pratio(epos)

	else:
		assert False
	
	population_output_cdf(epos)
		
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
	ax.set_title('transiting planet sample'+suffix)
	_set_axes(ax, epos, Trim=True)
	if Parametric or len(epos.groups)==1:
		ax.plot(sim['P'], sim['R'], ls='', marker='.', mew=0, ms=5.0, color='k')
	else:
		for k, sg in enumerate(epos.groups):
			subset= sim['i sg']==k
			ax.plot(sim['P'][subset], sim['R'][subset], ls='', marker='.', mew=0, ms=5.0, color=clrs[k % 4], label=sg['name'])

	
	#ax.legend(loc='lower left', shadow=False, prop={'size':14}, numpoints=1)
	aux.save(plt, '{}output/periodradius{}'.format(epos.plotdir,fsuffix))

def out_Pratio(epos, SNR=True, Parametric=False):

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
	ax.set_title('Period Ratio Outer Planet'.format(suffix))
	_set_axes(ax, epos, Trim=True)

	ax.set_ylabel('P2/P1')
	ax.set_ylim(0.1,100)

	ax.set_yscale('log')

	ax.axhspan(1.2, 4.0, facecolor='0.5', alpha=0.5)
	ax.axhline(2.2, color='0.5', ls='-')
		
	for k, sg in enumerate(epos.groups):
		subset= sim['i sg']==k
		ax.plot(sim['P'][subset], sim['dP'][subset], ls='', marker='.', mew=0, ms=5.0, color=clrs[k % 4], label=sg['name'])

	
	#ax.legend(loc='lower left', shadow=False, prop={'size':14}, numpoints=1)
	aux.save(plt, '{}output/out_Pratio{}'.format(epos.plotdir,fsuffix))

def hist_Pratio(epos, SNR=True, Parametric=False):

	sim=epos.synthetic_survey
	
	for k, sg in enumerate(epos.groups):
		f, ax = plt.subplots()
		ax.set_title('Period Ratio Outer Planet {}'.format(sg['name']))

		ax.set_xlabel('P2/P1')
		ax.set_xlim(1,10)

		ax.axvspan(1.2, 4.0, facecolor='0.5', alpha=0.5)
		ax.axvline(2.2, color='0.5', ls='-',label='Kepler')
		
		bins=np.linspace(0,10,21)
		
		# SNR
		subset= sim['i sg']==k
		aN= np.isfinite(sim['dP'])
		ax.hist(sim['dP'][subset&aN], bins=bins, color=clrs[k % 4], normed=True, label='Simulated')
		
		# all
		aN= np.isfinite(sg['all_Pratio'])
		ax.hist(sg['all_Pratio'][aN], bins=bins, color='k', fill=False, normed=True, label='Input')
		
		ax.legend(loc='upper right', shadow=False, prop={'size':14}, numpoints=1)
		aux.save(plt, '{}output/hist_Pratio.{}'.format(epos.plotdir,sg['name']))


def population_output_pdf(epos):
	#assert epos.populationtype is 'model'
	
	f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
	f.subplots_adjust(hspace=0, wspace=0)
	# 1 2
	# 3 4
	
	ax1.set_title('Synthetic Model Populations',loc='left')

	ax3.set_xlabel('Orbital Period [days]')
	if epos.Radius:	ax1.set_ylabel(r'Planet Radius [R$_\bigoplus$]')
	else:			ax1.set_ylabel(r'Planet Mass [M$_\bigoplus$]')
	
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
			R= sim['R'][subset]
			ax1.plot(P,R, ls='', marker='.', mew=0, ms=5.0, color=clrs[k % 4])
				
			xgrid= np.logspace(*np.log10(epos.xtrim))
			pdf= regression.sliding_window_log(P, None, xgrid) #, width=2. )
			ax3.plot(xgrid, pdf, ls='-', marker='', color=clrs[k % 4],label='{} x{:.3f}'.format(sg['name'], sg['weight']))

			ygrid= np.logspace(*np.log10(epos.ytrim))
			pdf= regression.sliding_window_log(R, None, ygrid) #, width=2. )
			ax2.plot(pdf, ygrid, ls='-', marker='', color=clrs[k % 4])
	else:
		# top left panel (P,R)
		ax1.plot(sim['P'], sim['R'], ls='', marker='.', mew=0, ms=5.0, color='k')
		
	# PDF, all combined, 2 panels
	xgrid= np.logspace(*np.log10(epos.xtrim))
	pdf= regression.sliding_window_log(sim['P'], None, xgrid) #, width=2. )
	weight= np.sum([sg['weight'] for sg in epos.groups]) if epos.populationtype is 'model' else epos.fac
	ax3.plot(xgrid, pdf, ls='-', marker='', color='k',label='combined x {:.3f}'.format(weight))

	ygrid= np.logspace(*np.log10(epos.ytrim))
	pdf= regression.sliding_window_log(sim['R'], None, ygrid) #, width=2. )
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
	
	aux.save(plt, epos.plotdir+'output/pdf.diag')

def population_output_cdf(epos):
	
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	f.set_size_inches(9, 7) # default 7, 5

	''' 
	top left: synthetic obsservation
	'''
	sim= epos.synthetic_survey	
	ax1.set_title('Synthetic ({})'.format(sim['P'].size))
	_set_axes(ax1, epos, Trim=True)
	ax1.plot(sim['P'], sim['R'], ls='', marker='.', mew=0, ms=5.0, color='r', alpha=0.5)
		
	if epos.Zoom:
		ax1.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
			epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
		
	''' 
	Top right: observed population 
	'''
	ax2.set_title('Observed ({})'.format(epos.obs_zoom['x'].size))
	_set_axes(ax2, epos, Trim=True)

	ax2.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='b')		
	if epos.Zoom:
		ax2.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
			epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )

	''' 
	cdf orbital period 
	'''
	ax3.set_title('Period, p={:.3g}'.format(epos.gof['p xvar']))
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
	ax4.set_title('Radius, p={:.3g}'.format(epos.gof['p yvar']))
	if epos.Radius:	ax4.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
	else:			ax4.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	ax4.set_ylabel('CDF')
	ax4.set_xscale('log')
	ax4.set_ylim([-0.05,1.05])
	ax4.set_xlim(*epos.yzoom)

	#model histogram x
	R= sim['R zoom']
	ax4.plot(np.sort(R), np.arange(R.size, dtype=float)/R.size, ls='-', marker='', color='r')

	#obs histogram x
	R= epos.obs_zoom['y']
	ax4.plot(np.sort(R), np.arange(R.size, dtype=float)/R.size, ls='-', marker='', color='b')		
		
	f.tight_layout()
	aux.save(plt, epos.plotdir+'output/cdf.diag')

''' parameter study '''
def para(epos):
	print '\nPlotting parameter study...'
	
	# identify indices to dimensions with arrays
	k_dim= []
	for k, dim in enumerate(epos.para['grid']): 
		if len(dim)>1: k_dim.append(k) 
	
	all_dim= epos.para['prob'].shape
	ndim= len(k_dim) 
	indices= range(len(all_dim)) 
				
	print '  {} dimensional; {}'.format(ndim, k_dim)
	if ndim>1:
		for dim in k_dim:
			# sum over all axes _except_ dim :S
			index= tuple(indices[x] for x in indices if (x is not dim))
			print dim, index
			p_1d= np.sum(epos.para['prob'], axis=index)
			print p_1d.shape
			para_1d(epos, p_1d, epos.para['grid'][dim], 'dim{}'.format(dim))

def para_1d(epos, pdf, xgrid, fname):
	f, ax = plt.subplots()
	
	ax.set_title('Marginalized PDF')
	ax.set_xlabel(fname)
	ax.set_ylabel('Probability')

	ax.plot(xgrid, pdf, ls='-', marker='+', color='b')
	
	aux.save(plt, '{}para/1d.{}'.format(epos.plotdir,fname))		

''' helpers '''

def _set_axes(ax, epos, Trim=False, Eff=False):
	ax.set_xlabel('Orbital Period [days]')
	if epos.Radius:	ax.set_ylabel(r'Planet Radius [R$_\bigoplus$]')
	else:			ax.set_ylabel(r'Planet Mass [M$_\bigoplus$]')

	if Trim:
		ax.set_xlim(epos.xtrim)
		ax.set_ylim(epos.ytrim)
	elif Eff:
		ax.set_xlim(epos.eff_xlim)
		ax.set_ylim(epos.eff_ylim)
	else:
		ax.set_xlim(epos.obs_xlim)
		ax.set_ylim(epos.obs_ylim)

	ax.set_xscale('log')
	ax.set_yscale('log')

	#xticks=[3, 10, 30, 100, 300]
	#ax.set_xticks(xticks)
	#ax.set_xticklabels(xticks)

	#yticks=[0.5, 1, 2, 4]
	#ax.set_yticks(yticks)
	#ax.set_yticklabels(yticks)
	
	pass