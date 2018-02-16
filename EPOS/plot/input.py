import matplotlib, warnings
#warnings.filterwarnings("ignore", category=UserWarning, module='matplotlib')
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from EPOS import regression
import survey, parametric, multi, massradius, helpers, model

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def all(epos, color=None, imin=1e-2):
	print '\nPlotting input...'

	survey.observed(epos, PlotBox=False)
	survey.completeness(epos, PlotBox=False)
	if not epos.RV: survey.completeness(epos, PlotBox=False, Transit=True)
	
	if epos.Range and epos.Zoom:
		survey.observed(epos, PlotBox=True)
		survey.completeness(epos, PlotBox=True)
		if not epos.RV: survey.completeness(epos, PlotBox=True, Transit=True)
	
	if epos.Parametric:
		parametric.oneD(epos)
		parametric.twoD(epos)
		parametric.panels(epos)
	else:
		model.panels_mass(epos, color=color)
		if hasattr(epos, 'func'):
			model.panels_mass(epos, Population=True, color=color)
		if hasattr(epos, 'occurrence'):
			if 'planet' in epos.occurrence:	
				model.panels_radius(epos, Occurrence=True, color=color)
			if 'model' in epos.occurrence:
				model.panels_radius(epos, Observation=True, color=color)
		if 'inc' in epos.pfm:
			model.inclination(epos, color=color, imin=imin)
		if 'dP' in epos.pfm:
			model.periodratio(epos, color=color)
			
# 		input(epos, PlotBox=False)
# 		input_diag(epos, PlotBox=False)
# 		if 'all_Pratio' in epos.groups[0]:
# 			Pratio(epos)
# 		if 'all_inc' in epos.groups[0]:
# 			inclination(epos)

	if epos.Multi:
		if epos.Parametric:
			multi.periodradius(epos, MC=False)
			multi.periodradius(epos, MC=False, Nth=True)
			multi.multiplicity(epos, MC=False)
			multi.periodratio(epos, MC=False)
			multi.periodratio_cdf(epos, MC=False)
			multi.periodinner(epos, MC=False)
			multi.periodinner_cdf(epos, MC=False)
		else:
			pass
		
	if epos.MassRadius: 
		massradius.massradius(epos, MC=False)
		massradius.massradius(epos, MC=False, Log=True)


def input(epos, PlotBox=True):
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

		helpers.save(plt, epos.plotdir+'input/model_planets.{}'.format(sg['name']))

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
	helpers.save(plt, epos.plotdir+'input/model_planets.all')

def input_diag(epos, PlotBox=True):
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
	
	helpers.save(plt, epos.plotdir+'input/model_planets.diag')


