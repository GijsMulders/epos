import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colorbar as clrbar
from mpl_toolkits.mplot3d import Axes3D

from EPOS import regression
import helpers

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def readme():
	print '\nPlots the exoplanet survey: observed planets and completeness'
	print 

def observed(epos, PlotBox=True):
	assert epos.Observation
	f, ax = plt.subplots()
	ax.set_title('Observed Population')
	
	helpers.set_axes(ax, epos, Trim=True)
	
	ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='k')
	# add multis? 
	
	if PlotBox:
		fname='.box'
		assert epos.Range
		ax.add_patch(patches.Rectangle( (epos.xtrim[0],epos.ytrim[0]), 
			epos.xtrim[1]-epos.xtrim[0], epos.ytrim[1]-epos.ytrim[0],fill=False, zorder=1, ls='--') )
		if epos.Zoom:
			ax.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
				epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
	else: fname=''

	helpers.save(plt, epos.plotdir+'survey/planets'+fname)
	
def completeness(epos, PlotBox=False, Transit=False):
	assert epos.DetectionEfficiency

	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	#f.set_size_inches(7, 5)

	ax.set_title('Detection Efficiency [%]' if Transit else 'Completeness [%]')
	helpers.set_axes(ax, epos, Trim=False, Eff=True)
	
	toplot= epos.eff_2D if Transit else epos.completeness
	with np.errstate(divide='ignore'): log_completeness=  np.log10(toplot)
	
	''' color map and ticks'''
	cmap = 'PuRd'
	vmin, vmax= -5, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	levels=np.linspace(vmin, vmax, 256)

	cs= ax.contourf(epos.eff_xvar, epos.eff_yvar,log_completeness.T, cmap=cmap, levels=levels)

	cbar= f.colorbar(cs, cax=axb, ticks=ticks)
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='out')
	axb.set_title('%')
	
	''' Plot the zoom box or a few black contours'''	
	if PlotBox:
		fname='.box'
		assert epos.Range
		ax.add_patch(patches.Rectangle( (epos.xtrim[0],epos.ytrim[0]), 
			epos.xtrim[1]-epos.xtrim[0], epos.ytrim[1]-epos.ytrim[0],fill=False, zorder=1, ls='--') )
		if epos.Zoom:
			ax.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
				epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
	else: 
		if epos.RV:
			levels= [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
			labels= ['10','20','30','40','50','60','70','80','90']
			cs= ax.contour(epos.eff_xvar, epos.eff_yvar, toplot.T,\
					levels= levels, colors = 'k', linewidths=2.)
		else:
			levels= [1e-4, 1e-3, 0.01, 0.1, 0.5, 0.9] if Transit else 10.**ticks
			cs= ax.contour(epos.eff_xvar, epos.eff_yvar, toplot.T, 
				colors='k', levels=levels)
			fmt_percent= lambda x: '{:g} %'.format(100.*x)
			plt.clabel(cs, cs.levels, inline=True, fmt=fmt_percent)

		fname=''

	helpers.save(plt, epos.plotdir+'survey/'+('efficiency' if Transit else 'completeness')+ \
				fname)
