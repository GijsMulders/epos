#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
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
	f, ax = plt.subplots()
	ax.set_title('Detection Efficiency [%]' if Transit else 'Completeness [%]')
	helpers.set_axes(ax, epos, Trim=False, Eff=True)

	# contour with labels (not straightforward)
	# missing: vmin, vmx, a log scale, contour linestyles
	labels=['0.001','0.01', '0.1','1.0','10','100']
	n=len(labels)-1
	levels_log= np.linspace(-n,0,10*n+1)
	ticks_log= np.linspace(-n,0,n+1)
	ticks= np.logspace(-n,0, n+1)
	
	toplot= epos.eff_2D if Transit else epos.completeness
	
	cmap = plt.cm.get_cmap('PuRd')
	with np.errstate(divide='ignore'): log_completeness=  np.log10(toplot)
	cs= ax.contourf(epos.eff_xvar, epos.eff_yvar,log_completeness.T, cmap=cmap, levels=levels_log) #vmin=-3., vmax=0.)	# vmax keywords don't work on colorbar
	cbar= f.colorbar(cs, ax=ax, shrink=0.95, ticks=ticks_log)
	cbar.ax.set_yticklabels(labels)
		
	if PlotBox:
		fname='.box'
		assert epos.Range
		ax.add_patch(patches.Rectangle( (epos.xtrim[0],epos.ytrim[0]), 
			epos.xtrim[1]-epos.xtrim[0], epos.ytrim[1]-epos.ytrim[0],fill=False, zorder=1, ls='--') )
		if epos.Zoom:
			ax.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
				epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )
	else: 
		#ax.contour(epos.eff_xvar, epos.eff_yvar, epos.eff_2D.T,levels= ticks, colors = ['k','k','k','k'], linewidths=3.)
		if epos.RV:
			levels= [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
			labels= ['10','20','30','40','50','60','70','80','90']
			cs= ax.contour(epos.eff_xvar, epos.eff_yvar, toplot.T,\
					levels= levels, colors = 'k', linewidths=2.)
		else:
			for tick, label in zip(ticks,labels):
				cs= ax.contour(epos.eff_xvar, epos.eff_yvar, toplot.T,\
					levels= [tick], colors = ['k'], linewidths=3.)
				plt.clabel(cs, inline=1, fmt=label+'%%')

		fname=''

	helpers.save(plt, epos.plotdir+'survey/'+('efficiency' if Transit else 'completeness')+ \
				fname)
