'''
Plots the exoplanet survey: observed planets and completeness
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colorbar as clrbar

import helpers

def all(epos):
	print '\nPlotting survey...'

	observed(epos, PlotBox=False)
	completeness(epos, PlotBox=False)
	if not epos.RV: 
		completeness(epos, PlotBox=False, Transit=True)
		if hasattr(epos, 'vetting'):
			completeness(epos, PlotBox=False, Transit=True, Vetting=False)
		
	if hasattr(epos, 'vetting'): 
		vetting(epos, PlotBox=False)
		completeness(epos, PlotBox=False, Vetting=False)
	
	if epos.Range and epos.Zoom:
		observed(epos, PlotBox=True)
		completeness(epos, PlotBox=True)
		if not epos.RV: 
			completeness(epos, PlotBox=True, Transit=True)
			if hasattr(epos, 'vetting'):
				completeness(epos, PlotBox=True, Transit=True)
		if hasattr(epos, 'vetting'): 
			vetting(epos, PlotBox=True)
			completeness(epos, PlotBox=True, Vetting=False)

	if hasattr(epos,'obs_score'):
		observed(epos, PlotBox=False, PlotScore=True)

def observed(epos, PlotBox=True, PlotScore=False):
	assert epos.Observation

	if PlotScore:
		f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
		f.subplots_adjust(wspace=0)
	else:
		f, ax = plt.subplots()
	ax.set_title('Observed Population')
	
	helpers.set_axes(ax, epos, Trim=epos.Range)
	
	if PlotScore:
		''' color scale? '''
		cmap='plasma' # viridis, plasma, inferno, magma, spring, cool
		vmin, vmax= 0, 1
		#ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
		clrs, norm= helpers.color_array(epos.obs_score,
			vmin=vmin,vmax=vmax, cmap=cmap)
		ax.scatter(epos.obs_xvar, epos.obs_yvar, color=clrs, s=3)

		# colorbar?
		cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm,
	                                orientation='vertical') # horizontal
		axb.tick_params(axis='y', direction='out')

		fname='.score'
	else:
		ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='k')
		fname=''
		# add multis? 
	
	if PlotBox:
		fname='.box'
		assert epos.Range
		ax.add_patch(patches.Rectangle( (epos.xtrim[0],epos.ytrim[0]), 
			epos.xtrim[1]-epos.xtrim[0], epos.ytrim[1]-epos.ytrim[0],fill=False, zorder=1, ls='--') )
		if epos.Zoom:
			ax.add_patch(patches.Rectangle( (epos.xzoom[0],epos.yzoom[0]), 
				epos.xzoom[1]-epos.xzoom[0], epos.yzoom[1]-epos.yzoom[0],fill=False, zorder=1) )

	helpers.save(plt, epos.plotdir+'survey/planets'+fname)
	
def completeness(epos, PlotBox=False, Transit=False, Vetting=True):
	assert epos.DetectionEfficiency

	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	#f.set_size_inches(7, 5)

	ax.set_title('Detection Efficiency' if Transit else 'Survey Completeness')
	helpers.set_axes(ax, epos, Trim=False, Eff=True)
	
	if Transit:
		toplot= epos.eff_2D
		if Vetting and hasattr(epos, 'vetting'): toplot*= epos.vetting
	else:
		toplot= epos.completeness if Vetting else epos.completeness_novet
	
	with np.errstate(divide='ignore'): log_completeness=  np.log10(toplot)
	
	''' color map and ticks'''
	cmap = 'YlOrBr' if Transit else 'PuRd'
	#cmap = 'rainbow' if Transit else 'PuRd'

	vmin, vmax= -4, 0
	if Transit: vmin=-3
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	levels=np.linspace(vmin, vmax, 256)

	cs= ax.contourf(epos.eff_xvar, epos.eff_yvar,log_completeness.T, cmap=cmap, levels=levels, vmin=vmin, vmax=vmax)

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
		if hasattr(epos, 'xtrim'):
			ax.set_xlim(*epos.xtrim)
			ax.set_ylim(*epos.ytrim)
		
		if epos.RV:
			levels= [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
			labels= ['10','20','30','40','50','60','70','80','90']
			cs= ax.contour(epos.eff_xvar, epos.eff_yvar, toplot.T,\
					levels= levels, colors = 'k', linewidths=2.)
		else:
			levels= [1e-3, 0.01, 0.1, 0.5, 0.9] if Transit else 10.**ticks
			cs= ax.contour(epos.eff_xvar, epos.eff_yvar, toplot.T, 
				colors='k', levels=levels)
			fmt_percent= lambda x: '{:g} %'.format(100.*x)
			plt.clabel(cs, cs.levels, inline=True, fmt=fmt_percent)

		fname=''

	if not Vetting: fname+= '.novet'
	helpers.save(plt, epos.plotdir+'survey/'+('efficiency' if Transit else 'completeness')+ \
				fname)

def vetting(epos, PlotBox=False):
	assert hasattr(epos, 'vetting')

	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	#f.set_size_inches(7, 5)

	ax.set_title('Vetting Efficiency')
	helpers.set_axes(ax, epos, Trim=False, Eff=True)
	

	''' color map and ticks'''
	cmap = 'plasma'
	levels=np.linspace(0,1, 256)
	ticks= np.array([0.0,0.25,0.5,0.75,1.0])
	ticks_percent=['0','25','50','75','100']

	cs= ax.contourf(epos.eff_xvar, epos.eff_yvar,epos.vetting.T, cmap=cmap, levels=levels)

	cbar= f.colorbar(cs, cax=axb, ticks=ticks)
	axb.set_yticklabels(ticks_percent)
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
	
		cs= ax.contour(epos.eff_xvar, epos.eff_yvar, epos.vetting.T, 
			colors='k', levels=ticks)
		fmt_percent= lambda x: '{:g} %'.format(100.*x)
		plt.clabel(cs, cs.levels, inline=True, fmt=fmt_percent)

		fname=''
		if hasattr(epos, 'xtrim'):		
			ax.set_xlim(*epos.xtrim)
			ax.set_ylim(*epos.ytrim)

	helpers.save(plt, epos.plotdir+'survey/vetting'+ fname)
