#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colorbar as clrbar
import numpy as np

import helpers

clrs= ['r','g','b','m'] # in epos.prep
#fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def all(epos):
	assert epos.Observation
	assert hasattr(epos, 'occurrence')
	colored(epos)
	colored(epos, Bins=True)
	#binned(epos)

def colored(epos, Bins=False):
	
	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	
	ax.set_title('Planet Completeness [%]') # if Bins
	
	helpers.set_axes(ax, epos, Trim=True)
	
	#ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='k')

	''' color scale? '''
	cmap='jet'
	vmin, vmax= -4, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	clrs, norm= helpers.color_array(np.log10(epos.occurrence['planet']['completeness']),
		vmin=vmin,vmax=vmax, cmap=cmap)
	ax.scatter(epos.obs_xvar, epos.obs_yvar, color=clrs)
	
	# colorbar?
	cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm, ticks=ticks,
                                orientation='vertical') # horizontal
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='both')
	
	''' bins?'''
	if Bins:
		occbin= epos.occurrence['bin']
		for k, (xbin, ybin, n, inbin, occ) in enumerate(
				zip(occbin['x'],occbin['y'],occbin['n'],occbin['i'], occbin['occ'])
				):
			clr= clrs[k%4]
		
			# colored dots
			#ax.plot(epos.obs_xvar[inbin], epos.obs_yvar[inbin], 
			#	ls='', marker='.', mew=0, ms=5.0, color=clr, zorder=1)
		
			# box
			ax.add_patch(patches.Rectangle( (xbin[0],ybin[0]), 
				xbin[1]-xbin[0], ybin[1]-ybin[0],
				fill=False, zorder=2, ls='-', color='k') )
		
			xnudge=1.01
			ynudge=1.02
			ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'{:.1%}'.format(occ), va='top')
			ax.text(xbin[1]/xnudge,ybin[0]*ynudge,'n={}'.format(n), ha='right')
	
		helpers.save(plt, epos.plotdir+'occurrence/bins')
	else:
		helpers.save(plt, epos.plotdir+'occurrence/colored')
	
def	binned(epos):
	f, ax = plt.subplots()
	ax.set_title('Occurrence per Planet')
	occbin= epos.occurrence['bin']
	
	helpers.set_axes(ax, epos, Trim=True)
	
	''' color scale? '''	
	ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='0.7', zorder=0)

	# loop over bins
	for k, (xbin, ybin, n, inbin, occ) in enumerate(
			zip(occbin['x'],occbin['y'],occbin['n'],occbin['i'], occbin['occ'])
		):
		clr= clrs[k%4]
		
		# colored dots
		#ax.plot(epos.obs_xvar[inbin], epos.obs_yvar[inbin], 
		#	ls='', marker='.', mew=0, ms=5.0, color=clr, zorder=1)
		
		# box
		ax.add_patch(patches.Rectangle( (xbin[0],ybin[0]), 
			xbin[1]-xbin[0], ybin[1]-ybin[0],fill=False, zorder=2, ls='-', color='k') )
		
		ax.text(xbin[0],ybin[1],'{:.1%}'.format(occ), va='top')
		ax.text(xbin[1],ybin[0],'n={}'.format(n), ha='right')
	
	helpers.save(plt, epos.plotdir+'occurrence/binned')
	
