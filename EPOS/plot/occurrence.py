#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

import helpers

clrs= ['r','g','b','m'] # in epos.prep
#fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def all(epos):
	assert epos.Observation
	assert hasattr(epos, 'occurrence')
	colored(epos)
	binned(epos)

def colored(epos):
	
	f, ax = plt.subplots()
	ax.set_title('Occurrence per Planet')
	
	helpers.set_axes(ax, epos, Trim=True)
	
	''' color scale? '''
	cmap = plt.cm.get_cmap('PuRd')
	#cmap=cmap	
	#epos.occurrence['planet']
	
	ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='k')
	
	# colorbar?
	#cbar= f.colorbar(cs, ax=ax, shrink=0.95, ticks=ticks_log)
	#cbar.ax.set_yticklabels(labels)
	
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
	
