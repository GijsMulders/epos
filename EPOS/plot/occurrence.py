#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colorbar as clrbar
import numpy as np

import helpers, parametric
from EPOS.population import periodradius

clrs= ['r','g','b','m'] # in epos.prep
#fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def all(epos):
	assert epos.Observation
	if hasattr(epos, 'occurrence'):
		
		if 'planet' in epos.occurrence:
			colored(epos)

		if 'model' in epos.occurrence:
			model(epos)

		
		if 'bin' in epos.occurrence:
			colored(epos, Bins=True)
			#binned(epos)
	
			if 'eta0' in epos.occurrence['bin']:
				integrated(epos)
				if 'eta' in epos.occurrence['bin']:
					integrated(epos, MCMC=True)
					integrated(epos, MCMC=True,Planets=True)

		if 'xzoom' in epos.occurrence:
			if epos.Parametric:
				parametric.oneD(epos, Occ=True)
				if hasattr(epos, 'chain'):
					 parametric.oneD(epos, Occ=True, MCMC=True)
	else:
		print '\nNo occurrence to plot, did you run EPOS.occurrence.all()? \n'
	
def colored(epos, Bins=False):
	
	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	
	name= 'Survey Completeness'
	if epos.name in ['dr25_F','dr25_G','dr25_K','dr25_M','dr25_GK']: name+= ' ('+epos.name[5:]+')'
	ax.set_title(name)
	
	helpers.set_axes(ax, epos, Trim=True)
	#if epos.plot['']
	
	#ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='k')

	''' color scale? '''
	cmap='magma' # viridis, plasma, inferno, magma, spring, cool
	cmap='viridis'
	vmin, vmax= -4, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	clrs, norm= helpers.color_array(np.log10(epos.occurrence['planet']['completeness']),
		vmin=vmin,vmax=vmax, cmap=cmap)
	ax.scatter(epos.obs_xvar, epos.obs_yvar, color=clrs, s=4)
	
	# colorbar?
	cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm, ticks=ticks,
                                orientation='vertical') # horizontal
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='out')
	
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
			 
			size=16 if not 'textsize' in epos.plotpars else epos.plotpars['textsize'] 
				# 12 fit in box, 16 default
			ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'{:.1%}'.format(occ), va='top', 
				size=size)
			ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'\n$\pm${:.1f}'.format(
				occbin['err'][k]*100), va='top', size=size)
			ax.text(xbin[1]/xnudge,ybin[0]*ynudge,'n={}'.format(n), ha='right',
				size=size)
	
		helpers.save(plt, epos.plotdir+'occurrence/bins')
	else:
		helpers.save(plt, epos.plotdir+'occurrence/colored')
		
def integrated(epos, MCMC=False, Planets=False):
	
	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	
	sy= 'M' if (epos.MassRadius or epos.RV) else 'R'
	ax.set_title('Occurrence'+ (' (dln'+sy+' dlnP)' if MCMC else ' (Initial Guess)'))
	
	helpers.set_axes(ax, epos, Trim=True, In=epos.MassRadius)
	
	''' color scale? '''
	cmap='jet' # cool, spring
	vmin, vmax= -5, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	levels=np.linspace(vmin, vmax, 256)
	
	''' 2D pdf '''
	pps, pdf, _, _= periodradius(epos, Init=not MCMC)
	pdflog= np.log10(pdf) # in %
	cs= ax.contourf(epos.X_in, epos.Y_in, pdflog, cmap=cmap, levels=levels)
	cbar= f.colorbar(cs, cax=axb, ticks=ticks)
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='out')
	axb.set_title('%')
	
	''' integrated occurrence per bin'''
	occbin= epos.occurrence['bin']
	key = 'eta' if MCMC else 'eta0'
	for k, (xbin, ybin, n, inbin, occ) in enumerate(
			zip(occbin['x'],occbin['y in'],occbin['n'],occbin['i'], occbin[key])
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
		size=16 if not 'textsize' in epos.plotpars else epos.plotpars['textsize'] 
				# 12 fit in box, 16 default
		ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'{:.1%}'.format(occ), va='top',size=size)
		if MCMC:
			ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'\n +{:.1%}\n  -{:.1%}'.format(
				occbin['eta+'][k],occbin['eta-'][k]
				), va='top',size=size)

	''' overplot planets '''
	if Planets:
		ax.plot(epos.obs_xvar, epos.obs_yvar, 
			ls='', marker='.', mew=0, ms=5, alpha=1, color='k')

	fname= 'posterior' if MCMC else 'integrated'
	if Planets: fname+= '.planets'
	helpers.save(plt, epos.plotdir+'occurrence/'+fname)

def model(epos):
	
	f, ax = plt.subplots()
	
	name= 'Model Occurrence Rate, $\eta={:.2g}$'.format(epos.occurrence['model']['eta'])
	ax.set_title(name)
	
	helpers.set_axes(ax, epos, Trim=True)
	
	ax.plot(epos.pfm['P'], epos.pfm['R'], ls='', marker='o', mew=0, ms=4, color='C0', zorder=0)
	
	''' bins'''
	occbin= epos.occurrence['model']['bin']
	for k, (xbin, ybin, n, inbin, occ) in enumerate(
			zip(occbin['x'],occbin['y'],occbin['n'],occbin['i'], occbin['occ'])
			):
	
		# box
		ax.add_patch(patches.Rectangle( (xbin[0],ybin[0]), 
			xbin[1]-xbin[0], ybin[1]-ybin[0],
			fill=False, zorder=2, ls='-', color='k') )
	
		xnudge=1.01
		ynudge=1.02
		 
		size=16 if not 'textsize' in epos.plotpars else epos.plotpars['textsize'] 
			# 12 fit in box, 16 default
		ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'{:.1%}'.format(occ), va='top', 
			size=size)
		ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'\n$\pm${:.1f}'.format(
			occbin['err'][k]*100), va='top', size=size)
		ax.text(xbin[1]/xnudge,ybin[0]*ynudge,'n={}'.format(n), ha='right',
			size=size)

	helpers.save(plt, epos.plotdir+'occurrence/model')
