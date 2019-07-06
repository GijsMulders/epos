#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colorbar as clrbar
import matplotlib.colors
import numpy as np

import helpers, parametric
from EPOS.population import periodradius

clrs= ['r','g','b','m'] # in epos.prep
#fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

def all(epos, color=None, alpha_fac=None):
	assert epos.Observation
	if hasattr(epos, 'occurrence'):
		
		if 'planet' in epos.occurrence:
			colored(epos)

		if 'model' in epos.occurrence:
			model(epos, color=color)
			if alpha_fac is not None:
				model(epos, color=color, alpha_fac=alpha_fac)
			#if Fade:
			model(epos, color=color, Gradient=True)

		if 'poly' in epos.occurrence:
			colored(epos, Poly=True)
			if 'model' in epos.occurrence:
				model(epos, color=color, alpha_fac=alpha_fac, Poly=True)
				if 'labels' in epos.occurrence['poly']:
					# only callable with models right now
					poly_only(epos)

		
		if 'bin' in epos.occurrence:
			colored(epos, Bins=True)
			if 'model' in epos.occurrence:
				model(epos, color=color, alpha_fac=alpha_fac, Bins=True)
	
			if 'eta0' in epos.occurrence['bin']:
				integrated(epos)
				if 'eta' in epos.occurrence['bin']:
					integrated(epos, MCMC=True)
					integrated(epos, MCMC=True,Planets=True)

		if 'xzoom' in epos.occurrence:
			if epos.Parametric:
				parametric.oneD(epos, Occ=True)
				if not epos.MonteCarlo and epos.Msini:
					parametric.oneD_y(epos, Occ=True, Convert=True)

				if hasattr(epos, 'chain'):
					 parametric.oneD(epos, Occ=True, MCMC=True)
					 if epos.Msini:
					 	parametric.oneD_y(epos, Occ=True, MCMC=True, Convert=True)
	else:
		print '\nNo occurrence to plot, did you run EPOS.occurrence.all()? \n'
	
def colored(epos, Bins=False, Poly=False):
	
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
	elif Poly:
		occpoly= epos.occurrence['poly']
		for k, (xc, yc, coords, n, inbin, occ, err) in enumerate(
				zip(occpoly['xc'],occpoly['yc'],occpoly['coords'],
					occpoly['n'],occpoly['i'], occpoly['occ'], occpoly['err'])
				):
				
			# box
			ax.add_patch(matplotlib.patches.Polygon(coords,
				fill=False, zorder=2, ls='-', color='k') )
			 
			size=16 if not 'textsize' in epos.plotpars else epos.plotpars['textsize'] 
				# 12 fit in box, 16 default
			ax.text(xc,yc,'{:.1%}\n$\pm${:.1f}'.format(occ, err*100), ha='center', va='center', 
				size=size)
			#ax.text(xbin[1]/xnudge,ybin[0]*ynudge,'n={}'.format(n), ha='right',
			#	size=size)
	
		helpers.save(plt, epos.plotdir+'occurrence/poly')
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

def model(epos, color='C0', alpha_fac=None, Bins=False, Poly=False, Gradient=False):
	
	f, ax = plt.subplots()
	
	name= '{}, $\eta={:.2g}$'.format(epos.name, epos.occurrence['model']['eta'])
	ax.set_title(name)
	
	helpers.set_axes(ax, epos, Trim=True)
	
	# set transparency / color gradient
	if Gradient:
		suffix= '.gradient'

		weigths= epos.occurrence['model']['completeness']
		cmin, cmax= 0.001, 0.1

		weigths= np.maximum(np.minimum(weigths,cmax), cmin)

		cmap='copper_r'
		#ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
		clrs, norm= helpers.color_array(np.log10(weigths),
			vmin=np.log10(cmin),vmax=np.log10(cmax), cmap=cmap)
		
		ax.scatter(epos.pfm['P'], epos.pfm['R'], 
			marker='o', s=13, lw=0, color=clrs,zorder=0)

		# colorbar?
		# cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm, ticks=ticks,
		#                             orientation='vertical') # horizontal
		# axb.set_yticklabels(100*10.**ticks)
		# axb.tick_params(axis='y', direction='out')


	elif alpha_fac is not None:
		suffix= '.alpha'

		weigths= epos.occurrence['model']['completeness']*alpha_fac #*epos.nstars
		alpha= np.maximum(np.minimum(weigths,1.), 0.0) # 0.2?

		if True:
			# color issues with  to_rgba_array
			clr_rgba = np.empty((len(alpha), 4), float)
			for i, a in enumerate(alpha):
				clr_rgba[i] = matplotlib.colors.to_rgba(color, a)

		else:
			clr= np.full_like(weigths,color,dtype=str)
			clr_rgba= matplotlib.colors.to_rgba_array(clr) # alpha
			#print clr_rgba[0,:]
			clr_rgba[:,3]= alpha

		ax.scatter(epos.pfm['P'], epos.pfm['R'], 
			marker='o', s=13, lw=0, color=clr_rgba,zorder=0)
	else:
		suffix=''
		clr= matplotlib.colors.to_rgba(color)
		ax.plot(epos.pfm['P'], epos.pfm['R'], ls='', marker='o', mew=0, ms=4, 
			color=clr, zorder=0)
	
	''' bins'''
	if Bins:
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
		helpers.save(plt, epos.plotdir+'occurrence/model_bins'+suffix)
	elif Poly:
		occpoly= epos.occurrence['model']['poly']
		for k, (xc, yc, coords, n, inbin, occ, err) in enumerate(
				zip(occpoly['xc'],occpoly['yc'],occpoly['coords'],
					occpoly['n'],occpoly['i'], occpoly['occ'], occpoly['err'])
				):
				
			# box
			ax.add_patch(matplotlib.patches.Polygon(coords,
				fill=False, zorder=2, ls='-', color='k') )
			 
			size=16 if not 'textsize' in epos.plotpars else epos.plotpars['textsize'] 
				# 12 fit in box, 16 default
			ax.text(xc,yc,'{:.1%}\n$\pm${:.1%}'.format(occ, err), ha='center', va='center', 
				size=size)
		helpers.save(plt, epos.plotdir+'occurrence/model_poly'+suffix)
	else:
		helpers.save(plt, epos.plotdir+'occurrence/model'+suffix)

def poly_only(epos):
	
	f, ax = plt.subplots()
	
	ax.set_title('Planet Classes')
	
	helpers.set_axes(ax, epos, Trim=True)
	
	# coordinates are from model routine
	occpoly= epos.occurrence['model']['poly']

	for k, (xc, yc, coords, label) in enumerate(
			zip(occpoly['xc'],occpoly['yc'],occpoly['coords'],
				epos.occurrence['poly']['labels'])
			):
			
		# box
		ax.add_patch(matplotlib.patches.Polygon(coords,
			fill=False, zorder=2, ls='-', color='k') )
		 
		size=16 if not 'textsize' in epos.plotpars else epos.plotpars['textsize'] 
			# 12 fit in box, 16 default
		ax.text(xc,yc,label, ha='center', va='center', 
			size=size)

	helpers.save(plt, epos.plotdir+'occurrence/poly_only')

