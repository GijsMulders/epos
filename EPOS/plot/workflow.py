import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from EPOS.plot import helpers


def periodradius_2x2(epos):
	f, axes = plt.subplots(2, 2, sharex=True, sharey=True)
	((ax1, ax2), (ax3, ax4)) = axes
	f.set_size_inches(9, 7) # default 7, 5
	f.subplots_adjust(wspace=0.5, hspace=0.5)

	if False:
		#f.suptitle(epos.name)
		f.suptitle("Compare Planet Formation Models to Observed Exoplanets with epos")
	else:
		words = ['Compare','Planet Formation Models','to',
			'Observed Exoplanets','with','epos']
		colors = ['b', '0.5', 'k', 'r', 'k', 'g']

		helpers.rainbow_text(0.05, 0.95, words, colors, size=16, f=f, fudge=1.5)
		#helpers.rainbow_text(0.05, 1.1, words, colors, size=18, ax=axes[0,0])

	helpers.set_axes(axes[1,0], epos, Trim=True)

	''' Top left: model '''
	ax= axes[0,0]
	ax.set_title('Formation Model: {}'.format(epos.name))
	pfm= epos.pfm
	ax.plot(pfm['P'],pfm['R'], ls='', marker='.', ms=5.0, color='0.5')

	''' Top right: synthetic obs'''	
	ax= axes[0,1]
	ax.set_title('Simulated Observations')

	if epos.MonteCarlo:
		sim= epos.synthetic_survey
		ax.plot(sim['P'], sim['Y'], ls='', marker='.', mew=0, ms=5.0, color='r', alpha=0.5)
	else:
		levels= np.linspace(0,np.max(sim['pdf']))		
		ax.contourf(epos.MC_xvar, epos.MC_yvar, sim['pdf'].T, cmap='Reds', levels=levels)

	''' Bottom left: occurrnce rates'''
	ax= axes[1,0]
	ax.set_title('Occurrence Rates')
	
	occbin= epos.occurrence['bin']
	maxocc= np.max(occbin['occ'])
	for k, (xbin, ybin, n, inbin, occ) in enumerate(
			zip(occbin['x'],occbin['y'],occbin['n'],occbin['i'], occbin['occ'])
		):
		
		# box
		ax.add_patch(patches.Rectangle( (xbin[0],ybin[0]), 
				xbin[1]-xbin[0], ybin[1]-ybin[0],
				fill=True, ls='-', 
				fc=cm.binary(occ/maxocc) ))
		
		#xnudge=1.01
		#ynudge=1.02	
		#ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'{:.1%}'.format(occ), va='top', 
		#	size=8)

	# colored dots
	#ax.plot(epos.obs_xvar, epos.obs_yvar, 
	#	ls='', marker='.', mew=0, ms=2.0, color='k')

	''' Bottom right: observations '''
	ax= axes[1,1]
	ax.set_title('Observed Planets')
	ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='C3')		

	''' Draw arrows between plots'''
	props=dict(transform=f.transFigure, arrowstyle='simple', connectionstyle='arc3',
    	alpha = 0.3, fc = 'g', mutation_scale = 80.)
	xl, xr= 0.44, 0.57
	yt, yb= 0.75, 0.25
	arrow1 = patches.FancyArrowPatch( (xl, yt), (xr, yt), **props)
	arrow2 = patches.FancyArrowPatch( (xr, yb), (xl, yb), **props)

	f.text((xl+xr)/2,yt,'Bias',color='k', 
		ha='center', va='center', transform=f.transFigure)
	f.text((xl+xr)/2,yb,'Debias',color='k', 
		ha='center', va='center', transform=f.transFigure)

	# and top-bottom
	props['connectionstyle']= 'arc3,rad=0.9'
	props['mutation_scale']= 30.
	props['fc']='b'
	xl, xr= 0.27, 0.75
	xw= 0.03
	yt, yb= 0.55, 0.46

	arrow3 = patches.FancyArrowPatch( (xl-xw, yt), (xl-xw, yb), **props)
	arrow4 = patches.FancyArrowPatch( (xl+xw, yb), (xl+xw, yt), **props)
	for xt in [xl, xr]:
		f.text(xt,(yt+yb)/2,'Compare',color='b', 
			ha='center', va='center', transform=f.transFigure)

	arrow5 = patches.FancyArrowPatch( (xr-xw, yt), (xr-xw, yb), **props)
	arrow6 = patches.FancyArrowPatch( (xr+xw, yb), (xr+xw, yt), **props)

	f.patches.extend([arrow1,arrow2,arrow3,arrow4,arrow5,arrow6])


	#f.tight_layout()
	helpers.save(plt, epos.plotdir+'workflow/arrows_2x2')

def periodradius(epos, color='C4', alpha=1):
	gs = gridspec.GridSpec(nrows=9, ncols=3)
	f=plt.figure()
	f.set_size_inches(13, 7) # default 7, 5
	f.subplots_adjust(wspace=0.5, hspace=0.0)

	ax1= f.add_subplot(gs[3:6,0])
	ax2= f.add_subplot(gs[0:3,1], sharex=ax1, sharey=ax1)
	ax3= f.add_subplot(gs[6:9,1], sharex=ax1, sharey=ax1)
	ax4= f.add_subplot(gs[3:6,2], sharex=ax1, sharey=ax1)

	if True:
		#f.suptitle(epos.name)
		f.suptitle("Compare Planet Formation Models to Observed Exoplanets with epos", 
			bbox=dict(boxstyle='round', fc='w', ec='k'))
	else:
		words = ['Compare','Planet Formation Models','to',
			'Observed Exoplanets','with','epos']
		colors = ['C0', color, '0.5', 'C3', '0.5', 'C2']

		helpers.rainbow_text(0.1, 0.95, words, colors, size=20, f=f, fudge=1.5)
		#helpers.rainbow_text(0.05, 1.1, words, colors, size=18, ax=axes[0,0])

	

	''' Left: model '''
	ax=ax1
	helpers.set_axes(ax, epos, Trim=True)
	ax.set_xlabel('')

	ax.set_title('Formation Model: {}'.format(epos.name))
	pfm= epos.pfm
	ax.plot(pfm['P'],pfm['R'], ls='', marker='.', ms=5.0, color=color, alpha=alpha)

	''' Top middle: synthetic obs'''	
	ax=ax2
	ax.set_title('Simulated Observations')
	if epos.MonteCarlo:
		sim= epos.synthetic_survey
		ax.plot(sim['P'], sim['Y'], ls='', marker='.', mew=0, ms=5.0, color='C2', alpha=0.5)
	else:
		levels= np.linspace(0,np.max(sim['pdf']))		
		ax.contourf(epos.MC_xvar, epos.MC_yvar, sim['pdf'].T, cmap='Greens', levels=levels)

	''' Bottom middle: occurrnce rates'''
	ax=ax3
	ax.set_title('Occurrence Rates')
	ax.set_xlabel('Orbital Period [days]')
	
	occbin= epos.occurrence['bin']
	maxocc= np.max(occbin['occ'])
	for k, (xbin, ybin, n, inbin, occ) in enumerate(
			zip(occbin['x'],occbin['y'],occbin['n'],occbin['i'], occbin['occ'])
		):
		
		# box
		ax.add_patch(patches.Rectangle( (xbin[0],ybin[0]), 
				xbin[1]-xbin[0], ybin[1]-ybin[0],
				fill=True, ls='-', 
				fc=cm.Greens(occ/maxocc) ))
		
		#xnudge=1.01
		#ynudge=1.02	
		#ax.text(xbin[0]*xnudge,ybin[1]/ynudge,'{:.1%}'.format(occ), va='top', 
		#	size=8)

	# colored dots
	#ax.plot(epos.obs_xvar, epos.obs_yvar, 
	#	ls='', marker='.', mew=0, ms=2.0, color='k')

	''' right: observations '''
	ax= ax4
	ax.set_title('Observed Planets')
	ax.plot(epos.obs_xvar, epos.obs_yvar, ls='', marker='.', mew=0, ms=5.0, color='C3')		

	''' 
	Draw horizontal bars
	'''
	xl=0.07
	xw, dy= 0.85, 0.08
	yb, yt= 0.22, 0.75
	props=dict(transform=f.transFigure,alpha=0.3, zorder=-10)
	bar_fw= patches.Rectangle((xl, yt-dy), xw, 2*dy, color='C1', **props)
	bar_occ= patches.Rectangle((xl, yb-dy), xw, 2*dy, color='C6', **props)

	bbox=dict(boxstyle='round', fc='w', ec='k')
	props= dict(rotation=0, ha='left', va='center', 
		transform=f.transFigure, bbox=bbox)

	f.patches.extend([bar_fw, bar_occ])
	f.text(xl-0.03,yt,'Forward Model',color='k', **props)
	f.text(xl-0.03,yb,'Inverse Model',color='k', **props)

	''' Draw arrows between plots'''
	props=dict(transform=f.transFigure, arrowstyle='simple', 
		connectionstyle='arc3,rad=-0.3',
    	alpha = 0.3, fc = 'g', mutation_scale = 50.)
	xl, xr= 0.3, 0.62
	yt, yb= 0.7, 0.2
	xw, yw= 0.1, 0.1
	arrow1 = patches.FancyArrowPatch( (xl, yt), (xl+xw, yt+yw), **props)
	arrow2 = patches.FancyArrowPatch( (xr+xw, yb+yw), (xr, yb), **props)

	f.text(xl-0.05,yt+0.06,'Apply Survey\nDetection Bias',color='g', 
		ha='center', va='center', transform=f.transFigure)
	f.text(xr+xw+0.07,yb,'Account for\nSurvey Completeness',color='g', 
		ha='center', va='center', transform=f.transFigure)

	# and top-bottom
	props['connectionstyle']= 'arc3,rad=0.0'
	props['mutation_scale']= 30.
	props['fc']='b'
	#props['shape']= 'full'
	#props['arrowstyle']='darrow'

	xl, xr= 0.33, 0.63
	yt, yb= 0.7, 0.25

	xw, yw= 0.005, 0.01
	dy=0.04
	xs,ys= 0.05,0.03

	arrow3 = patches.FancyArrowPatch( (xl-xw, yb+dy), (xl-xw+xs, yb-dy), **props)
	arrow4 = patches.FancyArrowPatch( (xl+xw+xs, yb-dy+ys), (xl+xw, yb+dy+ys), **props)
	#for x,y in zip([xl-0.04, xr+xs+0.04],[yb-dy,yt+dy]):
	#	f.text(x,y,'Compare',color='b', 
	#		ha='center', va='center', transform=f.transFigure)
	f.text(xr+xs+0.1,yt+dy,'Compare Distribution\nof Detections',color='b', 
		ha='center', va='center', transform=f.transFigure)

	f.text(xl-0.07,yb-dy,'Compare Intrinsic\n Planet Population',color='b', 
		ha='center', va='center', transform=f.transFigure)

	arrow5 = patches.FancyArrowPatch( (xr-xw, yt+dy), (xr-xw+xs, yt-dy), **props)
	arrow6 = patches.FancyArrowPatch( (xr+xw+xs, yt-dy+ys), (xr+xw, yt+dy+ys), **props)

	f.patches.extend([arrow1,arrow2,arrow3,arrow4,arrow5,arrow6])

	''' Draw crossed out arrow '''
	props['mutation_scale']= 50.
	props['fc']='w'
	props['zorder']=-1
	xc, yc= 0.5,0.5
	xw= 0.15
	if False:
		arrow_left = patches.FancyArrowPatch( (xc, yc), (xc-xw, yc), **props)
		arrow_right = patches.FancyArrowPatch( (xc, yc), (xc+xw, yc), **props)
	else:
		yw=0.02
		arrow_left = patches.FancyArrowPatch( (xc+xw-0.01, yc+yw), (xc-xw, yc+yw), **props)
		arrow_right = patches.FancyArrowPatch( (xc-xw+0.01, yc-yw), (xc+xw, yc-yw), **props)

	f.patches.extend([arrow_left,arrow_right])

	dx, dy= 0.05,0.05
	redcross = plt.scatter(xc, yc, s=2000, c='red', 
		transform=f.transFigure, marker='x', lw=7, clip_on=False)
	ax.add_artist(redcross)

	#f.tight_layout()
	helpers.save(plt, epos.plotdir+'workflow/arrows')

def multi(epos, color='C4', nplanets=10, PlotObs=True):

	pars_in=['dM', 'Pin', 'dP', 'n', 'inc'] # inc / ecc
	pars_out= ['dR', 'Pin', 'dP', 'n']

	#pars= ['mass', 'dM', 'Pin', 'dP', 'inc','n']
	titles= dict(inc='Mutual\nInclination', dP='Period\n Ratio', 
	Pin='Inner\n Planet', dM='Mass\n Ratio', dR='Radius\n Ratio', n='# Planets\nper System')
	units=dict(inc='[$\\degree$]',dP='',Pin='[day]',dM='',dR='',n='')
	xlim= dict(inc=[1e-2,90], dP=[1,10], Pin=[1,100], dM=[0.1,10], dR=[0.1,10], 
		n=[0,nplanets])
	nbins= dict(mass=30, sma=30, inc=30, dP=20, Pin=20, dM=20, dR=20 , n=nplanets+1)
	xticks= dict(mass=[0.1,1,10,30], inc=[1e-2,0.1,10,90], 
		dP=[1,2,4,10], Pin=[1,10,100], dM=[0.1,1,10], dR=[0.1,1,10] ,
		n=[nn for nn in [0,5,10,20] if nn<=nplanets])


	# axis layout
	gs = gridspec.GridSpec(nrows=4, ncols=5, top=0.8) # top to make space for title
	f=plt.figure()

	f.set_size_inches(13, 7) # default 7, 5
	f.subplots_adjust(wspace=0.2, hspace=0.0)

	ax_in= [f.add_subplot(gs[0,i]) for i in range(len(pars_in))]

	ax_out= [f.add_subplot(gs[3,i], sharex=ax_in[i]) 
		for i in range(len(pars_out)) ]

	f.suptitle('Compare Model Planetary Systems to Observed Exoplanets with epos',
		bbox=dict(boxstyle='round', fc='w', ec='k'))

	''' Input '''
	model={}
	model['Pin']= epos.pfm['Pin']
	model['dP']= epos.pfm['dP'][epos.pfm['dP']>1]
	model['dM']= epos.pfm['dM'][epos.pfm['dP']>1]
	model['n']= epos.pfm['n']
	model['inc']= epos.pfm['inc']

	for ax, par in zip(ax_in, pars_in):
		ax.set_title(titles[par])
		ax.set_xlim(xlim[par])
		ax.tick_params(left=False, labelleft=False)
		ax.tick_params(which='minor', bottom=False)
		[s(xticks[par]) for s in [ax.xaxis.set_ticks, ax.xaxis.set_ticklabels]]

		if par not in ['n']:
			ax.set_xscale('log')
			makespace= np.geomspace
		else:
			makespace= np.linspace

		bins= makespace(num=nbins[par], *xlim[par])

		# data
		ax.hist(model[par], bins=bins, color=color)

	''' Output '''
	zoomed={}
	zoomed['n']= epos.synthetic_survey['multi']['n']
	zoomed['dP']= epos.synthetic_survey['multi']['Pratio']
	zoomed['Pin']= epos.synthetic_survey['multi']['Pinner']
	zoomed['dR']= epos.synthetic_survey['multi']['Rratio']

	for ax, par in zip(ax_out, pars_out):
		ax.set_title(titles[par])
		#ax.tick_params(which='minor', bottom=False)
		[s(xticks[par]) for s in [ax.xaxis.set_ticks, ax.xaxis.set_ticklabels]]
		ax.tick_params(which='both',left=False, labelleft=False)

		if par not in ['n']:
			#ax.set_xscale('log')
			makespace= np.geomspace
		else:
			makespace= np.linspace
			ax.set_yscale('log')


		bins= makespace(num=nbins[par], *xlim[par])

		# data
		ax.hist(zoomed[par], bins=bins, color=color)

	''' Observations '''
	if PlotObs:
		obs_zoomed={}
		obs_zoomed['n']= epos.obs_zoom['multi']['n']
		obs_zoomed['dP']= epos.obs_zoom['multi']['Pratio']
		obs_zoomed['Pin']= epos.obs_zoom['multi']['Pinner']
		obs_zoomed['dR']= epos.obs_zoom['multi']['Rratio']

		for ax, par in zip(ax_out, pars_out):
			ax2= ax.twinx() # twinx re-adds ticks to orignal axis?
			ax.tick_params(left=False, labelleft=False, right=False, labelright=False)

			if par not in ['n']:
				makespace= np.geomspace
			else:
				makespace= np.linspace
				ax2.set_yscale('log')
			bins= makespace(num=nbins[par], *xlim[par])

			# data
			ax2.hist(obs_zoomed[par], bins=bins, color='C3', histtype='step',ls='--')
			ax2.tick_params(left=False, labelleft=False, right=False, labelright=False)


	''' Arrows '''
	props=dict(arrowstyle='simple', connectionstyle='arc3,rad=0',
		alpha = 1, fc = 'w', ec='k', mutation_scale=30, zorder=0)
	
	xc= 0.5
	for ax, par in zip(ax_out, pars_out):
		arrow = patches.FancyArrowPatch( (xc, 2.5), (xc, 1.5), transform=ax.transAxes, **props)
		f.patches.append(arrow)
	
	if len(pars_in) > len(pars_out):
		arrow = patches.FancyArrowPatch((xc, -0.5),(xc-0.7, -1.5),
			transform=ax_in[-1].transAxes, **props)
		f.patches.append(arrow)


	''' Text '''
	bbox=dict(boxstyle='round', fc='w', ec='k')
	props= dict(ha='center', va='center', bbox=bbox)

	f.text(-0.5,0.5,'Formation\nModel',color=color, transform=ax_in[0].transAxes, **props)

	f.text(-0.5,0.7,'Observable\nPlanets',color=color, transform=ax_out[0].transAxes, **props)
	if PlotObs:
		f.text(-0.5,0.2,'Kepler\nObservations',color='r', 
			ha='center', va='center', transform=ax_out[0].transAxes)
	
	f.text(0,-1,'Simulate\nTransit\nSurvey',color='k', 
			ha='center', va='center', transform=ax_in[0].transAxes)

	helpers.save(plt, epos.plotdir+'workflow/multi', bbox_inches=None)
