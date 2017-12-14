import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar as clrbar
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec

import helpers
from EPOS import cgs
from EPOS.run import _pdf

def oneD(epos, PlotZoom=False, MCMC=False, Occ=False):
	if not epos.Range: epos.set_ranges()
	
	oneD_x(epos, PlotZoom=PlotZoom, MCMC=MCMC, Occ=Occ)
	oneD_y(epos, PlotZoom=PlotZoom, MCMC=MCMC, Occ=Occ)
	if epos.RV or epos.MassRadius:
		oneD_y(epos, PlotZoom=PlotZoom, MCMC=MCMC, PlotQ=True)

def oneD_x(epos, PlotZoom=False, MCMC=False, Occ=False):	
	# initial guess
	pps, _, pdf0_X, _= _pdf(epos, Init=True)
	if MCMC:
		# best-fit parameters
		pps, _, pdf_X, _= _pdf(epos)

	if Occ:
		fname= 'occurrence/posterior' if MCMC else 'occurrence/input'
	else:	
		fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	
	''' construct the posterior parameters '''
	if MCMC: plotsample= epos.samples[np.random.randint(len(epos.samples), size=100)]
	
	''' Orbital Period '''
	f, ax = plt.subplots()
	ax.set_title('Marginalized Distribution ({:.2f})'.format(pps))
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('Occurrence / {} P'.format(epos.plotpars['area']))
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.xtrim)
	ax.set_ylim([1e-3,1e1])
	
	if MCMC:
		for fpara in plotsample:
			_, _, xpdf, _= _pdf(epos, fpara=fpara)
			ax.plot(epos.MC_xvar, xpdf*epos.scale_x, color='b', alpha=0.1)
						
		ax.plot(epos.MC_xvar, pdf0_X*epos.scale_x, marker='',ls=':',color='k',
					label='starting guess')
		ax.plot(epos.MC_xvar, pdf_X*epos.scale_x, marker='',ls='-',color='k',
					label='best fit')
	else:
		if 'P break' in epos.fitpars.keys2d:
			ax.axvline(epos.fitpars.get('P break', Init=True), ls='-', color='gray')
		ax.plot(epos.MC_xvar, pdf0_X*epos.scale_x, marker='',ls='-',color='k')
    
	if epos.Zoom: 
		for zoom in epos.xzoom: ax.axvline(zoom, ls='--', color='k')
	
	if Occ:
		occbin= epos.occurrence['bin']
		ax.errorbar(occbin['xc'], occbin['occ']/occbin['dlnx'], yerr= occbin['err']/occbin['dlnx'])


	# lognormal inner disk edge?
	#from scipy.stats import norm
	#gauss= 2.*norm(loc=1., scale=0.4).pdf(np.log10(epos.MC_xvar))
	#ax.plot(epos.MC_xvar, gauss, ls='-', marker='', color='r')
    
	helpers.save(plt, epos.plotdir+fname+'_x')
	#print epos.plotdir+fname+'_x'

def oneD_y(epos, PlotZoom=False, MCMC=False, PlotQ=False, Occ=False):
	# initial guess
	pps, _, _, pdf0_Y= _pdf(epos, Init=True)
	if MCMC:
		# best-fit parameters
		pps, _, _, pdf_Y= _pdf(epos)

	fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	
	''' construct the posterior parameters '''
	if MCMC: plotsample= epos.samples[np.random.randint(len(epos.samples), size=100)]
	
	''' Planet Radius, Mass, or q '''
	# TODO: zip into previous block
	f, ax = plt.subplots()
	ax.set_title('Marginalized Distribution ({:.2f})'.format(pps))
	if PlotQ:
		ax.set_xlabel(r'Planet Mass Ratio')
		ax.set_ylabel('Occurrence / {} q d log a'.format(epos.plotpars['area']))
		area= 2./3.* np.log10(epos.MC_xvar[-1]/epos.MC_xvar[0])
		#print epos.MC_xvar[-1], epos.MC_xvar[0]
		print area
		yscale= epos.scale_in_y /area
	elif epos.RV or epos.MassRadius:	
		ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
		ax.set_ylabel('Occurrence / {} M'.format(epos.plotpars['area']))
		yscale= epos.scale_in_y
	else:
		ax.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
		ax.set_ylabel('Occurrence / {} R'.format(epos.plotpars['area']))
		yscale= epos.scale_in_y
	
	M_to_q= 1./(epos.Mstar*cgs.Msun/cgs.Mearth)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(np.array(epos.in_ytrim) * (M_to_q if PlotQ else 1.))
	ax.set_ylim([1e-3,1e1])
	
	yvar= epos.in_yvar * (M_to_q if PlotQ else 1.)
	
	if MCMC:
		for fpara in plotsample:
			_, _, _, ypdf= _pdf(epos, fpara=fpara)
			ax.plot(yvar, ypdf*yscale, color='b', alpha=0.1)
		ax.plot(yvar, pdf0_Y*yscale, marker='',ls=':',color='k')
		ax.plot(yvar, pdf_Y*yscale, marker='',ls='-',color='k')
	else:
		if 'R break' in epos.fitpars.keys2d:
			ax.axvline(epos.fitpars.get('R break', Init=True), ls='-', color='gray')

		ax.plot(yvar, pdf0_Y*yscale, marker='',ls='-',color='k')

	if epos.Zoom and not epos.MassRadius: 
		for zoom in epos.yzoom: ax.axvline(zoom, ls='--', color='k')
	
	if PlotQ and 'q Suzuki' in epos.plotpars:
		A, qbr, p1, p2= epos.plotpars['q Suzuki']
		ax.plot(yvar, A*(yvar/qbr)**np.where(yvar<qbr,p1,p2), color='g', \
			label='Suzuki+ 2016' )
		
		#ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
		ax.legend(loc='upper right')


	helpers.save(plt, epos.plotdir+fname+('_q' if PlotQ else '_y'))

def twoD(epos, PlotZoom=False, MCMC=False):
	
	# where does this go -> run.py
	assert epos.populationtype is 'parametric'
	if not epos.Range: epos.set_ranges()

	# pdf
	pps, pdf, _, _= _pdf(epos, Init= not MCMC)
	pdflog= np.log10(pdf*epos.scale) # in %
		
	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	
	ax.set_title('Occurrence [%] / d ln p d ln '+('M' if epos.MassRadius else 'R'))
	helpers.set_axes(ax, epos, Trim=True, In= epos.MassRadius)

	''' color scale? '''
	cmap='jet'
	vmin, vmax= -5, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	levels= np.linspace(vmin, vmax)
	ax.contourf(epos.X_in, epos.Y_in, pdflog, cmap=cmap, levels=levels)
	
	# colorbar?
	norm = Normalize(vmin=vmin, vmax=vmax)
	cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm, ticks=ticks,
                                orientation='vertical') # horizontal
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='out')
	
	fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	helpers.save(plt, epos.plotdir+fname)

def panels(epos, PlotZoom=False, MCMC=False):
	''' Initial distribution, panel layout'''
	gs = gridspec.GridSpec(2, 3,
                       width_ratios=[6, 20, 1],
                       height_ratios=[10, 4]
                       )
	f= plt.figure()
	f.subplots_adjust(wspace=0, hspace=0)
	
	ax = plt.subplot(gs[0, 1])	
	axb = plt.subplot(gs[0, 2])	

	axR = plt.subplot(gs[0, 0])
	axP = plt.subplot(gs[1, 1])

	# pdf
	pps, pdf, pdf_X, pdf_Y= _pdf(epos, Init= not MCMC)
	pdflog= np.log10(pdf*epos.scale_in) # in %
		
	ax.set_title('Occurrence [%] / d ln p d ln '+('M' if epos.MassRadius else 'R'))
	helpers.set_axes(ax, epos, Trim=True, In= epos.MassRadius)

	# Side panels
	axP.plot(epos.MC_xvar, pdf_X*epos.scale_x, marker='',ls='-',color='k')
	axR.plot(pdf_Y*epos.scale_in_y, epos.in_yvar, marker='',ls='-',color='k')

	helpers.set_axis_distance(axP, epos, Trim=True)
	helpers.set_axis_size(axR, epos, Trim=True, In= epos.MassRadius)

	axP.set_yscale('log')
	axP.set_ylim([2e-3,5])	
	axP.set_yticks([0.01,0.1,1])
	axP.set_yticklabels(['1%','10%','100%'])
	axP.yaxis.tick_right()
	axP.yaxis.set_ticks_position('both')
	axP.tick_params(axis='y', which='minor',left='off',right='off')

	axR.set_xscale('log')
	axR.set_xlim([2e-3,5])
	axR.set_xticks([0.01,0.1,1])
	axR.set_xticklabels(['1%','10%','100%'], rotation=70)
	axR.tick_params(axis='x', which='minor',top='off',bottom='off')
	axP.tick_params(axis='y', which='minor',left='off',right='off')

	''' color scale? '''
	cmap='jet'
	vmin, vmax= -5, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	levels= np.linspace(vmin, vmax)
	ax.contourf(epos.X_in, epos.Y_in, pdflog, cmap=cmap, levels=levels)
	
	# colorbar?
	norm = Normalize(vmin=vmin, vmax=vmax)
	cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm, ticks=ticks,
                                orientation='vertical') # horizontal
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='out')
	
	helpers.save(plt, epos.plotdir+'input/panels')
