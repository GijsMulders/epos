import numpy as np
import matplotlib.pyplot as plt
import helpers
from EPOS.run import _pdf

def oneD(epos, PlotZoom=False, MCMC=False):
	
	# where does this go? parametric.py?
	assert epos.populationtype is 'parametric'
	
	if not epos.Range: epos.set_ranges()
	
	# initial guess
	pps, _, pdf0_X, pdf0_Y= _pdf(epos, Init=True)
	if MCMC:
		# best-fit parameters
		pps, _, pdf_X, pdf_Y= _pdf(epos)

	fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	
	''' construct the posterior parameters '''
	if MCMC: plotsample= epos.samples[np.random.randint(len(epos.samples), size=100)]
	
	''' Orbital Period '''
	f, ax = plt.subplots()
	ax.set_title('Marginalized Distribution ({:.2f})'.format(pps))
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('Occurrence / d ln P')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.xtrim)
	ax.set_ylim([1e-3,1e1])
	if MCMC:
		for fpara in plotsample:
			_, _, xpdf, _= _pdf(epos, fpara=fpara)
			ax.plot(epos.MC_xvar, xpdf*epos.scale_x, color='b', alpha=0.1)
						
		ax.plot(epos.MC_xvar, pdf0_X*epos.scale_x, marker='',ls=':',color='k')
		ax.plot(epos.MC_xvar, pdf_X*epos.scale_x, marker='',ls='-',color='k')
	else:
		ax.plot(epos.MC_xvar, pdf0_X*epos.scale_x, marker='',ls='-',color='k')
    
	if epos.Zoom: 
		for zoom in epos.xzoom: ax.axvline(zoom, ls='--', color='k')

	# lognormal inner disk edge?
	#from scipy.stats import norm
	#gauss= 2.*norm(loc=1., scale=0.4).pdf(np.log10(epos.MC_xvar))
	#ax.plot(epos.MC_xvar, gauss, ls='-', marker='', color='r')
    
	helpers.save(plt, epos.plotdir+fname+'_x')

	''' Planet Radius'''
	# TODO: zip into previous block
	f, ax = plt.subplots()
	ax.set_title('Marginalized Distribution ({:.2f})'.format(pps))
	if epos.RV:	ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	else:		ax.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
	ax.set_ylabel('Occurrence / d ln R')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.ytrim)
	ax.set_ylim([1e-3,1e1])
	if MCMC:
		for fpara in plotsample:
			_, _, _, ypdf= _pdf(epos, fpara=fpara)
			ax.plot(epos.MC_yvar, ypdf*epos.scale_y, color='b', alpha=0.1)
		ax.plot(epos.MC_yvar, pdf0_Y*epos.scale_y, marker='',ls=':',color='k')
		ax.plot(epos.MC_yvar, pdf_Y*epos.scale_y, marker='',ls='-',color='k')
	else:
		ax.plot(epos.MC_yvar, pdf0_Y*epos.scale_y, marker='',ls='-',color='k')

	if epos.Zoom: 
		for zoom in epos.yzoom: ax.axvline(zoom, ls='--', color='k')

	helpers.save(plt, epos.plotdir+fname+'_y')

def twoD(epos, PlotZoom=False, MCMC=False):
	
	# where does this go -> run.py
	assert epos.populationtype is 'parametric'
	if not epos.Range: epos.set_ranges()

	# pdf
	pps, pdf, _, _= _pdf(epos, Init= not MCMC)
	pdflog= np.log10(pdf*epos.scale) # in %
		
	f, ax = plt.subplots()
	ax.set_title('Occurrence [%] / d ln p d ln R')
	helpers.set_axes(ax, epos, Trim=True)

	# contour with labels (not straightforward)
	# missing: vmin, vmx, a log scale, contour linestyles
	levels= np.linspace(-5,0,51)
	ticks= np.linspace(-5,0,6)
	labels=['0.001', '0.01', '0.1','1.0','10','100']
# 	n=len(labels)-1
# 	levels_log= np.linspace(-n,-1,10*n+1)
# 	ticks_log= np.linspace(-n,-1,n+1)
# 	print ticks_log
# 	ticks= np.logspace(-n,-1, n+1)
	
	
	cmap = plt.cm.get_cmap('jet')
	cs= ax.contourf(epos.X, epos.Y, pdflog, cmap=cmap, levels=levels) #vmin=-3., vmax=0.)	# vmax keywords don't work on colorbar
	cbar= f.colorbar(cs, ax=ax, shrink=0.95, ticks=ticks)
	cbar.ax.set_yticklabels(labels)

	fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	helpers.save(plt, epos.plotdir+fname)

