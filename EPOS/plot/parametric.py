import numpy as np
import matplotlib.pyplot as plt
import helpers

def oneD(epos, PlotZoom=False, MCMC=False):
	
	# where does this go? parametric.py?
	assert epos.populationtype is 'parametric'
	
	if not epos.Range: epos.set_ranges()
	
	# initial guess
	pdf0= epos.func(epos.X, epos.Y, *epos.p0)
	pdf0_X, pdf0_Y= np.sum(pdf0, axis=1), np.sum(pdf0, axis=0)
	pps_x, pps_y=  np.sum(pdf0_X), np.sum(pdf0_Y)
	scale_x, scale_y= pdf0_X.size, pdf0_Y.size
	if MCMC:
		# best-fit parameters
		pdf= epos.func(epos.X, epos.Y, *epos.pfit)
		pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
		pps_x, pps_y=  np.sum(pdf_X), np.sum(pdf_Y)
		#scale_x, scale_y= pdf_X.size, pdf_Y.size		

	
	fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	
	''' Orbital Period '''
	f, ax = plt.subplots()
	ax.set_title('Marginalized Distribution ({:.2f})'.format(pps_x))
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('Occurrence')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.xtrim)
	ax.set_ylim([1e-3,1e1])
	if MCMC:
		for para in epos.samples[np.random.randint(len(epos.samples), size=100)]:
			xpdf= np.sum(epos.func(epos.X, epos.Y,*para), axis=1)
			ax.plot(epos.MC_xvar, xpdf*scale_x, color='b', alpha=0.1)
		ax.plot(epos.MC_xvar, pdf0_X*scale_x, marker='',ls=':',color='k')
		ax.plot(epos.MC_xvar, pdf_X*scale_x, marker='',ls='-',color='k')
	else:
		ax.plot(epos.MC_xvar, pdf0_X*scale_x, marker='',ls='-',color='k')
    
	helpers.save(plt, epos.plotdir+fname+'_x')

	''' Planet Radius'''
	# TODO: zip into previous block
	f, ax = plt.subplots()
	ax.set_title('Marginalized Distribution ({:.2f})'.format(pps_y))
	if epos.Radius:	ax.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
	else:			ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	ax.set_ylabel('Occurrence')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.ytrim)
	ax.set_ylim([1e-3,1e1])
	if MCMC:
		for para in epos.samples[np.random.randint(len(epos.samples), size=100)]:
			ypdf= np.sum(epos.func(epos.X, epos.Y,*para), axis=0)
			ax.plot(epos.MC_yvar, ypdf*scale_y, color='b', alpha=0.1)
		ax.plot(epos.MC_yvar, pdf0_Y*scale_y, marker='',ls=':',color='k')
		ax.plot(epos.MC_yvar, pdf_Y*scale_y, marker='',ls='-',color='k')
	else:
		ax.plot(epos.MC_yvar, pdf0_Y*scale_y, marker='',ls='-',color='k')

	helpers.save(plt, epos.plotdir+fname+'_y')

def twoD(epos, PlotZoom=False, MCMC=False):
	
	# where does this go -> run.py
	assert epos.populationtype is 'parametric'
	if not epos.Range: epos.set_ranges()
	
	pdf= epos.func(epos.X, epos.Y, *(epos.p0 if not MCMC else epos.pfit))
	scale= pdf.size
	pdflog= np.log10(pdf*scale) # in %
		
	f, ax = plt.subplots()
	ax.set_title('Probability Density Function ({:.2f})'.format(np.sum(pdf)))
	helpers.set_axes(ax, epos, Trim=True)

	# contour with labels (not straightforward)
	# missing: vmin, vmx, a log scale, contour linestyles
	labels=['0.001', '0.01', '0.1','1.0','10','100']
	n=len(labels)-1
	levels_log= np.linspace(-n,2,10*n+1)
	ticks_log= np.linspace(-n,2,n+1)
	ticks= np.logspace(-n,2, n+1)
	
	cmap = plt.cm.get_cmap('jet')
	cs= ax.contourf(epos.X, epos.Y, pdflog, cmap=cmap, levels=levels_log) #vmin=-3., vmax=0.)	# vmax keywords don't work on colorbar
	cbar= f.colorbar(cs, ax=ax, shrink=0.95, ticks=ticks_log)
	cbar.ax.set_yticklabels(labels)

	fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	helpers.save(plt, epos.plotdir+fname)

