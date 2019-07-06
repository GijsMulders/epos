import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar as clrbar
from matplotlib.colors import Normalize

import helpers
from EPOS import cgs
from EPOS.population import periodradius

def oneD(epos, PlotZoom=False, MCMC=False, Occ=False):
	if not epos.Range: epos.set_ranges()
	
	oneD_x(epos, PlotZoom=PlotZoom, MCMC=MCMC, Occ=Occ)
	oneD_y(epos, PlotZoom=PlotZoom, MCMC=MCMC, Occ=Occ)
	if epos.MassRadius:
		# works with occ??
		oneD_y(epos, PlotZoom=PlotZoom, MCMC=MCMC, PlotQ=True)

def oneD_x(epos, PlotZoom=False, MCMC=False, Occ=False, Log=True):	

	if Occ:
		fname= 'occurrence/posterior' if MCMC else 'occurrence/input'
		ybin= epos.yzoom
		unit= r'$M_\bigoplus$' if epos.RV else r'$R_\bigoplus$'
		title= r'Planet Occurrence ({1:.1f}-{2:.0f} {0})'.format(unit, *epos.yzoom)
	else:	
		fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
		ybin=None
		title= r'Marginalized Distribution ({:.1f}-{:.0f} $R_\bigoplus$)'.format(*epos.ytrim)
	
	if not Log: 
		fname+= '.linear'

	# initial guess
	pps, _, pdf0_X, _= periodradius(epos, Init=True, ybin=ybin)
	if MCMC:
		# best-fit parameters
		pps, _, pdf_X, _= periodradius(epos, ybin=ybin)
	
	''' construct the posterior parameters '''
	if MCMC: plotsample= epos.samples[np.random.randint(len(epos.samples), size=100)]
		
	''' Orbital Period '''
	f, ax = plt.subplots()
	ax.set_title(title)
	ax.set_ylabel('Occurrence / {}P'.format(epos.plotpars['area']))

	#ax.set_xlabel('Orbital Period [days]')
	#ax.set_xscale('log')
	#ax.set_xlim(epos.xtrim)
	helpers.set_axis_distance(ax, epos, Trim=True)
	
	if Log:
		ax.set_yscale('log')
		if 'occrange' in epos.plotpars:
			ax.set_ylim(epos.plotpars['occrange'])
		else:	
			ax.set_ylim([1e-3,1e1])
	else:
		ax.set_ylim([0,0.45])
	
	if MCMC:
		for fpara in plotsample:
			_, _, xpdf, _= periodradius(epos, fpara=fpara, ybin=ybin)
			ax.plot(epos.MC_xvar, xpdf, color='b', alpha=0.1)
						
		ax.plot(epos.MC_xvar, pdf0_X, marker='',ls=':',color='k',
					label='Starting Guess')
		ax.plot(epos.MC_xvar, pdf_X, marker='',ls='-',color='k',
					label='Best Fit')
	else:
		if 'P break' in epos.fitpars.keys2d:
			ax.axvline(epos.fitpars.get('P break', Init=True), ls='-', color='gray')
		ax.plot(epos.MC_xvar, pdf0_X, marker='',ls='-',color='k')
    
	# plot posterior excluding low detection regions (arbitrary 3000 planets assumed)
	#if not (epos.RV or epos.MassRadius):
	if False:
		cens= np.where(epos.f_det<1./3000.,0,1.)
		_, _, cens_pdf_X, _ = periodradius(epos, ybin=ybin, fdet=cens)
		ax.plot(epos.MC_xvar, cens_pdf_X, marker='',ls='-',color='green', label='biased')


	if epos.Zoom: 
		for zoom in epos.xzoom: ax.axvline(zoom, ls='--', color='k')
	
	if Occ:
		occbin= epos.occurrence['yzoom']
		ax.errorbar(occbin['xc'], occbin['occ']/occbin['dlnx'], 
			yerr= occbin['err']/occbin['dlnx'], color='r', marker='_', ls='', capsize=3) #,capthick=2

		if 'occ_MLE' in occbin:
			ax.errorbar(occbin['xc'], occbin['occ_MLE']/occbin['dlnx'], 
			yerr= occbin['err_MLE']/occbin['dlnx'],color='g', 
				label='MLE', marker='x', ls=':', capsize=3)

	# lognormal inner disk edge?
	#from scipy.stats import norm
	#gauss= 2.*norm(loc=1., scale=0.4).pdf(np.log10(epos.MC_xvar))
	#ax.plot(epos.MC_xvar, gauss, ls='-', marker='', color='r')

	if not Log:
		if 'Pinner all' in epos.plotpars:
			xx=np.geomspace(20,400)
			a,b,c= epos.plotpars['Pinner all']
			ax.plot(xx, a* (xx/b)**c, marker='', ls='-', color='g', 
				label='All Planets')
		if 'Pinner ylim' in epos.plotpars:
			ax.set_ylim(epos.plotpars['Pinner ylim'])

	if MCMC:
		ax.legend(loc='upper left')
    
	helpers.save(plt, epos.plotdir+fname+'_x')
	#print epos.plotdir+fname+'_x'

def oneD_y(epos, PlotZoom=False, MCMC=False, PlotQ=False, Occ=False, Convert=False):
	if Occ:
		fname= 'occurrence/posterior' if MCMC else 'occurrence/input'
		xbin= epos.xzoom
		title= 'Planet occurrence ({:.2g}-{:.0f} days)'.format(*epos.xzoom)

	else:	
		fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
		xbin=None
		title= 'Marginalized Distribution ({:.2g}-{:.0f} days)'.format(*epos.xtrim)
	
	# initial guess
	pps, _, _, pdf0_Y= periodradius(epos, Init=True, xbin=xbin, Convert=Convert)
	if MCMC:
		# best-fit parameters
		pps, _, _, pdf_Y= periodradius(epos, xbin=xbin, Convert=Convert)
	elif Convert:
		pps_in, _, _, pdf0_Y_in= periodradius(epos, Init=True, xbin=xbin, Convert=False)
	
	''' construct the posterior parameters '''
	if MCMC: plotsample= epos.samples[np.random.randint(len(epos.samples), size=100)]
	
	''' Planet Radius, Mass, or q '''
	M_to_q= 1./(epos.Mstar*cgs.Msun/cgs.Mearth)
	# TODO: zip into previous block
	f, ax = plt.subplots()
	ax.set_title(title)
	if PlotQ:
		ax.set_xlabel(r'Planet Mass Ratio')
		ax.set_ylabel('Occurrence / {}q dloga'.format(epos.plotpars['area']))
		area= 2./3.* np.log10(epos.MC_xvar[-1]/epos.MC_xvar[0])
		yscale= 1. /area
		ytrim= epos.in_ytrim
		yvar= epos.in_yvar * M_to_q
	elif epos.MassRadius or (epos.Msini and not Convert):	
		ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
		ax.set_ylabel('Occurrence / {}M'.format(epos.plotpars['area']))
		yscale= 1.
		ytrim= epos.in_ytrim
		yvar= epos.in_yvar
	elif epos.RV:
		ax.set_xlabel(r'Planet Minimum Mass [M$_\bigoplus$]')
		ax.set_ylabel('Occurrence / {}Msini'.format(epos.plotpars['area']))
		yscale= 1.
		ytrim= epos.ytrim
		yvar= epos.MC_yvar
	else:
		ax.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
		ax.set_ylabel('Occurrence / {}R'.format(epos.plotpars['area']))
		yscale= 1.
		ytrim= epos.ytrim
		yvar= epos.MC_yvar
	
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(np.array(ytrim) * (M_to_q if PlotQ else 1.))
	if 'occrange' in epos.plotpars:
		ax.set_ylim(epos.plotpars['occrange'])
	else:	
		ax.set_ylim([1e-3,1e1])
	
	if MCMC:
		for fpara in plotsample:
			_, _, _, ypdf= periodradius(epos, fpara=fpara, xbin=xbin, Convert=Convert)
			ax.plot(yvar, ypdf*yscale, color='b', alpha=0.1)
		ax.plot(yvar, pdf0_Y*yscale, marker='',ls=':',color='k', label='starting guess')
		ax.plot(yvar, pdf_Y*yscale, marker='',ls='-',color='k', label='best-fit')
	else:
		if 'R break' in epos.fitpars.keys2d:
			ax.axvline(epos.fitpars.get('R break', Init=True), ls='-', color='gray')

		ax.plot(yvar, pdf0_Y*yscale, marker='',ls='-',color='k', label='M sin i' if Convert else None)
		if Convert:
			ax.plot(epos.in_yvar, pdf0_Y_in*yscale, marker='',ls='--',color='k', label='Intrinsic')

	# plot posterior excluding low detection regions (arbitrary 3000 planets assumed)
	#if not (epos.RV or epos.MassRadius):
	if False:
		cens= np.where(epos.f_det<1./3000.,0,1.)
		_, _, _, cens_pdf_Y= periodradius(epos, xbin=xbin, fdet=cens)
		ax.plot(yvar, cens_pdf_Y*yscale, marker='',ls='-',color='green', label='biased')


	if epos.Zoom and not epos.MassRadius: 
		for zoom in epos.yzoom: ax.axvline(zoom, ls='--', color='k')
	
	if PlotQ and 'q Suzuki' in epos.plotpars:
		A, qbr, p1, p2= epos.plotpars['q Suzuki']
		ax.plot(yvar, A*(yvar/qbr)**np.where(yvar<qbr,p1,p2), color='g', \
			label='Suzuki+ 2016' )
		
		#ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
	
	if Occ:
		occbin= epos.occurrence['xzoom']
		ax.errorbar(occbin['yc'], occbin['occ']/occbin['dlny'], 
			yerr= occbin['err']/occbin['dlny'], color='r', marker='_', ls='', capsize=3)

		if 'occ_MLE' in occbin:
			ax.errorbar(occbin['yc'], occbin['occ_MLE']/occbin['dlny'], 
				yerr= occbin['err_MLE']/occbin['dlny'], color='g', 
				label='MLE', marker='x', ls=':', capsize=3)
		
			
	# 	for y, occ, err in zip(occbin['yc'], occbin['occ']/occbin['dlny'], 
	# 						occbin['err']/occbin['dlny']):
	# 		print '{:.3g} {:.3g} {:.3g}'.format(y, occ, err)

	#if Occ or MCMC:
	if MCMC or (Occ and Convert):
		ax.legend(loc='upper right')

	helpers.save(plt, epos.plotdir+fname+('_q' if PlotQ else '_y')+('_convert' if Convert else''))

def twoD(epos, PlotZoom=False, MCMC=False, Convert=False):
	
	# where does this go -> run.py
	assert epos.Parametric
	if not epos.Range: epos.set_ranges()

	# pdf
	pps, pdf, _, _= periodradius(epos, Init= not MCMC, Convert=Convert)
	pdflog= np.log10(pdf) # in %
		
	f, (ax, axb) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[20, 1]})
	f.subplots_adjust(wspace=0)
	
	ax.set_title('Occurrence [%] / d ln p d ln '+('M' if (epos.MassRadius or epos.RV) else 'R'))
	helpers.set_axes(ax, epos, Trim=True, In= (epos.MassRadius or epos.RV) and not Convert)

	''' color scale? '''
	cmap='jet'
	vmin, vmax= -5, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	levels= np.linspace(vmin, vmax)
	if Convert:
		ax.contourf(epos.X, epos.Y, pdflog, cmap=cmap, levels=levels)
	else:
		ax.contourf(epos.X_in, epos.Y_in, pdflog, cmap=cmap, levels=levels)
	
	# colorbar?
	norm = Normalize(vmin=vmin, vmax=vmax)
	cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm, ticks=ticks,
                                orientation='vertical') # horizontal
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='out')
	
	fname= 'mcmc/posterior' if MCMC else 'input/parametric_initial'
	helpers.save(plt, epos.plotdir+fname+('_convert' if Convert else''))

def panels(epos, PlotZoom=False, MCMC=False, Convert=False):
	''' Initial distribution, panel layout'''
	f, (ax, axb, axR, axP)= helpers.make_panels_clrbar(plt)

	# pdf
	pps, pdf, pdf_X, pdf_Y= periodradius(epos, Init= not MCMC, Convert=Convert)
	pdflog= np.log10(pdf) # in %
		
	ax.set_title('Planet Occurrence / dlnP dln'+('M' if epos.MassRadius else 'R'))
	helpers.set_axes(ax, epos, Trim=True, In= (epos.MassRadius or epos.RV) and not Convert)

	# Side panels
	axP.plot(epos.MC_xvar, pdf_X, marker='',ls='-',color='k')
	axR.plot(pdf_Y, epos.MC_yvar if Convert else epos.in_yvar, marker='',ls='-',color='k')

	#helpers.set_axis_distance(axP, epos, Trim=True)
	#helpers.set_axis_size(axR, epos, Trim=True, In= epos.MassRadius)
	axP.set_xlabel(ax.get_xlabel())
	axR.set_ylabel(ax.get_ylabel())

	axP.set_yscale('log')
	axP.set_ylim([2e-3,5])	
	axP.set_yticks([0.01,0.1,1])
	axP.set_yticklabels(['1%','10%','100%'])
	axP.yaxis.tick_right()
	axP.yaxis.set_ticks_position('both')
	axP.tick_params(axis='y', which='minor',left=False,right=False)

	axR.set_xscale('log')
	axR.set_xlim([2e-3,5])
	axR.set_xticks([0.01,0.1,1])
	axR.set_xticklabels(['1%','10%','100%'], rotation=70)
	axR.tick_params(axis='x', which='minor',top=False,bottom=False)
	axP.tick_params(axis='y', which='minor',left=False,right=False)
	axP.tick_params(axis='y', which='major',left=False,right=True)

	''' color scale? '''
	cmap='jet'
	vmin, vmax= -5, 0
	ticks=np.linspace(vmin, vmax, (vmax-vmin)+1)
	levels= np.linspace(vmin, vmax)
	if Convert:
		ax.contourf(epos.X, epos.Y, pdflog, cmap=cmap, levels=levels)
	else:
		ax.contourf(epos.X_in, epos.Y_in, pdflog, cmap=cmap, levels=levels)
	
	# colorbar?
	norm = Normalize(vmin=vmin, vmax=vmax)
	cb1 = clrbar.ColorbarBase(axb, cmap=cmap, norm=norm, ticks=ticks,
                                orientation='vertical') # horizontal
	axb.set_yticklabels(100*10.**ticks)
	axb.tick_params(axis='y', direction='out')
	axb.set_title('%')
	
	helpers.save(plt, epos.plotdir+'input/panels'+('_convert' if Convert else''))
