import numpy as np
import matplotlib.pyplot as plt
import helpers

def oneD(epos, PlotZoom=False):
	# these are NOT integrated over other axes (TODO) -> use epos.pdf
	
	# where does this go? parametric.py?
	assert epos.populationtype is 'parametric'
	if not epos.Range: epos.set_ranges()
	ngrid= 100.
	xgrid= np.logspace(np.log10(epos.xtrim[0]),np.log10(epos.xtrim[1]),ngrid)
	ygrid= np.logspace(np.log10(epos.ytrim[0]),np.log10(epos.ytrim[1]),ngrid)
	parax= epos.norm0* epos.xfunc(xgrid, *epos.xp0)
	paray= epos.norm0* epos.yfunc(ygrid, *epos.yp0)
	
	f, ax = plt.subplots()
	ax.set_title('Starting guess planet population')
	ax.set_xlabel('Orbital Period [days]')
	ax.set_ylabel('Occurrence []')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.xtrim)
	ax.set_ylim([1e-3,1e2])
	ax.plot(xgrid, 100.* parax, marker='',ls=':',color='k')
	helpers.save(plt, epos.plotdir+'input/parametric_initial_x')

	f, ax = plt.subplots()
	ax.set_title('Starting guess planet population')
	if epos.Radius:	ax.set_xlabel(r'Planet Radius [R$_\bigoplus$]')
	else:			ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	ax.set_ylabel('Occurrence []')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(epos.ytrim)
	ax.set_ylim([1e-3,1e2])
	ax.plot(ygrid, 100.* paray, marker='',ls=':',color='k')
	helpers.save(plt, epos.plotdir+'input/parametric_initial_y')

def twoD(epos, PlotZoom=False):
	
	# where does this go -> run.py
	assert epos.populationtype is 'parametric'
	if not epos.Range: epos.set_ranges()
	ngrid= 100. # use native MCgrid instead
	xgrid= np.logspace(np.log10(epos.xtrim[0]),np.log10(epos.xtrim[1]),ngrid)
	ygrid= np.logspace(np.log10(epos.ytrim[0]),np.log10(epos.ytrim[1]),ngrid)
	X,Y= np.meshgrid(xgrid, ygrid)
	para2D= epos.norm0* epos.xfunc(X, *epos.xp0)* epos.yfunc(Y,*epos.yp0)
	paralog= np.log10(para2D)
	
	f, ax = plt.subplots()
	ax.set_title('Starting guess planet population')
	helpers.set_axes(ax, epos, Trim=True)

	# contour with labels (not straightforward)
	# missing: vmin, vmx, a log scale, contour linestyles
	labels=['0.001','0.01', '0.1','1.0','10','100']
	n=len(labels)-1
	levels_log= np.linspace(-n,0,10*n+1)
	ticks_log= np.linspace(-n,0,n+1)
	ticks= np.logspace(-n,0, n+1)
	
	cmap = plt.cm.get_cmap('jet')
	cs= ax.contourf(X, Y, paralog, cmap=cmap, levels=levels_log) #vmin=-3., vmax=0.)	# vmax keywords don't work on colorbar
	cbar= f.colorbar(cs, ax=ax, shrink=0.95, ticks=ticks_log)
	cbar.ax.set_yticklabels(labels)

	helpers.save(plt, epos.plotdir+'input/parametric_initial')
