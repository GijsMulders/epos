import numpy as np
import matplotlib.pyplot as plt
import helpers
from EPOS import cgs
	
def massradius(epos, MC=False, Log=False, color='C1', Mlim=[0,20], Rlim=[0,7], rho=None):
	assert epos.MassRadius
	# plot R(M)
	f, ax = plt.subplots()
	ax.set_title('mass-radius relation + dispersion')
	ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	ax.set_ylabel(r'Planet Radius [R$_\bigoplus$]')

	if Log:
		ax.set_xlim(epos.masslimits)
		ax.set_ylim(0.1, 30.) # 7.	
		ax.set_xscale('log')
		ax.set_yscale('log')
	else:
		ax.set_xlim(*Mlim) # 20.
		ax.set_ylim(*Rlim) # 7.
	
	# MC data
	if MC:
		ss=epos.synthetic_survey
		ax.plot(ss['M'], ss['R'], ls='', marker='.', mew=0, ms=5.0, color=color, zorder=1)
		
	xM= np.logspace(*np.log10(ax.get_xlim())) if Log else np.linspace(*ax.get_xlim()) #np.max(tr['M']))
	xR, dispersion= epos.MR(xM)
	#ax.plot(xM, epos.MR['fRock'](xM), ls='-', marker='', color='r', label='Rocky constraint', zorder=0)
	#xGas= epos.MR['fGas'](xM)
	#dispersion= epos.MR['dispersion']
	ax.plot(xM, xR, ls='-', marker='', color='b', label=epos.MR_label)
	ax.plot(xM, xR-dispersion, ls='--', marker='', color='b')
	ax.plot(xM, xR+dispersion, ls='--', marker='', color='b')

	# constant density used in N-body
	if rho is not None:
		rhoR= (xM*cgs.Mearth / (4./3.*np.pi*rho))**(1./3.) / cgs.Rearth
		ax.plot(xM, rhoR, ls='-', color='0.5', marker='', 
			label=r'$\rho={:.1f}\,g\,cm^{{-3}}$'.format(rho), zorder=-1)

	ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
	
	suffix= 'MC' if MC else 'input'
	suffix+='.log' if Log else ''	
	
	helpers.save(plt, '{}/massradius/{}'.format(epos.plotdir,suffix))		
