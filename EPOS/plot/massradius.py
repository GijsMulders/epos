import numpy as np
import matplotlib.pyplot as plt
import helpers
	
def massradius(epos, MC=False, Log=False):
	assert epos.RadiusMassConversion
	# plot R(M)
	f, ax = plt.subplots()
	ax.set_title('mass-radius relation + dispersion')
	ax.set_xlabel(r'Planet Mass [M$_\bigoplus$]')
	ax.set_ylabel(r'Planet Radius [R$_\bigoplus$]')

	if Log:
		ax.set_xlim(0.5, 1000.) # 20.
		ax.set_ylim(0.5, 20.) # 7.	
		ax.set_xscale('log')
		ax.set_yscale('log')
	else:
		ax.set_xlim(0, 20.) # 20.
		ax.set_ylim(0, 7.) # 7.
	
	# MC data
	if MC:
		ss=epos.synthetic_survey
		ax.plot(ss['M'], ss['R'], ls='', marker='.', mew=0, ms=5.0, color='k', zorder=1)
		
	xM= np.logspace(*np.log10(ax.get_xlim())) if Log else np.linspace(*ax.get_xlim()) #np.max(tr['M']))
	xR, dispersion= epos.RM(xM)
	#ax.plot(xM, epos.MR['fRock'](xM), ls='-', marker='', color='r', label='Rocky constraint', zorder=0)
	#xGas= epos.MR['fGas'](xM)
	#dispersion= epos.MR['dispersion']
	ax.plot(xM, xR, ls='-', marker='', color='b', label=epos.RM_label)
	ax.plot(xM, xR-dispersion, ls='--', marker='', color='b')
	ax.plot(xM, xR+dispersion, ls='--', marker='', color='b')

	ax.legend(loc='lower right', shadow=False, prop={'size':14}, numpoints=1)
	
	suffix= 'MC' if MC else 'input'
	suffix+='.log' if Log else ''	
	
	helpers.save(plt, '{}/massradius/{}'.format(epos.plotdir,suffix))		