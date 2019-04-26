import numpy as np
import matplotlib.pyplot as plt
import periodradius, massradius, multi, helpers

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

''' output '''
def all(epos, color='C1'):

	if hasattr(epos, 'synthetic_survey'):
		print '\nPlotting output...'
	
		if epos.MonteCarlo:
			periodradius.periodradius(epos, SNR=False)
			periodradius.periodradius(epos, SNR=True)
			if not epos.Parametric:
				periodradius.periodradius(epos, Model=True, SNR=False, color=color)
		else:
			# plot period-radius for non-MC
			pass
		
		periodradius.cdf(epos, color=color)
		periodradius.panels(epos, color=color)
	
		if epos.Multi:
			multi.periodradius(epos)
			multi.periodradius(epos, Nth=True)
	
			multi.multiplicity(epos, MC=True, color=color)
			multi.multiplicity(epos, MC=True, Planets=True, color=color)
			multi.multiplicity_cdf(epos, MC=True)
			multi.periodratio(epos, MC=True, color=color)
			if epos.Parametric and epos.spacing is not None:
				multi.periodratio(epos, MC=True, Input=True, color=color)
			multi.periodratio_cdf(epos, MC=True, color=color)
			multi.periodinner(epos, MC=True, color=color)
			if epos.Parametric and epos.spacing is not None:
				multi.periodinner(epos, MC=True, Input=True, color=color)
			multi.periodinner_cdf(epos, MC=True, color=color)
			# pdf per subgroup
			#periodradius.pdf(epos)
			#periodradius.pdf_3d(epos)

		if epos.MassRadius:
			massradius.massradius(epos, MC=True, color=color)
			massradius.massradius(epos, MC=True, Log=True, color=color)

		else:
			if epos.Parametric and epos.Multi and not epos.RV:
				multi.periodratio(epos, MC=True, N=True)
				if epos.spacing is not None:
					multi.periodratio(epos, MC=True, N=True, Input=True)
				multi.periodinner(epos, MC=True, N=True)

	else:
		print '\nNo output to plot, did you run EPOS.run.once()? \n'
