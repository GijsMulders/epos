import numpy as np
import matplotlib.pyplot as plt
import periodradius, massradius, multi, helpers

clrs= ['r','g','b','m'] # in epos.prep
fmt_symbol= {'ls':'', 'marker':'o', 'mew':2, 'ms':8,'alpha':0.6}

''' output '''
def all(epos):

	if hasattr(epos, 'synthetic_survey'):
		print '\nPlotting output...'
	
		if epos.MonteCarlo:
			periodradius.periodradius(epos, Parametric=epos.Parametric, SNR=False)
			periodradius.periodradius(epos, Parametric=epos.Parametric, SNR=True)
		else:
			# plot period-radius for non-MC
			pass
		
		periodradius.cdf(epos)
		periodradius.panels(epos)
	
		if epos.Multi:
			multi.periodradius(epos)
			multi.periodradius(epos, Nth=True)
	
			multi.multiplicity(epos, MC=True)
			multi.multiplicity(epos, MC=True, Planets=True)
			multi.multiplicity_cdf(epos, MC=True)
			multi.periodratio(epos, MC=True)
			if epos.Parametric and epos.spacing is not None:
				multi.periodratio(epos, MC=True, Input=True)
			multi.periodratio_cdf(epos, MC=True)
			multi.periodinner(epos, MC=True)
			if epos.Parametric and epos.spacing is not None:
				multi.periodinner(epos, MC=True, Input=True)
			multi.periodinner_cdf(epos, MC=True)
			# pdf per subgroup
			#periodradius.pdf(epos)
			#periodradius.pdf_3d(epos)

		if epos.MassRadius:
			massradius.massradius(epos, MC=True)
			massradius.massradius(epos, MC=True, Log=True)

		else:
			if epos.Parametric and epos.Multi and not epos.RV:
				multi.periodratio(epos, MC=True, N=True)
				if epos.spacing is not None:
					multi.periodratio(epos, MC=True, N=True, Input=True)
				multi.periodinner(epos, MC=True, N=True)

	else:
		print '\nNo output to plot, did you run EPOS.run.once()? \n'
