import matplotlib, warnings
#warnings.filterwarnings("ignore", category=UserWarning, module='matplotlib')
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

import parametric, multi, massradius, helpers, model

def all(epos, color=None, imin=1e-2):
	print '\nPlotting input...'
	
	if epos.Parametric:
		parametric.oneD(epos)
		parametric.twoD(epos)
		parametric.panels(epos)
	else:
		if 'R' in epos.pfm:
			model.panels_radius(epos, color=color)		
		if 'M' in epos.pfm:
			model.panels_mass(epos, color=color)
		
		model.period(epos)
			
		model.multiplicity(epos, color=color)
		if hasattr(epos, 'func'):
			if epos.MassRadius:
				model.panels_mass(epos, Population=True, color=color)
			else:
				model.panels_radius(epos, Population=True, color=color)

		if hasattr(epos, 'occurrence'):
			if 'planet' in epos.occurrence:	
				model.panels_radius(epos, Occurrence=True, color=color)
			if 'model' in epos.occurrence:
				model.panels_radius(epos, Observation=True, color=color)

		if 'tag' in epos.pfm:
			model.panels_radius(epos, Tag=True)
			model.period(epos, Tag=True)
			
		if 'inc' in epos.pfm:
			model.inclination(epos, color=color, imin=imin)
		if 'dP' in epos.pfm:
			model.periodratio(epos, color=color)
			if 'R' in epos.pfm:
				model.periodratio_size(epos, color=color)
			
# 		input(epos, PlotBox=False)
# 		input_diag(epos, PlotBox=False)
# 		if 'all_Pratio' in epos.groups[0]:
# 			Pratio(epos)
# 		if 'all_inc' in epos.groups[0]:
# 			inclination(epos)

	if epos.Multi:
		if epos.Parametric:
			multi.periodradius(epos, MC=False)
			multi.periodradius(epos, MC=False, Nth=True)
			multi.multiplicity(epos, MC=False)
			multi.periodratio(epos, MC=False)
			multi.periodratio_cdf(epos, MC=False)
			multi.periodinner(epos, MC=False)
			multi.periodinner_cdf(epos, MC=False)
		else:
			pass
		
	if epos.MassRadius: 
		massradius.massradius(epos, MC=False)
		massradius.massradius(epos, MC=False, Log=True)



