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
	
	#add decorators for zoom T/F?
	# color into epos?

	if epos.Parametric:
		parametric.oneD(epos)
		parametric.twoD(epos)
		parametric.panels(epos)

		# plot the input distribution converted from M->Msini or M->R
		if not epos.MonteCarlo and epos.Msini:
			parametric.oneD_y(epos, Convert=True) # no conversion for x (yet)
			parametric.twoD(epos, Convert=True)
			parametric.panels(epos, Convert=True)

	else:
		if 'R' in epos.pfm:
			model.panels_radius(epos, color=color)		
		if 'M' in epos.pfm:
			model.panels_mass(epos, color=color)
		
		model.period(epos, color=color)
			
		model.multiplicity(epos, color=color)
		if hasattr(epos, 'func'):
			if epos.MassRadius:
				model.panels_mass(epos, Population=True, color=color)
			else:
				model.panels_radius(epos, Population=True, Shade=False, color=color) #Zoom?

		if hasattr(epos, 'occurrence'):
			# model counts (debiases data)
			if 'planet' in epos.occurrence:
				model.panels_radius(epos, Occurrence=True, color=color)
				if 'R' in epos.pfm or epos.MassRadius:
					model.period(epos, Occurrence=True, color=color)
				if epos.Zoom:
					model.panels_radius(epos, Occurrence=True, color=color, Zoom=True)

			# planet counts (biases model)
			if 'model' in epos.occurrence:
				model.panels_radius(epos, Observation=True, color=color)
				if epos.Zoom:
					model.panels_radius(epos, Observation=True, color=color, Zoom=True)

		if 'tag' in epos.pfm:
			model.panels_radius(epos, Tag=True) # Zoom?
			model.period(epos, Tag=True)
			
		if 'inc' in epos.pfm:
			model.inclination(epos, Simple=True, color=color, imin=imin)
		if 'dP' in epos.pfm:
			model.periodratio(epos, color=color, Simple=True)
			if 'R' in epos.pfm:
				model.periodratio_size(epos, color=color)
			if 'inc' in epos.pfm:
				model.periodratio_inc(epos, color=color, imin=imin)

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



