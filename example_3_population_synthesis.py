#! /usr/bin/env ipython
import numpy as np
import sys

import EPOS

'''
This example shows how to use EPOS with a planet populations synthesis model.
The model used here if from the review by Mordasini et al. 2018, Handbook of Exoplanets

You can replace EPOS.pfmodel.mordasini with a module that returns a dictionary with 
planet properties as numpy arrays. Planets with the same starID are assumed to be in 
a multi-planet systems
sma:	semi-major axis
mass:	planet mass
radius:	planet radius, optional if using a mass-radius relation
starID:	stellar identifier
inc:	planet mutual inclination
tag:	tag for simulation initial conditions such as Fe/H, optional
'''

''' initialize the EPOS class '''
epos= EPOS.classes.epos(name='example_3')

''' Load the Kepler planet list and detection efficiency '''
obs, survey= EPOS.kepler.dr25()
epos.set_observation(**obs)
epos.set_survey(**survey) #)

''' Define a parametric population for quick comparison '''
epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)
epos.pdfpars.add('pps', 2.0)
epos.pdfpars.add('P break', 10., is2D=True)
epos.pdfpars.add('P1',1.5, is2D=True)
epos.pdfpars.add('P2',0.0, is2D=True)
epos.pdfpars.add('R break',3.3, is2D=True) 
epos.pdfpars.add('R1',-0.5, is2D=True)
epos.pdfpars.add('R2',-6.0, is2D=True)

''' Load the planet formation model '''
pfm= EPOS.pfmodel.mordasini(name='syntheticpopmordasini1MsunJ31', dir='files', 
	cutoff= 1.,Verbose=True)
epos.set_population('Mordasini 2018', **pfm)
epos.fitpars.add('eta', 0.2, isnorm=True) # fraction of stars with simulated planets
epos.fitpars.add('f_iso', 0) # fraction of systems with isotropic orbits
epos.fitpars.add('f_inc', 1.0) # fudge factor for the inclinations
epos.fitpars.add('f_dP', 1.0) # fudge factor for the period ratios
epos.fitpars.add('f_cor', 0.5, fixed=True) # correlated noise

''' Uncomment to use a planet mass-radius conversion '''
#epos.set_massradius(EPOS.massradius.WRF15, 'WRF15+Rock', masslimits=[0.1,100])
#epos.set_massradius(EPOS.massradius.CK17, 'Chen & Kipping 2017', masslimits=[0.1,100])

''' Print some basic system properties'''
#print EPOS.analytics.pfmodel(epos)

''' define the simulated range (trim) and the range compared to observations (zoom)'''
epos.set_ranges(xtrim=[1,730],ytrim=[0.3,20.],xzoom=[1,400],yzoom=[0.5,14])

''' Calculate model occurrence rates
NOTE: to be replaced by EPOS.occurrence.all(epos) '''
epos.occurrence={}
EPOS.occurrence.planets(epos)
EPOS.occurrence.models(epos)

''' Run a forward model'''
EPOS.run.once(epos)

''' Run the MCMC to constrain fit parameters'''
#EPOS.run.mcmc(epos, nMC=100, nwalkers=20, nburn=20, threads= 5) # 5 mins
#EPOS.run.mcmc(epos, nMC=500, nwalkers=50, nburn=200, threads= 20) # 5 mins
#EPOS.run.mcmc(epos, nMC=2000, nwalkers=500, nburn=200, threads= 20) # more walkers?
#EPOS.run.mcmc(epos, nMC=10000, nwalkers=100, nburn=200, threads= 20) # looong

''' Generate Plots '''
EPOS.plot.input.all(epos, imin=1e-8)
#EPOS.plot.mcmc.all(epos)
EPOS.plot.output.all(epos)