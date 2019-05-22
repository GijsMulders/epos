#! /usr/bin/env ipython
import numpy as np
import sys

import EPOS

'''
This example shows how to use EPOS with a planet formation/population synthesis model.

pfm is a dictionary with planet properties as numpy arrays. 
Planets with the same starID are assumed to be in a multi-planet systems
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
obs, survey= EPOS.kepler.dr25(Gaia=True, Vetting=True)
epos.set_observation(**obs)
epos.set_survey(**survey) #)

''' Define a parametric population for quick comparison '''
epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)
epos.pdfpars.add('pps', 2.6)
epos.pdfpars.add('P break', 11., is2D=True)
epos.pdfpars.add('P1',1.6, is2D=True)
epos.pdfpars.add('P2',0.1, is2D=True)
epos.pdfpars.add('R break',3.4, is2D=True) 
epos.pdfpars.add('R1',-0.1, is2D=True)
epos.pdfpars.add('R2',-6.7, is2D=True)

''' Load a planet formation model '''
# NOTE: load your own planet formation / population synthesis model here

# generate some random data for testing EPOS
n= 616
sma= 10.**np.random.uniform(-1.3,1,n)
mass= 10.** (3.*np.random.power(0.5, n))
radius= 10.** np.random.power(0.5, n)
inc= np.random.rayleigh(2, n)
starID= np.repeat(np.arange(n/8), 8)

pfm= {'sma':sma, 'mass':mass,'radius':radius, 'inc':inc, 'starID':starID}

# Load the planet formaton model into EPOS
epos.set_population('Planet Formation Model', **pfm)
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
epos.set_ranges(xtrim=[1,400],ytrim=[0.5,20.],xzoom=[1,400],yzoom=[0.5,20])

''' define rectangular occurrence rate bins'''
xbins= [[1,10], [20,300], [1, 20], [5,200]]
ybins= [[5,20], [4,15],   [0.7,2], [1.5,3.5]]
for xb, yb in zip(xbins, ybins):
	print 'area= {}'.format(np.log(yb[-1]/yb[0])*np.log(xb[-1]/xb[0]))
epos.set_bins(xbins=xbins,ybins=ybins)

''' Define polygonic occurrence rate bins'''
xmin, xmax= 1, 100
ymin, ymax= 0.7, 2.2
xb, yb= 3, 1.2
xyrock= [[xmin,ymin],[xmin,ymax],[xb,ymax],[xmax,yb],[xmax,ymin]]

xmin, xmax= 3, 200
ymin, ymax= 1.2, 4
xb, yb= 100, 2.2
xyMN= [[xmin,ymax],[xmax,ymax],[xmax,ymin], [xb,ymin], [xmin,yb]]

xyHJ= [[xbins[0][0],ybins[0][1]],[xbins[0][1],ybins[0][1]],
  [xbins[0][1],ybins[0][0]],[xbins[0][0],ybins[0][0]]]
xyTG= [[xbins[1][0],ybins[1][1]],[xbins[1][1],ybins[1][1]],
  [xbins[1][1],ybins[1][0]],[xbins[1][0],ybins[1][0]]]

epos.set_bins_poly([xyrock, xyMN, xyHJ, xyTG], 
  labels=['Hot\nEarths','Mini-\nNeptunes','Hot\nJupiters','Warm\nGiants'])

''' Run a forward model'''
EPOS.run.once(epos)

''' Calculate model occurrence rates '''
EPOS.occurrence.all(epos)

''' Save occurence rate to disk '''
EPOS.save.survey(epos)
EPOS.save.model(epos)
EPOS.save.synthetic_survey(epos)
EPOS.save.occurrence(epos)

''' Generate Plots '''
EPOS.plot.survey.all(epos)
EPOS.plot.input.all(epos, imin=1e-6, color='C8')
EPOS.plot.occurrence.all(epos, color='C8', alpha_fac=50.)
EPOS.plot.output.all(epos, color='C8')