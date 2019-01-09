#! /usr/bin/env ipython
''' 
This script is meant to reproduce some of the results/figures from
Fernandes et al. 2019 (AJ submitted), in particular the 3 and 4 parameter 
solutions in Table 1 and 2 and Figure 4, (9), and 10

Text output will appear in the terminal.

Symmetric:
./paper_fernandes_2019_broken_powerlaw.py sym

Best-fit values
  pps= 0.476 +0.0798 -0.0704
  P break= 1.54e+03 +948 -381
  p1= 0.652 +0.197 -0.17
  m1= -0.442 +0.0526 -0.0544

  posterior per bin
  eta= 26.6% +6.1% -5.4%

Asymmetric:
./paper_fernandes_2019_broken_powerlaw.py

Best-fit values
  pps= 0.462 +0.0804 -0.0701
  P break= 2.04e+03 +1.16e+03 -1.13e+03
  p1= 0.641 +0.272 -0.146
  p2= -1.29 +0.96 -1.17
  m1= -0.439 +0.0538 -0.0536

  posterior per bin
  eta= 27.4% +7.0% -5.4%

Plots will be generated in the png/ subfolder
png/fernandes2019_rv/occurrence/posterior_x.png 
png/fernandes2019_rv/occurrence/posterior_y.png 
png/fernandes2019_rv_sym/mcmc/triangle.png 

Note that results may vary depending on the random seed, version of the code used, 
and external dependencies. 
Future versions of the code may employ different default parameters.
This script is compatible with version 1.1 of EPOS

If you have trouble running this script or EPOS, please consult the online documentation, 
try running the test and example scripts, or contact the first author

'''
import sys
import EPOS
from EPOS import cgs

Symmetric= 'sym' in sys.argv
Single= 'single' in sys.argv

if Single: suffix='_single'
elif Symmetric: suffix='_sym'
else: suffix=''

''' initialize the EPOS class '''
#Msini = True converts the mass distirbution to an msini distribution
epos= EPOS.epos(name='fernandes2019_rv{}'.format('_sym' if Symmetric else ''), 
	RV = True, MC = False, Msini = True)

''' load the exoplanets and completeness from Mayor+ 2011'''
obs, survey= EPOS.rv.Mayor2011()
epos.set_observation(**obs)
epos.set_survey(**survey)

''' Define a double broken power-law as a planet population '''
if Single:
  epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D_yonly)
elif Symmetric:
	epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D_symmetric)
else:
	epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)

''' Parameter initial guess and fitting ranges
Note: 
	- brokenpowerlaw2D uses 6 parameters, indicated with the is2D keyword
	- 'pps' is a normalization factor for the planet occurrence rate (planets-per-star)
	- parameters a_M and b_M are not fitted
	- dx is the range in walker initial positions for parameters that change sign (+/-)
'''
epos.fitpars.add('pps',		1.0, 	min=1e-3)
if Single:
  epos.fitpars.add('p1',    1.0,  min=0,  max=3,  is2D=True)
else:
  epos.fitpars.add('P break', 1e3,  min=100,max=7e3,is2D=True)
  epos.fitpars.add('p1',    1.0,  min=0,  max=3,  is2D=True)
  if not Symmetric:
    epos.fitpars.add('p2',		-0.5,	min=-3, max=0,	is2D=True)
epos.fitpars.add('M break',	10.0,	fixed=True, 	is2D=True) 
epos.fitpars.add('a_M',		0.0,	fixed=True, 	is2D=True)
epos.fitpars.add('m1',		-0.5,	fixed=False, dx=0.1,	is2D=True)

''' define the simulated range (trim) and the range compared to observations (zoom)'''
epos.set_ranges(xtrim=[1,1e5],ytrim=[10,1e5],xzoom=[10,1e4],yzoom=[30,6000], 
	Occ=True, UnitTicks=True)

''' Run the simulation once '''
EPOS.run.once(epos)

''' run an MCMC chain on multiple cores or read in a previously saved run'''
EPOS.run.mcmc(epos, nMC=1000, nwalkers=50, nburn=200, threads=20, Saved=True)

''' define bins where occurrence is calculated '''
Pbin= [365.24 * au**1.5 for au in [0.1,100.]]
Mbin= [Mj * cgs.Mjup/cgs.Mearth for Mj in [0.1,13.]] 
epos.set_bins(xbins=[Pbin], ybins=[Mbin])

''' Calculate the occurrence rates '''
EPOS.occurrence.all(epos)

''' Adjust plot parameters'''
epos.plotpars['occrange']= [2e-4,2.]
#epos.xtrim[1]= 4e5

''' plot everything '''
EPOS.plot.survey.all(epos)
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.mcmc.all(epos)
EPOS.plot.occurrence.all(epos)