#! /usr/bin/env ipython
''' 
Example script fitting a power-law in period and M sin i to radial velocity data
'''
import EPOS

''' initialize the EPOS class '''
epos= EPOS.epos(name='example_5', RV=True, MC=False, Msini=True)

''' load the exoplanets and completeness from Mayor+ 2011, Fernandes+ 2018'''
obs, survey= EPOS.rv.Mayor2011()
epos.set_observation(**obs)
epos.set_survey(**survey)

''' Define a double broken power-law as a planet population '''
epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)

''' Parameter initial guess and fitting ranges
Note: 
	- brokenpowerlaw2D uses 6 parameters, indicated with the is2D keyword
	- 'pps' is a normalization factor for the planet occurrence rate (planets-per-star)
	- parameters a_M and b_M are not fitted
	- dx is the range in walker initial positions for parameters that change sign (+/-)
'''
epos.fitpars.add('pps',		1.0, 	min=1e-3)
epos.fitpars.add('P break',	1e3,	min=100,max=7e3,is2D=True)
epos.fitpars.add('a_P',		1.0, 	min=0,	max=3, 	is2D=True)
epos.fitpars.add('b_P',		-0.5,	min=-3, max=0,	is2D=True)
epos.fitpars.add('M break',	10.0,	fixed=True, 	is2D=True) 
epos.fitpars.add('a_M',		0.0,	fixed=True, 	is2D=True)
epos.fitpars.add('b_M',		-0.5,	dx=0.1,		 	is2D=True)

''' define the simulated range (trim) and the range compared to observations (zoom)'''
epos.set_ranges(xtrim=[1,1e5],ytrim=[1,1e5],xzoom=[10,1e4],yzoom=[50,1e4], Occ=True)

''' Run the Monte Carlo Simulation once '''
EPOS.run.once(epos)

''' run an MCMC chain on multiple cores or read in a previously saved run'''
EPOS.run.mcmc(epos, nMC=1000, nwalkers=100, nburn=200, threads=20, Saved=True)

''' define bins where occurrence is calculated'''
#epos.set_bins(xbins=[[2,400],[0.9*365,2.2*365]], ybins=[[1,6],[0.7,1.5]]) # eta_zoom, eta_earth

''' Calculate the occurrence rates '''
EPOS.occurrence.all(epos)

''' Adjust plot parameters'''
epos.plotpars['occrange']= [2e-4,2.]

''' plot everything '''
EPOS.plot.survey.all(epos)
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.mcmc.all(epos)
EPOS.plot.occurrence.all(epos)