#! /usr/bin/env ipython
''' 
Example script fitting a double broken power-law in period and radius to Kepler data
'''
import EPOS

''' initialize the EPOS class '''
epos= EPOS.epos(name='example_1')

''' load the kepler dr25 exoplanets and survey efficiency '''
obs, survey= EPOS.kepler.dr25(Huber=True, Vetting=True, score=0.9)
epos.set_observation(**obs)
epos.set_survey(**survey)

''' Define a double broken power-law as a planet population '''
epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)

''' Parameter initial guess and fitting ranges
Note: 
	- brokenpowerlaw2D uses 6 parameters, indicated with the is2D keyword
	- 'pps' is a normalization factor for the planet occurrence rate (planets-per-star)
	- parameter b_R is not fitted
	- dx is the range in walker initial positions for parameters that change sign (+/-)
'''
epos.fitpars.add('pps',		2.0, 	min=0, 			isnorm=True)
epos.fitpars.add('P break',	10.,	min=2,	max=50,	is2D=True)
epos.fitpars.add('a_P',		1.5, 	min=0,			is2D=True)
epos.fitpars.add('b_P',		0.0,	dx=0.1,			is2D=True)
epos.fitpars.add('R break',	3.0,	min=1.0,max=5, 	is2D=True) 
epos.fitpars.add('a_R',		0.0,	dx=0.1, 		is2D=True)
epos.fitpars.add('b_R',		-4.,	fixed=True, 	is2D=True)

''' define the simulated range (trim) and the range compared to observations (zoom)'''
epos.set_ranges(xtrim=[0,730],ytrim=[0.3,20.],xzoom=[2,400],yzoom=[1,6], Occ=True)

''' Run the Monte Carlo Simulation once '''
EPOS.run.once(epos)

''' run an MCMC chain on multiple cores or read in a previously saved run'''
EPOS.run.mcmc(epos, nMC=1000, nwalkers=100, nburn=200, threads=20, Saved=True)

''' define bins where occurrence is calculated'''
epos.set_bins(xbins=[[2,400],[0.9*365,2.2*365]], ybins=[[1,6],[0.7,1.5]]) # eta_zoom, eta_earth

''' Calculate the occurrence rates '''
EPOS.occurrence.all(epos)

''' Adjust axes'''
epos.plotpars['textsize']= 12
epos.xtrim[1]= 1000

''' plot everything '''
EPOS.plot.survey.all(epos)
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.mcmc.all(epos)
EPOS.plot.occurrence.all(epos)