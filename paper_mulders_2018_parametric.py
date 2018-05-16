#! /usr/bin/env ipython
''' 
This script is meant to reproduce some of the results/figures from
Mulders et al. 2018 (AJ accepted 14 May 2018), results from parametric mode:
- Table 1
- Figure 14, 22, 23
- (Figures similar to 2, 4, 5, 10, 15, 16 can also be generated)

Text output will appear in the terminal.
Plots will be generated in the png/ subfolder 

Note that results may vary depending on the random seed, version of the code used, 
and external dependencies. 
Future versions of the code may employ different default parameters.
This script is compatible with version 1.0 of EPOS and was added to the github repository May 15 2018

If you have trouble running this script or EPOS, please consult the online documentation, 
try running the test and example scripts, or contact the first author
'''
import EPOS

''' initialize the EPOS class, with the same random seed as used in the paper '''
epos= EPOS.epos(name='mulders2018_parametric')

''' load the kepler dr25 exoplanets and survey efficiency '''
obs, survey= EPOS.kepler.dr25(Huber=True, Vetting=True, score=0.9)
epos.set_observation(**obs)
epos.set_survey(**survey)

''' Define a double broken power-law as a planet population '''
epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)

''' Parameter initial guess and fitting ranges
'''
epos.fitpars.add('pps',		2.0, 	min=0, 			isnorm=True)
epos.fitpars.add('P break',	10.,	min=2,	max=50,	is2D=True)
epos.fitpars.add('a_P',		1.5, 	min=0,			is2D=True)
epos.fitpars.add('b_P',		0.0,	dx=0.1,			is2D=True)
epos.fitpars.add('R break',	3.0,	min=1.0,max=5, 	is2D=True) 
epos.fitpars.add('a_R',		0.0,	dx=0.1, 		is2D=True)
epos.fitpars.add('b_R',		-4.,	min=-10, max=0, is2D=True)

''' define the simulated range (trim) and the range compared to observations (zoom)'''
epos.set_ranges(xtrim=[0,730],ytrim=[0.3,20.],xzoom=[2,400],yzoom=[0.5,6], Occ=True)

''' Run the Monte Carlo Simulation once '''
EPOS.run.once(epos)

''' run an MCMC chain on multiple cores or read in a previously saved run'''
EPOS.run.mcmc(epos, nMC=5000, nwalkers=200, nburn=1000, threads= 20, Saved=True)

''' define bins where occurrence rates are calculated'''
epos.set_bins(xbins=[[2,400],[0.9*365,2.2*365]], ybins=[[0.5,6.0],[0.7,1.5]]) # eta_earth

''' Calculate the occurrence rates '''
EPOS.occurrence.all(epos)

''' Adjust axes'''
epos.plotpars['textsize']= 10 # smaller text in plot
epos.xtrim[1]= 1000 # adjust plot range to include habitable zone

''' Better labels for corner.py '''
epos.fitpars.latexkeys={
	'pps':'$\eta$', 
	'P break':'$P_{\\rm break}$', 'a_P':'$a_P$', 'b_P':'$b_P$',
	'R break':'$R_{\\rm break}$', 'a_R':'$a_R$', 'b_R':'$b_R$',
}

''' plot everything '''
EPOS.plot.survey.all(epos)
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.mcmc.all(epos)
EPOS.plot.occurrence.all(epos)