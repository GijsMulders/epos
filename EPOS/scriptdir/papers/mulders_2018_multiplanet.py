#! /usr/bin/env ipython
''' 
This script is meant to reproduce some of the results/figures from
Mulders et al. 2018 (AJ accepted 14 May 2018), results from multi-planet mode:
- Table 2
- Figure 17, 18, 24
- (Figures similar to 6, 7, 8, 9, 11, 12, 13 can also be generated)

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

''' initialize the EPOS class '''
epos= EPOS.epos(name='mulders2018_multi')

''' load the kepler dr25 exoplanets and survey efficiency '''
obs, survey= EPOS.kepler.dr25(Huber=True, Vetting=True, score=0.9)
epos.set_observation(**obs)
epos.set_survey(**survey)

''' Define a double broken power-law as a planet population '''
epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)

''' Parameter initial guess and fitting ranges '''
epos.fitpars.add('pps',		0.4, 	min=0, 			isnorm=True)
epos.fitpars.add('P break',	10.,	min=2,	max=50,	is2D=True)
epos.fitpars.add('a_P',		1.5, 	min=0,			is2D=True)
epos.fitpars.add('b_P',		-1,		max=1,	dx=0.1,	is2D=True)
epos.fitpars.add('R break',	3.3,	fixed=True, 	is2D=True) 
epos.fitpars.add('a_R',		-0.5,	fixed=True, 	is2D=True)
epos.fitpars.add('b_R',		-6.,	fixed=True, 	is2D=True)

''' Set parameters for the spacing between planets'''
epos.set_multi(spacing='dimensionless')
epos.fitpars.add('npl', 10, fixed=True) # planets per system
epos.fitpars.add('log D', -0.3)
epos.fitpars.add('sigma', 0.2, min=0)

''' Some additional fit parameters'''
epos.fitpars.add('dR', 0.01, fixed=True)	# Dispersion in planet radii
epos.fitpars.add('inc', 2.0)				# mode of mutual inclinations
epos.fitpars.add('f_iso', 0.4)				# Fraction of isotropic systems
epos.fitpars.add('f_cor', 0.5, fixed=True)	# Correlated noise

''' define the simulated range (trim) and the range compared to observations (zoom)'''
epos.set_ranges(xtrim=[0,730],ytrim=[0.3,20.],xzoom=[2,400],yzoom=[0.5,6])

''' Run the Monte Carlo Simulation once '''
EPOS.run.once(epos)

''' run an MCMC chain on multiple cores or read in a previously saved run'''
EPOS.run.mcmc(epos, nMC=2000, nwalkers=100, nburn=200, threads= 30)

''' Change the titles of the corner plots '''
epos.fitpars.latexkeys={
	'pps':'$\eta_s$', 
	'P break':'$P_{\\rm break}$', 'a_P':'$a_P$', 'b_P':'$b_P$',
	'R break':'$R_{\\rm break}$', 'a_R':'$a_R$', 'b_R':'$b_R$',
	'npl':'$n_p$', 'log D':'log $D$', 'sigma':'$\sigma$',
	'inc':'$i$', 'f_iso':'$f_{\\rm iso}$', 'f_cor':'$f_{\\rm cor}$'
}

''' plot everything '''
EPOS.plot.survey.all(epos)
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.mcmc.all(epos)
EPOS.plot.occurrence.all(epos)