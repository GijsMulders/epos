#! /usr/bin/env ipython
''' 
Example script fitting a population of multi-planet systems to Kepler data
'''
import EPOS

''' initialize the EPOS class '''
epos= EPOS.epos(name='example_2')

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
epos.fitpars.add('pps',		0.4, 	min=0)
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
epos.set_ranges(xtrim=[0,730],ytrim=[0.3,20.],xzoom=[2,400],yzoom=[1,6])

''' Run the Monte Carlo Simulation once '''
EPOS.run.once(epos)

''' run an MCMC chain on multiple cores or read in a previously saved run'''
EPOS.run.mcmc(epos, nMC=500, nwalkers=50, nburn=200, threads= 30, Saved=True) # ~20 mins
#EPOS.run.mcmc(epos, nMC=1000, nwalkers=100, nburn=200, threads=20, Saved=True) # ~5 hrs

''' Change the titles of the corner plots '''
epos.fitpars.latexkeys={
	'pps':'$\eta_s$', 
	'P break':'$P_{\\rm break}$', 'P1':'$a_P$', 'P2':'$b_P$',
	'P break':'$P_{\\rm break}$', 'R1':'$a_R$', 'R2':'$b_R$',
	'npl':'$n_p$', 'log D':'log $D$', 'sigma':'$\sigma$',
	'inc':'$i$', 'f_iso':'$f_{\\rm iso}$', 'f_cor':'$f_{\\rm cor}$'
}

''' plot everything '''
EPOS.plot.survey.all(epos)
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.mcmc.all(epos)
EPOS.plot.occurrence.all(epos)