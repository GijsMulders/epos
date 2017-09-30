#! /usr/bin/env ipython
import EPOS

''' This example demonstrates the MCMC fitting of EPOS '''

EPOS.readme.all()

# initialize the EPOS class
epos= EPOS.classes.epos(name='example2')

# load the observed exoplanets and survey efficiency
EPOS.kepler.readme()
obs, survey= EPOS.kepler.dr25()
epos.set_observation(**obs)
epos.set_survey(**survey)

# define the simulated range (trim) and the range compared to observations (zoom)
epos.set_ranges(xtrim=[0,300],ytrim=[0.5,12.],xzoom=[2,200],yzoom=[0.7,8])

# use a double broken power-law
epos.set_parametric(func=EPOS.fitfunctions.brokenpowerlaw2D, 
					pname=['norm','P break','P1','P2','R break','R1','R2'],
					p0=[0.01,10., 2., 0.2,2.5, 0.1, -3.0]) # starting guess

# one MC
EPOS.run.once(epos) # == prep

# run the MCM chain (single core)
#EPOS.run.mcmc(epos, nMC=100, nwalkers=20, nburn=50) # 2 minutes

# run this instead for a longer MCMC chain (8 threads)
#EPOS.run.mcmc(epos, nMC=1000, nwalkers=100, nburn=200, threads=8)

# plot the MCMC chain
EPOS.plot.mcmc.all(epos)

# plot stuff
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
