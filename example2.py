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
epos.set_ranges(xtrim=[0,730],ytrim=[0.5,12.],xzoom=[2,200],yzoom=[0.7,8])

# use a double broken power-law
epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)
epos.fitpars.add('pps', 2.0, min=0)
epos.fitpars.add('P break', 10., min=0, is2D=True)
epos.fitpars.add('P1',2.0, is2D=True)
epos.fitpars.add('P2',0.0, dx=0.1, is2D=True)
epos.fitpars.add('R break',3, min=0, is2D=True) 
epos.fitpars.add('R1',0, dx=0.1, is2D=True)
epos.fitpars.add('R2',-3., is2D=True)

# one MC
EPOS.run.once(epos) # == prep

# run the MCM chain (single core)
#EPOS.run.mcmc(epos, nMC=100, nwalkers=20, nburn=50) # 2 minutes

# run this instead for a longer MCMC chain (8 threads)
#EPOS.run.mcmc(epos, nMC=1000, nwalkers=100, nburn=200, threads=8)

# occurrence rate?
#epos.set_bins(xbins=[[300,700]], ybins=[[0.5,1.4]]) # habitable zone-ish
#EPOS.occurrence.all(epos)

# plot stuff
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.mcmc.all(epos)
EPOS.plot.occurrence.all(epos)