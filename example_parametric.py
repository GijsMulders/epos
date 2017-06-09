#! /usr/bin/env ipython

import EPOS

EPOS.readme.all()

# initialize the EPOS class
epos= EPOS.classes.epos(name='para') #, xzoom=[7,500], yzoom=[1.5,4])


# load the observed exoplanets and survey efficiency
EPOS.kepler.readme()
obs_P, obs_R, star_ID= EPOS.kepler.obs_Q16()
eff_P, eff_R, eff_2D= EPOS.kepler.eff_Q16()

epos.set_observation(obs_P, obs_R, star_ID, nstars=1.6862e5)
epos.set_survey(eff_P, eff_R, eff_2D) #) 
epos.set_ranges(xtrim=[0,300],ytrim=[0.5,12.],xzoom=[2,200],yzoom=[0.7,8])

epos.set_parametric(func=EPOS.fitfunctions.brokenpowerlaw2D, 
					pname=['norm','P break','P1','P2','R break','R1','R2'],
					p0=[0.005,10.5, 2.0, 0.35, 2.5, -0.2, -3.0]) # 'fit'
					#p0=[0.01,10., 2., 0.1,2.5, 0.1, -3.0]) # starting guess, bad when 0

# one MC
EPOS.run.once(epos) # == prep

# uncomment one of these for running an MCMC chain
#EPOS.run.mcmc(epos, nMC=100, nwalkers=20, nburn=50) # 2 minutes
#EPOS.run.mcmc(epos, nMC=500, nwalkers=100, nburn=200) # 1 hour...

# uncomment this to plot the MCMC chain
#EPOS.plot.mcmc.all(epos)

# plot stuff
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
