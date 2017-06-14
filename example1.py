#! /usr/bin/env ipython
import EPOS

''' This example demonstrates the core functionality of EPOS '''

EPOS.readme.all()

# initialize the EPOS class
epos= EPOS.classes.epos(name='example1')

# load the observed exoplanets and survey efficiency
EPOS.kepler.readme()
obs_P, obs_R, star_ID= EPOS.kepler.obs_Q16()
eff_P, eff_R, eff_2D= EPOS.kepler.eff_Q16()

epos.set_observation(obs_P, obs_R, star_ID, nstars=1.6862e5)
epos.set_survey(eff_P, eff_R, eff_2D) #) 

# define the simulated range (trim) and the range compared to observations (zoom)
epos.set_ranges(xtrim=[10,300],ytrim=[0.5,12.],xzoom=[50,200],yzoom=[0.7,3])

# define the type of planet population to fit
epos.set_parametric(func=EPOS.fitfunctions.powerlaw2D, 
					pname=['norm','P1','P2'],
					p0=[0.002, 0.3,-0.2]) # 

# uncomment for a broken power-law
#epos.set_ranges(xtrim=[0,300],ytrim=[0.5,12.],xzoom=[2,200],yzoom=[0.7,8])
#epos.set_parametric(func=EPOS.fitfunctions.brokenpowerlaw2D, 
#					pname=['norm','P break','P1','P2','R break','R1','R2'],
#					p0=[0.005,10.5, 2.0, 0.35, 2.5, -0.2, -3.0]) 

# one MC
EPOS.run.once(epos) # == prep

# plot stuff
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
