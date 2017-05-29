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
epos.set_survey(eff_P, eff_R, eff_2D)
epos.set_ranges(xtrim=[0,300],ytrim=[0.5,12.],xzoom=[2,200],yzoom=[0.7,8])

bp1d= EPOS.fitfunctions.brokenpowerlaw1D
epos.set_parametric(xfunc=bp1d, yfunc=bp1d, xparams=[10.5, 1.5, 0.35], yparams=[2.5, -0.2, -3.0], normalization=0.005)

EPOS.run.once(epos)

EPOS.plot.input(epos)
EPOS.plot.output(epos)
