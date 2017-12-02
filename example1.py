#! /usr/bin/env ipython
import EPOS

''' This example demonstrates the core functionality of EPOS '''

# initialize the EPOS class
epos= EPOS.epos(name='example1')

# load the observed exoplanets and survey efficiency
obs, survey= EPOS.kepler.dr25()
epos.set_observation(**obs)
epos.set_survey(**survey)

# define the simulated range (trim) and the range compared to observations (zoom)
epos.set_ranges(xtrim=[10,730],ytrim=[0.5,12.],xzoom=[50,200],yzoom=[0.7,3])

# define the parameteric distribution
epos.set_parametric(EPOS.fitfunctions.powerlaw2D)
epos.fitpars.add('pps', 2.0, min=0)
epos.fitpars.add('P1',0.3, is2D=True)
epos.fitpars.add('P2',-0.2, dx=0.1, is2D=True)

# uncomment for a broken power-law
#epos.set_ranges(xtrim=[0,300],ytrim=[0.5,12.],xzoom=[2,200],yzoom=[0.7,8])
# epos.set_parametric(EPOS.fitfunctions.brokenpowerlaw2D)
# epos.fitpars.add('pps', 2.0, min=0)
# epos.fitpars.add('P break', 10.5, min=0, is2D=True)
# epos.fitpars.add('P1',2.0, is2D=True)
# epos.fitpars.add('P2',0.35, dx=0.1, is2D=True)
# epos.fitpars.add('R break',2.5, min=0, is2D=True) 
# epos.fitpars.add('R1',-0.02, dx=0.1, is2D=True)
# epos.fitpars.add('R2',-3., is2D=True)

# one MC
EPOS.run.once(epos) # == prep

# occurrence rate
epos.set_bins(xbins=[[300,700]], ybins=[[0.5,1.4]]) # habitable zone-ish
EPOS.occurrence.all(epos)


# plot stuff
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)
EPOS.plot.occurrence.all(epos)