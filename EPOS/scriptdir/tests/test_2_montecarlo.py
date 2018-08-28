#! /usr/bin/env ipython
''' 
Test if EPOS can do a Monte Carlo simulation

Plots should appear in the directory
png/test_2/input/
png/test_2/output/
'''

import EPOS

''' initialize the EPOS class '''
epos= EPOS.epos(name='test_2')

''' load the kepler dr25 exoplanets and survey efficiency '''
obs, survey= EPOS.kepler.dr25(Huber=True, Vetting=True, score=0.9)
epos.set_observation(**obs)
epos.set_survey(**survey)

''' define the parameteric distribution, here a power-law in radius and period '''
epos.set_parametric(EPOS.fitfunctions.powerlaw2D)
epos.fitpars.add('pps', 2.0, min=0)
epos.fitpars.add('P1',0.3, is2D=True)
epos.fitpars.add('P2',-0.2, dx=0.1, is2D=True)

''' define the simulated range (trim) and the range compared to observations (zoom) '''
epos.set_ranges(xtrim=[10,730],ytrim=[0.5,12.],xzoom=[50,200],yzoom=[0.7,3])

''' Run the Monte Carlo Simulation once '''
EPOS.run.once(epos)

''' plot parameteric distribution and simulated output '''
EPOS.plot.input.all(epos)
EPOS.plot.output.all(epos)