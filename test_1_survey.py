#! /usr/bin/env ipython
''' 
Test if EPOS can read exoplanet data in files/

Plots should appear in the directory
png/test_1/survey/
'''

import EPOS

''' initialize the EPOS class '''
epos= EPOS.epos(name='test_1')

''' load the kepler dr25 exoplanets and survey efficiency '''
obs, survey= EPOS.kepler.dr25(Huber=True, Vetting=True, score=0.9)
epos.set_observation(**obs)
epos.set_survey(**survey)

''' plot survey '''
EPOS.plot.survey.all(epos)
