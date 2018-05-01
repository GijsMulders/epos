#! /usr/bin/env ipython
''' 
Test if EPOS can calculate occurrence rates

Plots should appear in the directory
png/test_5/occurrence/
'''

import EPOS

''' initialize the EPOS class '''
epos= EPOS.epos(name='test_5')

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
epos.set_ranges(xtrim=[10,730],ytrim=[0.5,4.],xzoom=[20,300],yzoom=[0.7,3], Occ=True)

''' define bins where occurrence is calculated'''
epos.set_bins(xbins=[[20,300],[0.9*365,2.2*365]], ybins=[[0.7,3],[0.7,1.5]]) # eta_zoom, eta_earth

''' Run the Monte Carlo Simulation once '''
EPOS.run.once(epos)

''' Calculate the occurrence rates '''
EPOS.occurrence.all(epos)

''' Adjust axes'''
epos.plotpars['textsize']= 12
epos.xtrim[1]= 1000

''' plot occurrence rates '''
EPOS.plot.occurrence.all(epos)