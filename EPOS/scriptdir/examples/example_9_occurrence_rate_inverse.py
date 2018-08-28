#! /usr/bin/env ipython
''' 
Estimate plane occurrence rates using the inverse detection efficiency method

Plots should appear in the directory
png/example_9/occurrence/

Binned occurrence rates are saved in 
json/example_9/

This example calculates the occurrence rates for different classes of planets, 
defined by a size and period-range
'''

import EPOS

''' initialize the EPOS class '''
epos= EPOS.epos(name='example_9')

''' load the kepler dr25 exoplanets and survey efficiency '''
obs, survey= EPOS.kepler.dr25(Huber=True, Vetting=True, score=0.9)
epos.set_observation(**obs)
epos.set_survey(**survey)

''' Set plotting ranges (zoom ranges are not used) '''
epos.set_ranges(xtrim=[0.5,830],ytrim=[0.3,20.],xzoom=[20,300],yzoom=[0.7,3])

''' define a set of period (x) and radius (y) bins to calculate occurrence rates'''
x_HJ,y_HJ= [1,10],[6,20] # Hot-Jupiters
x_SEMN, y_SEMN= [2,150],[1.0,4.0] # super-earths and mini-Neptunes

''' Pass the list of bins to EPOS'''
#epos.set_bins(xbins=[x_HJ,x_SEMN], 
#	ybins=[y_HJ,y_SEMN])

import numpy as np
epos.set_bins(xgrid=np.geomspace(10,640,7), 
   	ygrid=np.geomspace(0.67,17,9), Grid=True)

''' Calculate the occurrence rates '''
EPOS.occurrence.all(epos)

''' Save the occurrence rates per bin in '''
EPOS.save.occurrence(epos)

''' Adjust axes to plot habitable zone '''
epos.plotpars['textsize']= 8
epos.xtrim[1]= 1000

''' plot occurrence rates '''
EPOS.plot.occurrence.all(epos)