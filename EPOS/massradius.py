''' 
Helper function for setting mass-radius distributions
Mass, Radius in earth masses
'''
import numpy as np

def WRF15(mass):
	Mp=np.asarray(mass)
	
	Rgas= (Mp/2.7)**(1/1.3)
	# quick fit to Zeng (see PebbleAccretion/density.py)
	Rrock= Mp**(1./3.7) 

	isrock= Rrock > Rgas
	Rgas[Rgas>12]=12.

	mean=np.where(isrock,Rrock,Rgas)
	#dispersion= np.where(isrock,0.2*Rrock, 1.9/1.3/2.7) # WRF15 estimate
	dispersion=0.2*mean # ~0.5 Rearth at 2.5 Rearth seems reasonable
	
	return mean, dispersion

def CK17(mass):
	''' Chen & Kipping 2017, Forecasting'''
	Mp=np.asarray(mass)
	
	# transitions between rocky, neptune, gas giant
	Mrockmax= 2.0
	Mnepmax= 132 # 0.414 Mjup
	
	#irock= Mp<=Mrockmax
	inep= (Mrockmax<Mp) & (Mp<=Mnepmax)
	igas= Mp>Mnepmax
	
	# Rocky (also sets array)
	Rp= Mp**0.28
	disp= 0.04*Rp
	Rrockmax= Mrockmax**0.28
	
	# Neptunes
	Rp[inep]= Rrockmax * (Mp[inep]/Mrockmax)**0.59
	Rnepmax= Rrockmax * (Mnepmax/Mrockmax)**0.59
	disp[inep]= 0.15*Rp[inep]
	
	# Gas giants
	Rp[igas]= Rnepmax * (Mp[igas]/Mnepmax)**-0.044
	disp[igas]= 0.07*Rp[igas]
	
	return Rp, disp