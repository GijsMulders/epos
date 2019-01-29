'''This module contains some fitting functions for parametric fits to Kepler data
These functions work on numpy arrays ONLY, 
all functions return a non-normalized probability density distribution, 
assuming x and y are numpy arrays of equal length (i.e. coordinates)
'''

import numpy as np
import scipy.stats

def uniform(x,y): 
	''' Uniform distribution'''
	return np.ones_like(x)

def powerlaw2D(x, y, p1, p2):
	''' Power law in x and y'''
	return x**p1 * y**p2

def powerlaw2D_yonly(x, y, p1):
	''' Power law in y, uniform distribution in x'''
	return y**p1

def brokenpowerlaw2D(x, y, xp, p1, p2, yp, p3, p4):
	''' Broken powerlaw in x and y
	
	Args:
		x(np.array): x
		y(np.array): y
		xp(float): break in x
		p1(float): power law index at x<xp
		p2(float): power law index at x>xp
		yp(float): break in y
		p3(float): power law index at y<yp
		p4(float): power law index at y>yp
	'''
	return brokenpowerlaw1D(x, xp, p1, p2) * brokenpowerlaw1D(y, yp, p3, p4)

def brokenpowerlaw2D_xonly(x, y, xp, p1, p2, p3):
	return brokenpowerlaw1D(x, xp, p1, p2) * y**p3

def brokenpowerlaw2D_yonly(x, y, p1, yp, p3, p4):
	return x**p1 * brokenpowerlaw1D(y, yp, p3, p4)

def brokenpowerlaw2D_symmetric(x, y, xp, p1, yp, p3, p4):
	return brokenpowerlaw2D(x, y, xp, p1, -p1, yp, p3, p4)

def doublebrokenpowerlaw2D(x, y, xp, p1, p2, yp, p3, p4, yr, yrp, p5, p6):
	''' Broken powerlaw in x, double broken power-law in y
	
	Args:
		x(np.array): x
		y(np.array): y
		xp(float): break in x
		p1(float): power law index at x<xp
		p2(float): power law index at x>xp
		yp(float): break in y
		p3(float): power law index at y<yp
		p4(float): power law index at y>yp
		yr(float): ratio between power-law normalizations
		yrp(float): break in y
		p5(float): power law index at y<yp
		p6(float): power law index at y>yp		
	'''
	xfunc= brokenpowerlaw1D(x, xp, p1, p2)
	yfunc= brokenpowerlaw1D(y, yp, p3, p4) + yr* brokenpowerlaw1D(y, yrp, p5, p6) 
	return xfunc* yfunc
	
def brokenpowerlaw1D(x, xp, p1, p2):
	# normalized to 1 at x=xp
	return (x/xp)**np.where(x<xp,p1,p2)

def lognormal_size(x, y, xp, p1, p2, y0, dy):
	'''
	Lognormal distribution in planet size `y` and a broken power law in distance `x` 
	
	Args:
		x(np.array): x
		y(np.array): y
		xp(float): break in x
		p1(float): power law index at x<xp
		p2(float): power law index at x>xp
		y(float): mean of y
		dy(flat): dispersion of y, in dex	
	Note:		
		Only works with support for a non-separable function in x and y
	'''	
	return brokenpowerlaw1D(x, xp, p1, p2) * \
			scipy.stats.norm.pdf(np.log10(y), loc=np.log10(y0), scale=dy)

def lognormal_rise(x, y, y0, dy, yp):
	'''
	Lognormal distribution in planet size `y` that rises with distance `x` 
	
	Args:
		x(np.array): x
		y(np.array): y
		y(float): mean of y
		dy(flat): dispersion of y, in dex	
		yp(float): break in x
	Note:		
		Only works with support for a non-separable function in x and y
	'''	
	return scipy.stats.norm.pdf(np.log10(y), loc=np.log10(y0* x**yp), scale=dy)

def schechter1D(x, xp, p1):
	''' Schechter function normalized to p break'''
	return (x/xp)**p1 * np.exp(-x/xp)

def schechter_size(x, y, xp, p1, p2, yp, p3):
	'''
	Schechter distribution in planet size `y` and a broken power law in distance `x` 
	
	Args:
		x(np.array): x
		y(np.array): y
		xp(float): break in x
		p1(float): power law index at x<xp
		p2(float): power law index at x>xp
		y(float): mean of y
		dy(flat): dispersion of y, in dex	
	Note:		
		Only works with support for a non-separable function in x and y
	'''	
	return brokenpowerlaw1D(x, xp, p1, p2) * \
			schechter1D(y,yp,p3)

def bimodal2D(x, y, rp0, rxp, rp1, rp2, ry, ryw, gxp, gp1, gp2, gy, gyw):
	''' Bimodal distribution where each component (r,g) is a broken power law in `x` 
	and a lognormal in `y`
	
	Args:
		x(np.array): x
		y(np.array): y
		rp0(float): ratio of rocky `r` component versus gaseous `g` component
		rxp(float): break in x
		rp1(float): power law index at x<xp
		rp2(float): power law index at x>xp
		ry(float): mean of y
		ryw(flat): dispersion of y, in dex	
	Note:		
		Only works with support for a non-separable function in x and y
	'''
	rocky= rp0 * lognormal_size(x, y, rxp, rp1, rp2, ry, ryw)
	gas=   lognormal_size(x, y, gxp, gp1, gp2, gy, gyw)
	return rocky+gas
