import numpy as np
from scipy.stats import norm

def readme():
	print '\nThis module contains fitting functions for parametric fits to Kepler data'
	print
	print 'These functions work on numpy arrays ONLY'

def powerlaw2D(x, y, p0, p1, p2):
	return p0 * x**p1 * y**p2

def powerlaw2D_yonly(x, y, p0, p1):
	return p0 * y**p1

def brokenpowerlaw2D(x, y, p0, xp, p1, p2, yp, p3, p4):
	return p0 * brokenpowerlaw1D(x, xp, p1, p2) * brokenpowerlaw1D(y, yp, p3, p4)

def brokenpowerlaw2D_yonly(x, y, p0, p1, yp, p3, p4):
	return p0 * x**p1 * brokenpowerlaw1D(y, yp, p3, p4)

	
def brokenpowerlaw1D(x, xp, p1, p2):
	# normalized to 1 at x=xp
	return (x/xp)**np.where(x<xp,p1,p2)

def bimodal2D(x, y, rp0, rxp, rp1, rp2, ry, ryw, gp0, gxp, gp1, gp2, gy, gyw):
	rocky= rp0 * brokenpowerlaw1D(x, rxp, rp1, rp2) * \
			norm.pdf(np.log10(y), loc=np.log10(ry), scale=ryw)
	gas=   gp0 * brokenpowerlaw1D(x, gxp, gp1, gp2) * \
			norm.pdf(np.log10(y), loc=np.log10(gy), scale=gyw)
	return rocky+gas
