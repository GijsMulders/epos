import numpy as np

def readme():
	print '\nThis module contains fitting functions for parametric fits to Kepler data'
	print
	print 'These functions work on numpy arrays ONLY'

def brokenpowerlaw2D(x, y, p0, xp, p1, p2, yp, p3, p4):
	return p0 * brokenpowerlaw1D(x, xp, p1, p2) * brokenpowerlaw1D(y, yp, p3, p4)
	
def brokenpowerlaw1D(x, xp, p1, p2):
	# normalized to 1 at x=xp
	return (x/xp)**np.where(x<xp,p1,p2)
