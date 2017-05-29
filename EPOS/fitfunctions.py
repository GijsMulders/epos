import numpy as np

def readme():
	print '\nThis module contains fitting functions for parametric fits to Kepler data'
	print
	print 'These functions work on numpy arrays ONLY'
	
def brokenpowerlaw1D(x, xp, p1, p2):
	# normalized to 1 at x=xp
	return (x/xp)**np.where(x<xp,p1,p2)
