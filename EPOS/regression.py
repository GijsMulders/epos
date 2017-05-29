'''
simple regression helper function
'''
import numpy as np

# simple sliding window, width is factor width (2)
def sliding_window_log(x, y, xgrid, width=2.):
	if y is None: y=np.ones_like(x)
	f= np.zeros_like(xgrid)
	for i, (_x, _y) in enumerate(zip(x,y)):
		#f+= np.where( 
		f[np.abs(np.log10(xgrid/_x))< (np.log10(width)/2.) ]+= _y
	return f
		