import numpy as np
import scipy.stats # norm
import EPOS.fitfunctions
from scipy.interpolate import interp2d, RectBivariateSpline

def periodradius(epos, Init=False, fpara=None, fdet=None, 
			xbin=None, ybin=None, xgrid=None, ygrid=None, Convert=False):
	''' 
	Return the period-radius distribution
	in referse to the model parameter (Mass, Msini or radius)
	out reverse to the observed parameter (Msini or Radius)
	default is the input grid. 
	For MC=False, the distribution can also be Converted onto the ouput grid (M->Msini or M->R)
	'''
	if fpara is None:
		pps= epos.pdfpars.getpps(Init=Init)
		fpar2d= epos.pdfpars.get2d(Init=Init)
	else:
		pps= epos.pdfpars.getpps_fromlist(fpara)
		fpar2d= epos.pdfpars.get2d_fromlist(fpara)
		#print fpara
	
	pdf= epos.func(epos.X_in, epos.Y_in, *fpar2d)
	pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
	sum_pdf= np.sum(pdf)
	sum_pdf_X= np.sum(pdf_X)
	sum_pdf_Y= np.sum(pdf_Y)

	''' Convert M->Msini or M->R'''
	if Convert:
		if epos.RV:
			pdf= _convert_pdf_M_to_Msini(epos.MC_xvar, epos.in_yvar, epos.MC_yvar, pdf)
			pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
		elif epos.MR:
			pass

	# convolve with detection efficiency?
	if fdet is not None:
		det_pdf= pdf*fdet
		det_pdf_X, det_pdf_Y= np.sum(det_pdf, axis=1), np.sum(det_pdf, axis=0)
	
	# calculate pdf on different grid?
	if xbin is not None:
		xgrid= np.logspace(np.log10(xbin[0]),np.log10(xbin[-1]),5)
	if ybin is not None:
		ygrid= np.logspace(np.log10(ybin[0]),np.log10(ybin[-1]),5)
	
	# calculate pdf on provided grid
	if (xgrid is not None) or (ygrid is not None):
		#assert Convert==False, 'Not implemented'
		if xgrid is None:
			xgrid= epos.MC_xvar
		if ygrid is None:
			#ygrid= epos.MC_yvar
			ygrid= epos.in_yvar
	
		X,Y=np.meshgrid(xgrid, ygrid,indexing='ij')
		pdf= epos.func(X,Y, *fpar2d)
		
		# normalized per unit dlnxdlny
		pdf= pps* pdf/sum_pdf* epos.scale

		if Convert:
			assert np.all(ygrid==epos.MC_yvar) or np.all(ygrid==epos.in_yvar), 'cannot do custom y grid'
			pdf= _convert_pdf_M_to_Msini(xgrid, ygrid, epos.MC_yvar, pdf)
			#pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)

		if fdet is not None:
			#assert False, 'Where is this used?'
			#func_fdet= interp2d(xgrid, ygrid, fdet.T, kind='cubic')
			func_fdet= RectBivariateSpline(epos.MC_xvar, epos.MC_yvar, fdet) # was in_yvar
			_fdet= func_fdet(xgrid, ygrid)
			#print fdet.shape, _fdet.shape
			#print xgrid
			#print ygrid
			#print _fdet
			pdf *= _fdet

		dlnx= np.log(xgrid[-1]/xgrid[0])
		dlny= np.log(ygrid[-1]/ygrid[0])

		pdf_X= np.average(pdf, axis=1) * dlny
		pdf_Y= np.average(pdf, axis=0) * dlnx

	else:
		# calculate pdf on default grid
		# normalized in units dlnx, dlny, and dlnxdlny
		if fdet is None:
			pdf_X= pps* pdf_X/sum_pdf_X * epos.scale_x
			pdf_Y= pps* pdf_Y/sum_pdf_Y * epos.scale_in_y
			pdf= pps* pdf/sum_pdf* epos.scale
		else:
			pdf_X= pps* det_pdf_X/sum_pdf_X * epos.scale_x
			pdf_Y= pps* det_pdf_Y/sum_pdf_Y * epos.scale_in_y
			pdf= pps* det_pdf/sum_pdf* epos.scale
	
	return pps, pdf, pdf_X, pdf_Y

def periodratio(epos, Pgrid=None, fpara=None, Init=False):
	
	if fpara is None:
		if epos.spacing == 'powerlaw':
			dPbreak= epos.fitpars.get('dP break')
			dP1= epos.fitpars.get('dP 1')
			dP2= epos.fitpars.get('dP 2')
			pdf= EPOS.fitfunctions.brokenpowerlaw1D(Pgrid, dPbreak, dP1, dP2)
		elif epos.spacing=='dimensionless':
			logD=  epos.fitpars.get('log D')
			sigma= epos.fitpars.get('sigma')

			with np.errstate(divide='ignore'): 
				Dgrid= np.log10(2.*(Pgrid**(2./3.)-1.)/(Pgrid**(2./3.)+1.))
			Dgrid[0]= -2
			#print Dgrid
			pdf= scipy.stats.norm(logD,sigma).pdf(Dgrid)
			
		cdf= np.cumsum(pdf)
		cdf/=cdf[-1]
	else:
		# copy-pasta from run.py
		raise ValueError('test needed')
	
	return pdf, cdf

''' Helper function '''
def _int_fmsini(x):
	''' Integral of x/sqrt(1.-x**2.) to avoid infinite value at x=1'''
   	return -np.sqrt(1.-x**2.)

def _convert_pdf_M_to_Msini(xgrid, in_ygrid, out_ygrid, pdf):
	'''Generate the shape of the Msini distribution'''
	fint1d= _int_fmsini(in_ygrid[1:]/in_ygrid[-1])-_int_fmsini(in_ygrid[:-1]/in_ygrid[-1])
	fint2d= np.tile(fint1d, (xgrid.size,1))
	#print 'fint 1d {}'.format(fint1d.shape)
	#print 'fint 2d {}'.format(fint2d.shape)
	#print 'pdf {}'.format(pdf.shape)

	''' Convolve the the pdf with the msini distribution through shift and add'''	
	pdf_ext= np.zeros_like(pdf) # extended grid
	for i in range(in_ygrid.size-1):
		fobsi= pdf[:,i][:, np.newaxis]*fint2d[:,-(i+1):] # matrix multiplication.  
		pdf_ext[:, :(i+1)]+= fobsi

	''' Truncate the pdf to the mass range of Msini'''
	pdf_out= pdf_ext[:,:out_ygrid.size] # this in on out grid
	#print 'pdf before {}, after {}'.format(pdf.shape, pdf_ext[:,:out_ygrid.size].shape)
	return pdf_out 