import numpy as np
import scipy.stats # norm
import EPOS.fitfunctions
from scipy.interpolate import interp2d, RectBivariateSpline

def periodradius(epos, Init=False, fpara=None, fdet=None, 
			xbin=None, ybin=None, xgrid=None, ygrid=None):
	''' return the period-radius distribution'''
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
	
	if fdet is not None:
		det_pdf= epos.func(epos.X_in, epos.Y_in, *fpar2d)*fdet
		det_pdf_X, det_pdf_Y= np.sum(pdf*fdet, axis=1), np.sum(pdf*fdet, axis=0)
	
	# calculate pdf on different grid?
	if xbin is not None:
		xgrid= np.logspace(np.log10(xbin[0]),np.log10(xbin[-1]),5)
	if ybin is not None:
		ygrid= np.logspace(np.log10(ybin[0]),np.log10(ybin[-1]),5)
	
	# calculate pdf on provided grid
	if (xgrid is not None) or (ygrid is not None):
		if xgrid is None:
			xgrid= epos.MC_xvar
		if ygrid is None:
			#ygrid= epos.MC_yvar
			ygrid= epos.in_yvar
	
		X,Y=np.meshgrid(xgrid, ygrid,indexing='ij')
		pdf= epos.func(X,Y, *fpar2d)
		#pdf_X, pdf_Y= np.sum(pdf, axis=1), np.sum(pdf, axis=0)
		
		# normalized per unit dlnxdlny
		pdf= pps* pdf/sum_pdf* epos.scale
		if fdet is not None:
			#func_fdet= interp2d(xgrid, ygrid, fdet.T, kind='cubic')
			func_fdet= RectBivariateSpline(epos.MC_xvar, epos.in_yvar, fdet)
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
			# normalized in total counts
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