import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb

	#calibrationLamp = [calRed, calBlue]
	#calPix = [calPixRed]
	#calAng = [calAngRed]
def AngPerPix(CalibrationLamp, calPix, calAng, show_plt=0):
	#cal pix and cal ang are values that determine, exactly,
	#what pixels were found to be contributing to certain 
	#wavelengths, as determined by an outside source
	


	openM = pyfits.open(CalibrationLamp)
	dataM = openM[0].data

	dataRowM = dataM[200,:]
	nPixel=openM[0].header['NAXIS1']

	xPixel = np.arange(0,nPixel) + 1

	#=======================Calibrate ==========================

	#calPix
	#calAng

	a, b, c = np.polyfit(calPix, calAng,2)
	print a, b, c
	pv = np.poly1d([a,b,c])
	fit = pv(xPixel)
	#fit = a*np.power(xPixel,2) + b*xPixel + c
	


	print 'halpha', pv(509.7)
	Lambda_x = fit
	
	
	
	
	#=========Calibration plot==============

	figCal = plt.figure()
	plt.plot(xPixel, dataRowM)
	plt.title("Spectrum Plot, ARC LAMP")
	plt.xlabel("pixels")
	plt.ylabel("Flux")
	if show_plt ==1:
		figCal.show()
		plt.savefig("calDM.eps")
		raw_input()
	plt.close()
	
	#=========Fit plot==============
	figCal = plt.figure()
	plt.plot(calPix, calAng, '.')
	plt.plot(xPixel, fit)
	plt.title("Poly Fit, 2nd order, ARC LAMP")
	plt.xlabel("pixels")
	plt.ylabel("wavelength")
	if show_plt == 1:
		plt.savefig("2_order_fit_arcLamp.eps")
		figCal.show()
		raw_input()
	plt.close()

	#=========spectrum plot==============
	figSpectrum = plt.figure()
	plt.plot(Lambda_x, dataRowM)
	plt.title("Spectrum, ARC LAMP")
	plt.xlabel("angstroms")
	if show_plt == 1:
		figSpectrum.show()
		plt.savefig("Spectrum_ArcLamp.eps")
		raw_input()
	plt.close()
	
	
	
	return Lambda_x

