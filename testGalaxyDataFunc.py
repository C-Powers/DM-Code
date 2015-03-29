import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb
import astropy.io.fits as fits
import scipy.optimize


def gauss_off(x, Z, A, x0, sigma):
   	return Z + A*np.exp(- (x-x0)**2 / (2.*sigma**2) )
   	
def gaussianFit(dataRow, Lambda_x, indices):
	# Guess
	Aguess = np.max(dataRow[indices])
	sguess = 2.
	xguess = np.median(Lambda_x[indices])
	Z = np.median(dataRow[indices])
	pguess = [Z, Aguess, xguess, sguess]

	# Fit
	popt, pcov = scipy.optimize.curve_fit(gauss_off, Lambda_x[indices], dataRow[indices], p0=pguess)
	#pdb.set_trace()
	Z, A, x0, sigma = popt
	mn = np.min(Lambda_x[indices])
	mx = np.max(Lambda_x[indices])
	xval = np.arange(mn,mx,0.1)
	gaussFit = gauss_off(xval, Z, A, x0, sigma)

	indexData = dataRow[indices]

	return gaussFit, indexData, x0
	


def galaxyData(Galaxy, Lambda_x, show_plt=1):
	
	openD = pyfits.open(Galaxy)
	dataD = openD[0].data
	dataRowD = dataD[143,:]
	


	
	
	figGalaxy = plt.figure()
	plt.plot(Lambda_x,dataRowD)
	plt.title("Single Row of Galaxy, Sky Subtracted Data")
	plt.xlabel("Wavelength (Angstroms)")
	plt.ylabel("Flux")
	if show_plt==1:
		figGalaxy.show()
		raw_input()
	plt.close()
	
	

	
	###################
	#data rows
	###################
	
	run_row = np.arange(80.0,172.0)
	

	dataRowsArray = np.zeros((len(run_row), 1232))
	for i in range(0, len(run_row)):
		rowy = dataD[run_row[i],:]
		#dataRowsArray = np.append(dataRowsArray,rowy)	
		dataRowsArray[i,:] = rowy

	dataRows= np.zeros((len(dataRowsArray[0]), len(dataRowsArray[1])))
	for i in range(0, len(dataRowsArray-3)):
		start=3*i
		sum = np.sum(dataRowsArray[start:3+start],1)
	#	dataRows=np.append(dataRows, sum)
		dataRows[:,i] = sum
	
		
	pdb.set_trace()
		
	#========================Fit the data rows to a gaussian===================
	# Define a Gaussian plus a floor offset
 
 	indices = np.arange(505,517)
 	mn = np.min(Lambda_x[indices])
	mx = np.max(Lambda_x[indices])
	xval = np.arange(mn,mx,0.1)
	
 	 	
 	gaussFits = []
 	dataIndexed = []
 	mu = []
 	
 	for i in range(0, len(dataRows)):
 		gaussFit, indexData, x0 = gaussianFit(dataRows[i], Lambda_x, indices)
 		gaussFits.append(gaussFit)
 		dataIndexed.append(indexData)
 		mu.append(x0)

	colors = ["maroon", "red", "orange", "green", "cyan", "blue", "purple", "black"]
	figFit = plt.figure()
	for i in range(0, len(dataIndexed)):
		plt.plot(Lambda_x[indices], dataIndexed[i], 'x')
		plt.plot(xval, gaussFits[i])
		plt.axvline(mu[i])
	plt.title("Gauss Fit")
	plt.xlabel("Wavelength")
	plt.ylabel("Flux/50")
	if show_plt==1:
		plt.savefig("Gaussian_fits_regions.eps")
		figFit.show()
		raw_input()
	plt.close()
	
	figReg5 = plt.figure()
	plt.plot(Lambda_x[indices], dataIndexed[4], 'x')
	plt.plot(xval, gaussFits[4])
	plt.title("Gaussian Fit, Region 5")
	plt.xlabel("Wavelength")
	plt.ylabel("Flux")
	if show_plt==1:
		plt.savefig("Gaussian_fit5_region5.eps")
		figReg5.show()
	raw_input()
	plt.close()
	
	
	
	
	##########################################
	#============FIND REDSHIFT============
	##########################################
	
	
	#++++find two lines that we can know, for sure, we have
	#left pair
	a1 = [5776.34, 5785.7]
	
	#H_alpha and comrade (middle pair)
	a2 = [6593.86, 6612.78]
	a2Lab = 6562.80
	
	#right pair (sulfer?)
	a3 = [6747.45, 6761.63]
	a3Lab = [6716.47, 6730.17]
	
	deltaData = a3[0]-a3[1]
	aveData = (a3[0]+a3[1])/2.0
	ratioData = deltaData/aveData
	
	deltaLab = a3Lab[0]-a3Lab[1]
	aveLab = (a3Lab[0]+a3Lab[1])/2.0
	ratioLab = deltaLab/aveLab
	print ratioLab

	# ratio = ratioData/ratioLab ++ comes out to be .922
	
	###++++++++++++++++ NOW CALCULATE REDSHIFT ++++++++++++++++++
	
	Z_sulfer = (a3[0]/a3Lab[0])-1
	print Z_sulfer, "redshift, sulfer"
	
	Z_HA = (a2[0]/a2Lab)-1
	print Z_HA, "redshift, Hydrogen Alpha"
	
	
	
	
	
	
	##############################################################
	#Get a velocity distribution function
	##############################################################



	
	
	#length of slit, 2"
	#assume galaxy is basically the length of slit.
	#galaxy extends from [64.0, 192.0] on the CCD detector
	
	detRange = [64.0, 192.0]
	detRangeVal = 192.0-64.0
	print detRangeVal, "det range"
	medDetRange = np.median(detRange)
	print medDetRange, "med"
	
	#middle of the galaxy, where there should be no observable rotation, is at row 128.0
	
	#+++++arcsecs as a function of row
	nPixelRow=openD[0].header['NAXIS2']
	xRow = np.arange(0,nPixelRow) + 1
	
	ArcSec = (.78)*xRow - 49.92
	
	#+++++Now, get the arcsecond values for rows 1-8
	

	arcSec_rows_array = []
	for i in range(0, len(run_row)):
		asRow = ArcSec[run_row[i]]
		arcSec_rows_array.append(asRow)
	
	arcSec_rows = arcSec_rows_array
	
	#+++++Now, calculate the redshift for our mu values, assuming that row 128.0 is our central value
	
	Z_dop=[]
	velRot_mps=[]
	for i in range(0, len(mu)):
		z_doppler = mu[i]/mu[50] - 1
		vel_rot = (3*np.power(10,8))*z_doppler
		Z_dop.append(z_doppler)
		velRot_mps.append(vel_rot)
		
	velRot_kmps = velRot_mps/(np.power(10,3))
	
	
	#pdb.set_trace()
	
	figRotCurve = plt.figure()
	for i in range(0, len(velRot_kmps)):
		plt.plot(arcSec_rows[i], velRot_kmps[i], 'o', color="black")
	plt.title("Rotation Curve")
	plt.xlabel("Distance (arcsec)")
	plt.ylabel("Velocity (km/s)")
	if show_plt==1:
		plt.savefig("rotCurve.eps")
		figRotCurve.show()
	raw_input()
	plt.close()
	
	
	
	
	
	###############################################################
	#GET A VEL DIST FUNC BY TAKING EVERY ROW
	###############################################################
	
	

	
