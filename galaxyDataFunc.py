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
	
	
	
	#now, look at how the rows correspond to different wavelengths
	#dataRow1 = dataD[167,:]
	#dataRow2 = dataD[160,:]
	#dataRow3 = dataD[151,:]
	#dataRow4 = dataD[143,:]
	#dataRow5 = dataD[135,:]
	#dataRow6 = dataD[111,:]
	#dataRow7 = dataD[100,:]
	#dataRow8 = dataD[85,:]
	
	
	#============================ Create data (accumulation of rows) to use for SCIENCE===============
	
	#+++DATA ROW 1 (around 172-170)
	indRow1 = 171
	dataRow11 = dataD[172,:]
	dataRow12 = dataD[171,:]
	dataRow13 = dataD[170,:]
	
	rows1 = [dataRow11, dataRow12, dataRow13]
	
	#+++DATA ROW 2 (162-159)
	indRow2 = 161
	dataRow21 = dataD[162,:]
	dataRow22 = dataD[161,:]
	dataRow23 = dataD[160,:]
	dataRow24 = dataD[159,:]
	
	rows2 = [dataRow21, dataRow22, dataRow23, dataRow24]

	
	#+++DATA ROW 3 (153-150)
	indRow3 = 152
	dataRow31 = dataD[153,:]
	dataRow32 = dataD[152,:]
	dataRow33 = dataD[151,:]
	dataRow34 = dataD[150,:]
	
	rows3 = [dataRow31, dataRow32, dataRow33, dataRow34]

	
	#+++DATA ROW 4 (146-142)
	indRow4=144
	dataRow41 = dataD[146,:]
	dataRow42 = dataD[145,:]
	dataRow43 = dataD[144,:]
	dataRow44 = dataD[143,:]
	dataRow45 = dataD[142,:]
	
	rows4 = [dataRow41, dataRow42, dataRow43, dataRow44, dataRow45]

	
	
	#+++DATA ROW 5 (138-133)
	indRow5 = 136
	dataRow51 = dataD[138,:]
	dataRow52 = dataD[137,:]
	dataRow53 = dataD[136,:]
	dataRow54 = dataD[135,:]
	dataRow55 = dataD[134,:]
	dataRow56 = dataD[133,:]
	
	rows5 = [dataRow51, dataRow52, dataRow53, dataRow54, dataRow55, dataRow56]

	
	
	#+++DATA ROW 6 (113-108)
	indRow6 = 111
	dataRow61 = dataD[113,:]
	dataRow62 = dataD[112,:]
	dataRow63 = dataD[111,:]
	dataRow64 = dataD[110,:]
	dataRow65 = dataD[109,:]
	dataRow66 = dataD[108,:]
	
	rows6 = [dataRow61, dataRow62, dataRow63, dataRow64, dataRow65, dataRow66]

	
	
	#+++DATA ROW 7 (102-97)
	indRow7=100
	dataRow71 = dataD[102,:]
	dataRow72 = dataD[101,:]
	dataRow73 = dataD[100,:]
	dataRow74 = dataD[99,:]
	dataRow75 = dataD[98,:]
	dataRow76 = dataD[97,:]
	
	rows7 = [dataRow71, dataRow72, dataRow73, dataRow74, dataRow75, dataRow76]

	
	
	#+++DATA ROW 8 (87-83)
	indRow8=85
	dataRow81 = dataD[87,:]
	dataRow82 = dataD[86,:]
	dataRow83 = dataD[85,:]
	dataRow84 = dataD[84,:]
	dataRow85 = dataD[83,:]
	
	rows8 = [dataRow81, dataRow82, dataRow83, dataRow84, dataRow85]
	
	

	#+++++++++ add rows+++++++++++
	dataRow1 = np.sum(rows1,0)
	dataRow2 = np.sum(rows2,0)
	dataRow3 = np.sum(rows3,0)
	dataRow4 = np.sum(rows4,0)
	dataRow5 = np.sum(rows5,0)
	dataRow6 = np.sum(rows6,0)
	dataRow7 = np.sum(rows7,0)
	dataRow8 = np.sum(rows8,0)
	
	dataRows = [dataRow1, dataRow2, dataRow3, dataRow4, dataRow5, dataRow6, dataRow7, dataRow8]

	#pdb.set_trace()
	
	
	
	figGalaxyRows = plt.figure()
	add=[100,200,300,400,500,600,700,800]
	colors = ["red", "orange", "yellow", "green", "black", "blue", "purple", "black"]
	for i in range(0, len(dataRows)):
		plt.plot(Lambda_x, dataRows[i]/50+add[i], color=colors[i])
	plt.title("Galaxy - Multiple Rows")
	plt.xlabel("Wavelength")
	plt.ylabel("Flux/50")
	if show_plt==1:
		figGalaxyRows.show()
		raw_input()
	plt.close()
	
		
	#========================Fit the data rows to a gaussian===================
	# Define a Gaussian plus a floor offset
 
 	indices = np.arange(505,518)
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




	figFit = plt.figure()
	for i in range(0, len(dataIndexed)):
		plt.plot(Lambda_x[indices], dataIndexed[i]/50+add[i], 'x', color = colors[i])
		plt.plot(xval, gaussFits[i]/50+add[i], color = colors[i])
		plt.axvline(mu[i], color=colors[i])
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
	
	#use dataRow5
	figSpecPlt = plt.figure()
	plt.plot(Lambda_x, dataRow5)
	plt.axvline(6747.45, color="red")
	plt.axvline(6761.63, color="red")
	plt.axvline(6716.47, color="green")
	plt.axvline(6730.17, color="green")
	plt.title("Spec Plot, row5")
	plt.xlabel("Wavelength")
	plt.ylabel("Flux/50")
	plt.savefig("Spectral_plot_region5.eps")
	figSpecPlt.show()
	raw_input()
	plt.close()
	
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
	
	print "ratioLab", ratioLab
	print "ratioData", ratioData

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
	

	row_indices = [indRow1, indRow2, indRow3, indRow4, indRow5, indRow6, indRow7, indRow8]
	arcSec_rows_array = []
	for i in range(0, len(row_indices)):
		asRow = ArcSec[row_indices]
		arcSec_rows_array.append(asRow)
	
	arcSec_rows = arcSec_rows_array[0]
	
	#+++++Now, calculate the redshift for our mu values, assuming that row 128.0 is our central value
	
	Z_dop=[]
	velRot_mps=[]
	for i in range(0, len(mu)):
		z_doppler = mu[i]/mu[4] - 1
		vel_rot = (3*np.power(10,8))*z_doppler
		Z_dop.append(z_doppler)
		velRot_mps.append(vel_rot)
		
	velRot_kmps = velRot_mps/(np.power(10,3))
	
	figRotCurve = plt.figure()
	for i in range(0, len(velRot_kmps)):
		plt.plot(arcSec_rows[i], velRot_kmps[i], 'o', color=colors[i])
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
	
	

	
