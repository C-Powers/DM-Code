import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb
import astropy.io.fits as fits
import scipy.optimize
import DM_errors
import scipy.special

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
	
	perr = np.sqrt(np.diag(pcov))
	
	print "----------PERR", perr
	
	return gaussFit, indexData, x0, perr
	


def galaxyData(Galaxy, Lambda_x, subtractedErrors, show_plt=1):
	
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
	
	
	#+++DATA ROW 2 (around 169-167)
	indRow2=168
	dataRow21=dataD[169,:]
	dataRow22=dataD[168,:]
	dataRow23=dataD[167,:]
	rows2 = [dataRow21, dataRow22, dataRow23]

	
	#+++DATA ROW 3 (around 166-163)
	indRow3=165
	dataRow31=dataD[166,:]
	dataRow32=dataD[165,:]
	dataRow33=dataD[164,:]
	dataRow34=dataD[163,:]
	rows3 = [dataRow31, dataRow32, dataRow33, dataRow34]

	
	#+++DATA ROW 4 (162-159)
	indRow4 = 161
	dataRow41 = dataD[162,:]
	dataRow42 = dataD[161,:]
	dataRow43 = dataD[160,:]
	dataRow44 = dataD[159,:]
	rows4 = [dataRow41, dataRow42, dataRow43, dataRow44]

	#+++DATA ROW 5 (155-153)
	indRow5 = 154
	dataRow51 = dataD[155,:]
	dataRow52 = dataD[154,:]
	dataRow53 = dataD[153,:]
	rows5 = [dataRow51, dataRow52, dataRow53]
	
	

	#+++DATA ROW 6 (153-150)
	indRow6 = 152
	dataRow61 = dataD[153,:]
	dataRow62 = dataD[152,:]
	dataRow63 = dataD[151,:]
	dataRow64 = dataD[150,:]
	rows6 = [dataRow61, dataRow62, dataRow63, dataRow64]
	
	
	#+++DATA ROW 7 (146-144)
	indRow7 = 145
	dataRow71 = dataD[146,:]
	dataRow72 = dataD[145,:]
	dataRow73 = dataD[144,:]
	rows7 = [dataRow71, dataRow72, dataRow73]
	
	
	#+++DATA ROW 8 (144-141)
	indRow8 = 143
	dataRow81 = dataD[144,:]
	dataRow82 = dataD[143,:]
	dataRow83 = dataD[142,:]
	dataRow84 = dataD[141,:]
	rows8 = [dataRow81, dataRow82, dataRow83, dataRow84]

	
	
	#+++DATA ROW 9 (144-141)
	indRow9 = 143
	dataRow91 = dataD[144,:]
	dataRow92 = dataD[143,:]
	dataRow93 = dataD[142,:]
	dataRow94 = dataD[141,:]
	rows9 = [dataRow91, dataRow92, dataRow93, dataRow94]
	
	
	#+++DATA ROW 10 (138-136)
	indRow10 = 137
	dataRow101 = dataD[138,:]
	dataRow102 = dataD[137,:]
	dataRow103 = dataD[136,:]
	rows10 = [dataRow101, dataRow102, dataRow103]
	
	#+++DATA ROW 11 (136-133)
	indRow11 = 134
	dataRow111 = dataD[136,:]
	dataRow112 = dataD[135,:]
	dataRow113 = dataD[134,:]
	dataRow114 = dataD[133,:]
	rows11 = [dataRow111, dataRow112, dataRow113, dataRow114]
	
	
	#+++DATA ROW 12 (112-109)
	indRow12 = 111
	dataRow122 = dataD[112,:]
	dataRow123 = dataD[111,:]
	dataRow124 = dataD[110,:]
	dataRow125 = dataD[109,:]
	rows12 = [dataRow122, dataRow123, dataRow124, dataRow125]

	
	
	#+++DATA ROW 13 (101-99)
	indRow13=100
	dataRow132 = dataD[101,:]
	dataRow133 = dataD[100,:]
	dataRow134 = dataD[99,:]
	dataRow135 = dataD[98,:]
	rows13 = [dataRow132, dataRow133, dataRow134, dataRow135]

	
	#+++DATA ROW 14 (87-83)
	indRow14=85
	dataRow141 = dataD[87,:]
	dataRow142 = dataD[86,:]
	dataRow143 = dataD[85,:]
	dataRow144 = dataD[84,:]
	dataRow145 = dataD[83,:]
	rows14 = [dataRow141, dataRow142, dataRow143, dataRow144, dataRow145]
	
	

	#+++++++++ add rows+++++++++++
	dataRow1 = np.sum(rows1,0)
	dataRow2 = np.sum(rows2,0)
	dataRow3 = np.sum(rows3,0)
	dataRow4 = np.sum(rows4,0)
	dataRow5 = np.sum(rows5,0)
	dataRow6 = np.sum(rows6,0)
	dataRow7 = np.sum(rows7,0)
	dataRow8 = np.sum(rows8,0)
	dataRow9 = np.sum(rows9,0)
	dataRow10 = np.sum(rows10,0)
	dataRow11 = np.sum(rows11,0)
	dataRow12 = np.sum(rows12,0)
	dataRow13 = np.sum(rows13,0)
	dataRow14 = np.sum(rows14,0)
	
	
	dataRows = [dataRow1, dataRow2, dataRow3, dataRow4, dataRow5, dataRow6, dataRow7, dataRow8, dataRow9, dataRow10, dataRow11, dataRow12, dataRow13, dataRow14]

	#pdb.set_trace()
	
	
	
	figGalaxyRows = plt.figure()
	add=[100,200,300,400,500,600,700,800, 900, 1000, 1100, 1200, 1300, 1400]
	colors = ["red", "orange", "yellow", "green", "black", "blue", "purple", "black"]
	for i in range(0, len(dataRows)):
		plt.plot(Lambda_x, dataRows[i]/50+add[i])
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
 	errorGauss =[]
 	
 	

	for i in range(0, len(dataRows)):
		gaussFit, indexData, x0, sigma = gaussianFit(dataRows[i], Lambda_x, indices)
		gaussFits.append(gaussFit)
		dataIndexed.append(indexData)
		mu.append(x0)
		errorGauss.append(sigma)
	
 	print "sigma", errorGauss
 	print "mu", mu


	#errors_use=
	#pdb.set_trace()
	#ERRORS IN LAMBDA = .02
	lambdaErrorX = [.02]
	lambdaError = lambdaErrorX*len(Lambda_x)
	lambdaErrorIndexed = lambdaErrorX*len(indices)
	
	
	
	figFit = plt.figure()
	for i in range(0, len(dataIndexed)):
		yErrorsL = subtractedErrors[i]
		yErrors=np.array(yErrorsL)
		plt.plot(Lambda_x[indices], dataIndexed[i]/50+add[i], 'x')
		plt.errorbar(Lambda_x[indices], dataIndexed[i]/50+add[i], yerr=yErrors[indices]/50, xerr=lambdaErrorIndexed, fmt=None)
		plt.plot(xval, gaussFits[i]/50+add[i], color="red")
		#plt.axvline(mu[i])
	plt.title("Gauss Fit")
	plt.xlabel("Wavelength")
	plt.ylabel("Flux/50")
	if show_plt==1:
		plt.savefig("Gaussian_fits_regions.eps")
		figFit.show()
		raw_input()
	plt.close()
	
	figReg5 = plt.figure()
	yErrorsL = subtractedErrors[4]
	yErrors=np.array(yErrorsL)
	plt.plot(Lambda_x[indices], dataIndexed[4], 'x')
	plt.plot(xval, gaussFits[4])
	plt.errorbar(Lambda_x[indices], dataIndexed[4], yerr=yErrors[indices], fmt=None)
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
	plt.xlim(5500,7000)
	plt.ylim(-1000,3000)
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
	a2 = [6593.9, 6612.8]
	a2Lab = 6562.80
	aH = 6593.9
	
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
	print Z_sulfer, "redshift, sulfer1"
	
	Z_sulfer2 = (a3[1]/a3Lab[1])-1
	print Z_sulfer2, "redshift, sulfer2"
	
	Z_HA = (aH/a2Lab)-1
	print Z_HA, "redshift, Hydrogen Alpha"
	
	
	Zs = [Z_sulfer, Z_sulfer2, Z_HA]
	
	Z_med = np.median(Zs)
	
	print "Z med", Z_med
	
	Z_err = np.sqrt(np.power((Z_sulfer-Z_med),2) + np.power((Z_sulfer2-Z_med),2) + np.power(Z_HA-Z_med,2))
	
	print "Z_err", Z_err	
	
	
	
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
	

	row_indices = [indRow1, indRow2, indRow3, indRow4, indRow5, indRow6, indRow7, indRow8, indRow9, indRow10, indRow11, indRow12, indRow13, indRow14]
	arcSec_rows_array = []
	for i in range(0, len(row_indices)):
		asRow = ArcSec[row_indices]
		arcSec_rows_array.append(asRow)
	
	arcSec_rows = arcSec_rows_array[0]
	
	mpc_rows = .1*arcSec_rows-5
	
	#+++++Now, calculate the redshift for our mu values, assuming that row 128.0 is our central value
	
	Z_dop=[]
	velRot_mps=[]
	for i in range(0, len(mu)):
		z_doppler = mu[i]/mu[10] - 1
		z_use = z_doppler - Z_med
		vel_rot = (3*np.power(10,8))*z_doppler
		Z_dop.append(z_doppler)
		velRot_mps.append(vel_rot)
		
	velRot_kmps = velRot_mps/(np.power(10,3))
	
	
	
	
	
	
	
	
	######################################
	##Get ERRORS
	######################################
	#vErrors= []
	#errors in lambda are systematic, .02 angstroms
	#for i in range(0, len(yErrors)
	
	###Z errors
	z_error=[]
	zErrors = DM_errors.z_errors(lambdaError, mu)
	
	gError=[]
	for i in range(0, len(errorGauss)):
		valNow = errorGauss[i][3]
		gError.append(valNow)
		
		
	print "gError", gError
	
	
	vErrors = DM_errors.v_errors(zErrors)
	
	vErrors_kmpsX = vErrors * np.power(10,3)
	vErrors_kmps=vErrors_kmpsX[0]

	figRotCurve = plt.figure()
	for i in range(0, len(velRot_kmps)):
		plt.plot(mpc_rows[i], velRot_kmps[i], '.', color="black")
		plt.errorbar(mpc_rows[i], velRot_kmps[i], yerr=10*gError[i], ecolor="red",  fmt=None)
	plt.title("Rotation Curve")
	plt.xlabel("Distance (Mpc)")
	plt.ylabel("Velocity (km/s)")
	if show_plt==1:
		plt.savefig("rotCurve.eps")
		figRotCurve.show()
	raw_input()
	plt.close()
	
	
	figRotCurveHalf = plt.figure()
	for i in range(0, len(velRot_kmps)):
		plt.plot(velRot_kmps[i], mpc_rows[i], '.', color="black")
		plt.errorbar(velRot_kmps[i], mpc_rows[i], xerr=3*gError[i], ecolor="red", fmt=None)
	plt.title("Rotation Curve")
	plt.ylabel("Distance (Mpc)")
	plt.xlabel("Velocity (km/s)")
	plt.ylim(0, 4)
	if show_plt==1:
		plt.savefig("rotCurveHalf.eps")
		figRotCurveHalf.show()
	raw_input()
	plt.close()
	
	
	###save text files
	np.savetxt("velRot_kmps.txt", velRot_kmps)
	np.savetxt("mpc_rows.txt", mpc_rows)
	np.savetxt("error_vel.txt", gError)
	
	
	
	
	
	
	#######################
	### BEssel function stuff?
	########################
	R=mpc_rows
	
	G = 6.67/(10**11)
	Rd = np.power(3.08*np.power(10,16),3.5)
	E = (2*10**40)/(2*np.pi*(Rd)**2)
	
	I0 = scipy.special.iv(0, R/(2*Rd))
	I1 = scipy.special.iv(1, R/(2*Rd))
	
	
	K0 = scipy.special.kn(0, R/(2*Rd))
	K1 = scipy.special.kn(1, R/(2*Rd))

		
	v2 = -2*np.pi*G*E*Rd*np.power(R/(2*Rd),2)*(I0*K1 - I1*K0)
	
	v = np.sqrt(v2)
	
	print v, "vcurve"
	
	pdb.set_trace()
	
	figRotCurveFit = plt.figure()
	for i in range(0, len(velRot_kmps)):
		plt.plot(velRot_kmps[i], mpc_rows[i], '.', color="black")
		plt.errorbar(velRot_kmps[i], mpc_rows[i], xerr=3*gError[i], ecolor="red", fmt=None)
		plt.plot(mpc_rows[i], v[i], "o", color="blue")
	plt.title("Rotation Curve")
	plt.ylabel("Distance (Mpc)")
	plt.xlabel("Velocity (km/s)")
	plt.ylim(0, 4)
	if show_plt==1:
		plt.savefig("rotCurveFit.eps")
		figRotCurveFit.show()
	raw_input()
	plt.close()
	
	
	
	

	
