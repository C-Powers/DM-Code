import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb
import astropy.io.fits as fits
import DM_errors


def skySubtraction(dataFile, show_plt=0):


	open = pyfits.open(dataFile)
	data = open[0].data
	
	
	nPixel=open[0].header['NAXIS1']
	xPixel = np.arange(0,nPixel) + 1


	#create a few rows to then average and subtract out of the data file
	dataRow1 = data[60,:]
	dataRow2 = data[65,:]
	dataRow3 = data[223,:]
	dataRow4 = data[225,:]

	#dataRow = [dataRow1, dataRow2, dataRow3, dataRow4]
	#USE THESE, IT SEEMS AS THOUGH THE GALAXY EXTENDS TO THE BOTTOM OF THE IMAGE
	dataRow = [dataRow3, dataRow4]

	dataRowMedian = np.median(dataRow, 0)


		
	figCal = plt.figure()
	plt.plot(xPixel, dataRow1)
	plt.plot(xPixel, dataRow2)
	plt.plot(xPixel, dataRow3)
	plt.plot(xPixel, dataRow4)
	plt.title("Spectrum Plot of Individual Rows, 4 chosen rows")
	plt.xlabel("pixels")
	plt.ylabel("Flux")
	if show_plt ==1:
		figCal.show()
		raw_input()
	plt.close()
	
	
	figCal1 = plt.figure()
	plt.plot(xPixel, dataRowMedian)
	plt.title("Spectrum Plot of Individual Rows (median)")
	plt.xlabel("pixels")
	plt.ylabel("Flux")
	if show_plt ==1:
		figCal1.show()
		raw_input()
	plt.close()
	
	
	
	
	#=========SKY SUBTRACT===========
	#dataToSu= data[166,:]
	dataToSubtract = np.tile(dataRowMedian, (250,1))
	skySubtract = data - dataToSubtract

		
	
	#print "ALMOST DONE"
	fits.writeto("skySubtract.fits", skySubtract, clobber=True)
	#print "DONE"
	
	
	#===============================================
	# Get errors on non-sky subtracted data
	#===============================================


	#============================ Create data (accumulation of rows) to use for SCIENCE===============
	
	#+++DATA ROW 1 (around 172-170)
	indRow1 = 171
	dataRow11 = data[172,:]
	dataRow12 = data[171,:]
	dataRow13 = data[170,:]
	
	rows1 = [dataRow11, dataRow12, dataRow13]
	
	
	#+++DATA ROW 2 (around 169-167)
	indRow2=168
	dataRow21=data[169,:]
	dataRow22=data[168,:]
	dataRow23=data[167,:]
	rows2 = [dataRow21, dataRow22, dataRow23]

	
	#+++DATA ROW 3 (around 166-163)
	indRow3=165
	dataRow31=data[166,:]
	dataRow32=data[165,:]
	dataRow33=data[164,:]
	dataRow34=data[163,:]
	rows3 = [dataRow31, dataRow32, dataRow33, dataRow34]

	
	#+++DATA ROW 4 (162-159)
	indRow4 = 161
	dataRow41 = data[162,:]
	dataRow42 = data[161,:]
	dataRow43 = data[160,:]
	dataRow44 = data[159,:]
	rows4 = [dataRow41, dataRow42, dataRow43, dataRow44]

	#+++DATA ROW 5 (155-153)
	indRow5 = 154
	dataRow51 = data[155,:]
	dataRow52 = data[154,:]
	dataRow53 = data[153,:]
	rows5 = [dataRow51, dataRow52, dataRow53]
	
	

	#+++DATA ROW 6 (153-150)
	indRow6 = 152
	dataRow61 = data[153,:]
	dataRow62 = data[152,:]
	dataRow63 = data[151,:]
	dataRow64 = data[150,:]
	rows6 = [dataRow61, dataRow62, dataRow63, dataRow64]
	
	
	#+++DATA ROW 7 (146-144)
	indRow7 = 145
	dataRow71 = data[146,:]
	dataRow72 = data[145,:]
	dataRow73 = data[144,:]
	rows7 = [dataRow71, dataRow72, dataRow73]
	
	
	#+++DATA ROW 8 (144-141)
	indRow8 = 143
	dataRow81 = data[144,:]
	dataRow82 = data[143,:]
	dataRow83 = data[142,:]
	dataRow84 = data[141,:]
	rows8 = [dataRow81, dataRow82, dataRow83, dataRow84]

	
	
	#+++DATA ROW 9 (144-141)
	indRow9 = 143
	dataRow91 = data[144,:]
	dataRow92 = data[143,:]
	dataRow93 = data[142,:]
	dataRow94 = data[141,:]
	rows9 = [dataRow91, dataRow92, dataRow93, dataRow94]
	
	
	#+++DATA ROW 10 (138-136)
	indRow10 = 137
	dataRow101 = data[138,:]
	dataRow102 = data[137,:]
	dataRow103 = data[136,:]
	rows10 = [dataRow101, dataRow102, dataRow103]
	
	#+++DATA ROW 11 (136-133)
	indRow11 = 134
	dataRow111 = data[136,:]
	dataRow112 = data[135,:]
	dataRow113 = data[134,:]
	dataRow114 = data[133,:]
	rows11 = [dataRow111, dataRow112, dataRow113, dataRow114]
	
	
	#+++DATA ROW 12 (112-109)
	indRow12 = 111
	dataRow122 = data[112,:]
	dataRow123 = data[111,:]
	dataRow124 = data[110,:]
	dataRow125 = data[109,:]
	rows12 = [dataRow122, dataRow123, dataRow124, dataRow125]

	
	
	#+++DATA ROW 13 (101-99)
	indRow13=100
	dataRow132 = data[101,:]
	dataRow133 = data[100,:]
	dataRow134 = data[99,:]
	dataRow135 = data[98,:]
	rows13 = [dataRow132, dataRow133, dataRow134, dataRow135]

	
	#+++DATA ROW 14 (87-83)
	indRow14=85
	dataRow141 = data[87,:]
	dataRow142 = data[86,:]
	dataRow143 = data[85,:]
	dataRow144 = data[84,:]
	dataRow145 = data[83,:]
	rows14 = [dataRow141, dataRow142, dataRow143, dataRow144, dataRow145]




	allRows = [rows1, rows2, rows3, rows4, rows5, rows6, rows7, rows8, rows9, rows10, rows11, rows12, rows13, rows14]


	########################################
	#finally do errores
	########################################




	sigma1 = []
	for i in range(0, len(rows1)):
		errors1 = DM_errors.DM_errors(rows1[i])
		sigma1.append(errors1)
	#print "sigma1", sigma1
	
	sigma2 = []
	for i in range(0, len(rows2)):
		errors2 = DM_errors.DM_errors(rows2[i])
		sigma2.append(errors2)
	#print "sigma2", sigma2

	sigma3 = []
	for i in range(0, len(rows3)):
		errors3 = DM_errors.DM_errors(rows3[i])
		sigma3.append(errors3)
	#print "sigma3", sigma3
	
	sigma4 = []
	for i in range(0, len(rows4)):
		errors4 = DM_errors.DM_errors(rows4[i])
		sigma4.append(errors4)
	#print "sigma4", sigma4
	
	sigma5 = []
	for i in range(0, len(rows5)):
		errors5 = DM_errors.DM_errors(rows5[i])
		sigma5.append(errors5)
	#print "sigma5", sigma5
	
	sigma6 = []
	for i in range(0, len(rows6)):
		errors6 = DM_errors.DM_errors(rows6[i])
		sigma6.append(errors6)
	#print "sigma6", sigma6
	
	sigma7 = []
	for i in range(0, len(rows7)):
		errors7 = DM_errors.DM_errors(rows7[i])
		sigma7.append(errors7)
	#print "sigma7", sigma7
	
	sigma8 = []
	for i in range(0, len(rows8)):
		errors8 = DM_errors.DM_errors(rows8[i])
		sigma8.append(errors8)
	#print "sigma8", sigma8
	
	sigma9 = []
	for i in range(0, len(rows9)):
		errors9 = DM_errors.DM_errors(rows9[i])
		sigma9.append(errors9)
	#print "sigma9", sigma9
	

	sigma10 = []
	for i in range(0, len(rows10)):
		errors10 = DM_errors.DM_errors(rows10[i])
		sigma10.append(errors10)
	#print "sigma10", sigma10
	
	sigma11 = []
	for i in range(0, len(rows11)):
		errors11 = DM_errors.DM_errors(rows11[i])
		sigma11.append(errors11)
	#print "sigma11", sigma11
	
	sigma12 = []
	for i in range(0, len(rows12)):
		errors12 = DM_errors.DM_errors(rows12[i])
		sigma12.append(errors12)
	#print "sigma12", sigma12
	

	sigma13 = []
	for i in range(0, len(rows13)):
		errors13 = DM_errors.DM_errors(rows13[i])
		sigma13.append(errors13)
	#print "sigma13", sigma13
	
	sigma14 = []
	for i in range(0, len(rows14)):
		errors14 = DM_errors.DM_errors(rows14[i])
		sigma14.append(errors14)
	#print "sigma14", sigma14
	
	##
	sigmaSubtract = []
	for i in range(0, len(dataRow)):
		errorsSkySub = DM_errors.DM_errors(dataRow[i])
		sigmaSubtract.append(errorsSkySub)
	#print "sigmaSubtract", sigmaSubtract
	
	medError = np.median(sigmaSubtract,0)
	
	medError.tolist()



	errorsSqrt = [sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, sigma8, sigma9, sigma10, sigma11, sigma12, sigma13, sigma14]
	
	#subtractedError = []
	#subtractedErrorRun = DM_errors.subtract_errors(sigma1,medError)
	#subtractedError.append(subtractedErrorRun)
	
	
	
	
	for i in range(0, len(sigma1)):
		finalE1=sigma1[i]-medError
		
		
	for i in range(0, len(sigma2)):
		finalE2=sigma2[i]-medError
		
	for i in range(0, len(sigma3)):
		finalE3=sigma3[i]-medError
		
	for i in range(0, len(sigma4)):
		finalE4=sigma4[i]-medError
		
	for i in range(0, len(sigma5)):
		finalE5=sigma5[i]-medError
		
	for i in range(0, len(sigma6)):
		finalE6=sigma6[i]-medError
		
	for i in range(0, len(sigma7)):
		finalE7=sigma7[i]-medError
		
	for i in range(0, len(sigma8)):
		finalE8=sigma8[i]-medError
		
	for i in range(0, len(sigma9)):
		finalE9=sigma9[i]-medError
		
	for i in range(0, len(sigma10)):
		finalE10=sigma10[i]-medError
		
	for i in range(0, len(sigma11)):
		finalE11=sigma11[i]-medError
		
	for i in range(0, len(sigma12)):
		finalE12=sigma12[i]-medError
		
	for i in range(0, len(sigma13)):
		finalE13=sigma13[i]-medError
		
	for i in range(0, len(sigma14)):
		finalE14=sigma14[i]-medError
		

	
	subtractedError=[finalE1, finalE2, finalE3, finalE4, finalE5, finalE6, finalE7, finalE8, finalE9, finalE10, finalE11, finalE12, finalE13, finalE14]
	
	
	error_rows1=DM_errors.add_errors(finalE1)
	error_rows2=DM_errors.add_errors(finalE2)
	error_rows3=DM_errors.add_errors(finalE3)
	error_rows4=DM_errors.add_errors(finalE4)
	error_rows5=DM_errors.add_errors(finalE5)
	error_rows6=DM_errors.add_errors(finalE6)
	error_rows7=DM_errors.add_errors(finalE7)
	error_rows8=DM_errors.add_errors(finalE8)
	error_rows9=DM_errors.add_errors(finalE9)
	error_rows10=DM_errors.add_errors(finalE10)
	error_rows11=DM_errors.add_errors(finalE11)
	error_rows12=DM_errors.add_errors(finalE12)
	error_rows13=DM_errors.add_errors(finalE13)
	error_rows14=DM_errors.add_errors(finalE14)

	
	finalError = [error_rows1, error_rows2, error_rows3, error_rows4, error_rows5, error_rows6, error_rows7, error_rows8, error_rows9, error_rows10, error_rows11, error_rows12, error_rows13, error_rows14]

	
	
	
	
	return finalError