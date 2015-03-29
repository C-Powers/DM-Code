import AngPerPixFunc
import skySubtractionFunc
import testGalaxyDataFunc
import galaxyDataFunc2
###NOTE THAT ERRORS MUST BE CALCULATED AS A SQUARE ROOT OF THE COUNTS "BEFORE" SKY SUBTRACTION

#====================  GET LAMBDA X VALUES ====================
calPix = [421.999, 429.988, 473.989, 484.963, 546.969, 653.0, 696.944]
calAng = [6382.99, 6402.25, 6506.53, 6532.88, 6678.2, 6929.47, 7032.41]

#546.966, 562.991, 867.989]
#6678.2, 6717.04, 7245.17]

blue = "b1009.fits"
red = "corrr1013.fits"

#Find wavelength as a function of pixels
Lambda_x = AngPerPixFunc.AngPerPix(red, calPix,calAng)

#Subtract Sky Lines
subtractedErrors = skySubtractionFunc.skySubtraction("corrdata.fits")

#run galaxy function
Galaxy = "skySubtract.fits"
galaxyDataFunc2.galaxyData(Galaxy, Lambda_x, subtractedErrors)