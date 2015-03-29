import AngPerPixFunc
import skySubtractionFunc



#====================  GET LAMBDA X VALUES ====================
calPix = [421.999, 429.988, 473.989, 546.966, 562.991, 867.989]
calAng = [6382.99, 6402.25, 6506.53, 6678.2, 6717.04, 7245.17]
blue = "b1009.fits"
red = "r1013.fits"

Lambda_x = AngPerPixFunc.AngPerPix(red, calPix,calAng)

skySubtract = skySubtractionFunc.skySubtraction("r1104.fits")


