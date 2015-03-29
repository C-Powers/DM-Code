#error function
import numpy as np

def DM_errors(inpArr):

	errors=[]
	for i in range(0, len(inpArr)):
		if inpArr[i].all < 0:
			inpArr[i]==0
		if inpArr[i].all >= 0:
			sqrt = np.sqrt(inpArr[i])
		errors.append(sqrt)	
		
	for i in range(0, len(errors)):
		if np.isnan(errors[i]):
		#if errors[i] == "nan":
			errors[i]=0
			
			
	
	#print "INP ARRAY", inpArr
	print "ERRORS", errors
		
		
	return errors
	
	
	
def subtract_errors(array_to_sub, sub_array):

	final_error=[]
	for i in range(0, len(array_to_sub)):
		errorQuant = np.sqrt(np.power(array_to_sub[i],2) + np.power(sub_array[i], 2))
		final_error.append(errorQuant)
	return final_error
	
def add_errors(array1):
	final_array=[]
	for i in range(0, len(array1)):
		errorQuant = np.sqrt(np.sum(np.power(array1[i],2)))
		final_array.append(errorQuant)
	return final_array
	
	
def z_errors(lambdaError, lambdas):
	final_array=[]
	for i in range(0, len(lambdas)):
		sigma_z = np.sqrt( np.power((1/lambdas[10]),2) * np.power(lambdaError[10],2) + (1/4)*np.power(lambdas[10]/np.power(lambdas[i],2),2)*np.power(lambdaError[i],2))
		final_array.append(sigma_z)
	return final_array
	
	
def v_errors(inpArray):
	final_array=[]
	for i in range(0, len(inpArray)):
		c = 3*np.power(10,5) #km/s
		errorV = c* np.sqrt(np.power(inpArray,2))
		final_array.append(errorV)
	return final_array
