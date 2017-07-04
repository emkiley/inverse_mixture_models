import numpy as np
#from pylab import *

# Mixture Components-------------------------------------------------------
# Names of mixture components
n_mix1 = ['G', 'H', 'I', 'J', 'K', 'L']
# Properties of mixture components
epsR_mix1 = np.r_[62.4,64.0,64.8,66.0,68.0,70.0]
epsI_mix1 = np.r_[3.5,2.7,0.5,0.7,2.5,1.1]

# Target values-----------------------------------------------------------
target_names = ['Target']
target_epsR = [65.0]
target_epsI = [2.08]


# END OF USER INPUT --------------------------------------------------------

def get_fracs(e1,e2,e11,e21,e12,e22,e13,e23,alpha):
	# Raise inputs to alpha power
	e1=e1**(alpha)
	e2=e2**(alpha)
	e11=e11**(alpha)
	e21=e21**(alpha)
	e12=e12**(alpha)
	e22=e22**(alpha)
	e13=e13**(alpha)
	e23=e23**(alpha)
	
	# Compute volume fractions
	v1 = ((e22-e23)*(e1-e13) + (e13-e12)*(e2-e23))/((e11-e13)*(e22-e23)-(e21-e23)*(e12-e13))
	v2 = ((e23-e21)*(e1-e13) + (e11-e13)*(e2-e23))/((e11-e13)*(e22-e23)-(e21-e23)*(e12-e13))

#	Check the values against linear system solution
#	A = np.array([[e11-e13,e12-e13],[e21-e23,e22-e23]])
#	b = np.r_[e1-e13,e2-e23]
#	v = linalg.solve(A,b)

#	print v1-v[0]
#	print v2-v[1]
	
	return v1,v2,1-v1-v2


l = np.size(epsR_mix1)
m = np.size(target_epsR)

if np.size(target_epsI)!=m:
	print('Make sure the target data comes in *pairs*')

if np.size(epsI_mix1)!=l:
	print('Make sure the component data comes in *pairs*')

tfile = open("volumefractions.txt", "w+")

for n in range(0,m):
	tfile.write("\n\n\n\n\n-----------"+str(target_names[n])+"--------------------\n\n\n\n")
	
	tablecount = 0
	for i in range(0, l):
	    for j in range(i+1, l):
	        for k in range(j+1, l):
			v1,v2,v3 = get_fracs(target_epsR[n],target_epsI[n],epsR_mix1[i],epsI_mix1[i],epsR_mix1[j],epsI_mix1[j],epsR_mix1[k],epsI_mix1[k],1.0)
			if (v1>0 and (v2>0 and v3>0)):
				tablecount = tablecount +1
				printstring = target_names[n]+"   |  " + n_mix1[i] + "   |  " + n_mix1[j] + "    |  " +n_mix1[k] + "\n------------------------------------------------------------\n"
				tfile.write(printstring)
				printstring = "alpha = 1   | " + str(v1) + " | " + str(v2) + " | " + str(v3) + "\n"
				tfile.write(printstring)
	
				v1,v2,v3 = get_fracs(target_epsR[n],target_epsI[n],epsR_mix1[i],epsI_mix1[i],epsR_mix1[j],epsI_mix1[j],epsR_mix1[k],epsI_mix1[k],1.0/2.0)
				printstring = "alpha = 1/2 | " + str(v1) + " | " + str(v2) + " | " + str(v3) + "\n"
				tfile.write(printstring)
	
				v1,v2,v3 = get_fracs(target_epsR[n],target_epsI[n],epsR_mix1[i],epsI_mix1[i],epsR_mix1[j],epsI_mix1[j],epsR_mix1[k],epsI_mix1[k],1.0/3.0)
				printstring = "alpha = 1/3 | " + str(v1) + " | " + str(v2) + " | " + str(v3) + "\n"
				tfile.write(printstring)

				v1,v2,v3 = get_fracs(target_epsR[n],target_epsI[n],epsR_mix1[i],epsI_mix1[i],epsR_mix1[j],epsI_mix1[j],epsR_mix1[k],epsI_mix1[k],-1.0)
				printstring = "alpha = -1 | " + str(v1) + " | " + str(v2) + " | " + str(v3) + "\n" + "\n"
				tfile.write(printstring)

	tfile.write(target_names[n]+' had '+str(tablecount)+' many tables')

tfile.close()
