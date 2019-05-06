import numpy as np
import csv

file = "COLLV2.dat"
with open(file)as fn:
	content = fn.readlines()

ni = np.zeros((5))
an = np.zeros((4))
cn = np.zeros((4))
nu = np.zeros((8))
rho = np.zeros((9))	
coll = np.zeros((5,4,4,8,9))
ncoll = np.zeros((5,4,4,8,9))

count = 0
for i in range(5):
	for j in range(4):
		for k in range(4):
			for l in range(8):
				for m in range(9):
				
					x = content[count].rstrip().split()
					#y = content[count+1].rstrip().split()
					
					ni[i]=float(x[0])
					an[j]=float(x[1])
					cn[k]=float(x[2])
					nu[l]=float(x[3])
					rho[m]=float(x[4])
					coll[i,j,k,l,m]=float(x[5])
					ncoll[i,j,k,l,m]=float(x[6])
				
					#rho[m]=float(y[0])
					#coll[i,j,k,l,m]=float(y[1])
					#ncoll[i,j,k,l,m]=float(y[2])
					
					count = count + 1
					#count = count + 2



