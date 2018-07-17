import numpy as np

dr=np.zeros((100))
r=np.zeros((100))

for n in range(100):
	if(n<9):dr[n]=1.
	elif(n>=9 and n<18):dr[n]=10
	else:dr[n]=100
	
r[0] = 1
for n in range(1,100):
	r[n] = r[n-1]+dr[n-1]