import numpy as np
import random as r
import time 
import pdb

## main loop timed
start = time.time()

n 		= 27
NSTEP 	= 10001
vmax 	= 2.4
rcut 	= 3
L 		= 30  		# V = L ^ 3 	V  = L^3, L = 3*Sigma; we need D = N/V; 0.001 = (27)/(L)^3  
dt 		= 0.001

x = np.zeros([n,NSTEP]) # make array here
y = np.zeros([n,NSTEP]) # make array here
z = np.zeros([n,NSTEP]) # make array here

max_n = int(n**(1/3.0))

# for i in range(0,max_n):
# 	print isu
# 	z.append([i,i,i])
# j = np.array([item for sublist in z for item in sublist])
# print np.vstack((j,j))

c = 0

for i in range(0,max_n):
	for j in range(0,max_n):
		for k in range(0,max_n):
			x[c,0]= i
			y[c,0]= j
			z[c,0]= k
			c += 1

vx = np.zeros([n,NSTEP]) # make array here
vy = np.zeros([n,NSTEP]) # make array here
vz = np.zeros([n,NSTEP]) # make array here

for i in range(0,27):
	vx[i,0]=vmax*(2*r.random()-1)
	vy[i,0]=vmax*(2*r.random()-1)
	vz[i,0]=vmax*(2*r.random()-1)

# subtract centre of mass velocity
vxcm = vx[0:,0].mean()
vycm = vy[0:,0].mean()
vzcm = vz[0:,0].mean()

fx = np.zeros([n,NSTEP]) # make array here
fy = np.zeros([n,NSTEP]) # make array here
fz = np.zeros([n,NSTEP]) # make array here

# timing
test_start = time.time()

k = 0

# FORCE (k)
for i in range(0,n-1):
	for j in range(i+1,n):

		xmin = (x[i,k]-x[j,k])-L*round( (x[i,k]-x[j,k]) / L )
		ymin = (y[i,k]-y[j,k])-L*round( (y[i,k]-y[j,k]) / L )
		zmin = (z[i,k]-z[j,k])-L*round( (z[i,k]-z[j,k]) / L )

		rmin2 = xmin**2 + ymin**2 + zmin**2

		if rmin2 < rcut**2:
			f = (48/rmin2**7)-(24/rmin2**4)

			fx[i,k] += f*xmin
			fy[i,k] += f*ymin
			fz[i,k] += f*zmin
			fx[j,k] += -f*xmin
			fy[j,k] += -f*ymin
			fz[j,k] += -f*zmin

# MAIN 
for k in range(1,NSTEP):

	# MOVER (k-1)
	for i in range(0,n):
		x[i,k] = x[i,k-1] + dt*vx[i,k-1] + ((dt**2)/2)*fx[i,k-1]
		y[i,k] = y[i,k-1] + dt*vy[i,k-1] + ((dt**2)/2)*fy[i,k-1]
		z[i,k] = z[i,k-1] + dt*vz[i,k-1] + ((dt**2)/2)*fz[i,k-1]

	# FORCE (k)
	for i in range(0,n-1):
		for j in range(i+1,n):
			xmin = (x[i,k]-x[j,k])-L*round( (x[i,k]-x[j,k]) / L )
			ymin = (y[i,k]-y[j,k])-L*round( (y[i,k]-y[j,k]) / L )
			zmin = (z[i,k]-z[j,k])-L*round( (z[i,k]-z[j,k]) / L )

			rmin2 = xmin**2 + ymin**2 + zmin**2

			if rmin2 < rcut**2:
				f = (48/rmin2**7)-(24/rmin2**4)
				fx[i,k] += f*xmin
				fy[i,k] += f*ymin
				fz[i,k] += f*zmin
				fx[j,k] += -f*xmin
				fy[j,k] += -f*ymin
				fz[j,k] += -f*zmin

	# MOVEV (k-1)
	for i in range(0,n):
		vx[i,k]=vx[i,k-1]+(dt/2)*(fx[i,k]+fx[i,k-1])
		vy[i,k]=vy[i,k-1]+(dt/2)*(fy[i,k]+fy[i,k-1])
		vz[i,k]=vz[i,k-1]+(dt/2)*(fz[i,k]+fz[i,k-1])

test_end  = time.time()
test_time = test_end-test_start
#print "test: "+str(test_end-test_start)

# Mean square displacement -- diffusion coefficient estimate 

print "The Diffusion CO-efficient accurate to +/- 25 orders of magnitude"
#print ((x[0:,1000-1]-x[0:,NSTEP-1])**2+(y[0:,1000-1]-y[0:,NSTEP-1])**2+(z[0:,1000-1]-z[0:,NSTEP-1])**2).sum()/(6*dt*NSTEP*n)
print ((x[0:,9-1]-x[0:,NSTEP-1])**2+(y[0:,9-1]-y[0:,NSTEP-1])**2+(z[0:,9-1]-z[0:,NSTEP-1])**2).sum()/(6*dt*NSTEP*n)

end = time.time()

#print end - start

total_time = end - start

print "Proportional cost of test: "+str((test_time/total_time)*100)+"%"

