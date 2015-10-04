import numpy as np
import math as m
import random as r
import time 
import pdb
import pickle

start = time.time() ## Timer Start Point

n 		= 27
NSTEP 	= 500001
vmax 	= 2.7
rcut 	= 3	
L 		= 30
dt 		= 0.001

print "Initializing coordinate array..."

x = np.zeros([n,NSTEP]) # make array here
y = np.zeros([n,NSTEP]) # make array here
z = np.zeros([n,NSTEP]) # make array here

c = 0

dist  = int(round(L/((n**(1/3.0)))))
max_n = int(round(n**(1/3.0)))*dist
for i in range(0,max_n,dist):
	for j in range(0,max_n,dist):
		for k in range(0,max_n,dist):
			x[c,0]= i
			y[c,0]= j
			z[c,0]= k
			c += 1

print "Initializing random velocities..."

# d = { 	"x":x,
# 		"y":y,\
# 		"z":z,\
# 		"vx":np.zeros([n,NSTEP]),\
# 		"vy":np.zeros([n,NSTEP]),\
# 		"vz":np.zeros([n,NSTEP]),\
# 		"fx":np.zeros([n,NSTEP]),\
# 		"fy":np.zeros([n,NSTEP]),\
# 		"fz":np.zeros([n,NSTEP]),\
# 		"params":{"n":n,"NSTEP":NSTEP,"vmax":vmax,"rcut":rcut,"L":L,"dt":dt}}

for i in range(0,n):
	vx[i,0]=vmax*(2*r.random()-1)
	vy[i,0]=vmax*(2*r.random()-1)
	vz[i,0]=vmax*(2*r.random()-1)

# subtract centre of mass velocity
vxcm = vx[0:,0].mean()
vycm = vy[0:,0].mean()
vzcm = vz[0:,0].mean()

print "Initializing forces..."

fx = np.zeros([n,NSTEP])
fy = np.zeros([n,NSTEP])
fz = np.zeros([n,NSTEP])

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
			f = (48/(rmin2**7))-(24/(rmin2**4))

			fx[i,k] += f*xmin
			fy[i,k] += f*ymin
			fz[i,k] += f*zmin

			fx[j,k] += -f*xmin
			fy[j,k] += -f*ymin
			fz[j,k] += -f*zmin

print "Running main loop:..."

for k in range(1,NSTEP):

	if k%(NSTEP/100)==0:
		print "\t"+str(k/(NSTEP/100))+" percent complete..."

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

# Mean square displacement -- diffusion coefficient estimate 
DC=((x[0:,1000-1]-x[0:,NSTEP-1])**2+(y[0:,1000-1]-y[0:,NSTEP-1])**2+(z[0:,1000-1]-z[0:,NSTEP-1])**2).sum()/(6*dt*NSTEP*n)
print "The Diffusion Coefficient:", DC

end = time.time()
total_time = end - start

# save the data to file to be re-loaded next time
d = {"x":x,"y":y,"z":z,"vx":vx,"vy":vy,"vz":vz,"fx":fx,"fy":fy,"fz":fz,
		"params":{"n":n,"NSTEP":NSTEP,"vmax":vmax,"rcut":rcut,"L":L,"dt":dt}}
fname = "cache_n"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+"_t"+str(time.time())+".p"
pickle.dump(d,open( fname, "wb" ))


print "Test timer in seconds:", end - start    ##," ","Proportional cost of test: "+str((test_time/total_time)*100)+"%"
