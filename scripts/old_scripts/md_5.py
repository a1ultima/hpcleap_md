import numpy as np
import math as m
import random as r
import time 
import pdb
import pickle

start = time.time() ## Timer Start Point

n 		= 27
NSTEP 	= 10000
vmax 	= 2.7
rcut 	= 3	
L 		= 30
dt 		= 0.001

print "Initializing coordinate array..."

data = { 	"x":np.zeros([n,NSTEP]),\
			"y":np.zeros([n,NSTEP]),\
			"z":np.zeros([n,NSTEP]),\
			"vx":np.zeros([n,NSTEP]),\
			"vy":np.zeros([n,NSTEP]),\
			"vz":np.zeros([n,NSTEP]),\
			"fx":np.zeros([n,NSTEP]),\
			"fy":np.zeros([n,NSTEP]),\
			"fz":np.zeros([n,NSTEP]),\
			"params":{	"n":n,\
						"NSTEP":NSTEP,\
						"vmax":vmax,\
						"rcut":rcut,\
						"L":L,\
						"dt":dt}}

c = 0

dist  = int(round(L/((n**(1/3.0)))))
max_n = int(round(n**(1/3.0)))*dist
for i in range(0,max_n,dist):
	for j in range(0,max_n,dist):
		for k in range(0,max_n,dist):
			data["x"][c,0]= i
			data["y"][c,0]= j
			data["z"][c,0]= k
			c += 1

print "Initializing random velocities..."

for i in range(0,n):
	data["vx"][i,0]=vmax*(2*r.random()-1)
	data["vy"][i,0]=vmax*(2*r.random()-1)
	data["vz"][i,0]=vmax*(2*r.random()-1)

# pdb.set_trace()

# subtract centre of mass velocity
vxcm = data["vx"][0:,0].mean()
vycm = data["vy"][0:,0].mean()
vzcm = data["vz"][0:,0].mean()

# data["vx"] = data["vx"]-vxcm
# data["vy"] = data["vy"]-vycm
# data["vz"] = data["vz"]-vzcm


print "Initializing forces..."

test_start = time.time()

k = 0

# FORCE (k)
for i in range(0,n-1):
	for j in range(i+1,n):

		xmin = (data["x"][i,k]-data["x"][j,k])-L*round( (data["x"][i,k]-data["x"][j,k]) / L )
		ymin = (data["y"][i,k]-data["y"][j,k])-L*round( (data["y"][i,k]-data["y"][j,k]) / L )
		zmin = (data["z"][i,k]-data["z"][j,k])-L*round( (data["z"][i,k]-data["z"][j,k]) / L )

		rmin2 = xmin**2 + ymin**2 + zmin**2

		if rmin2 < rcut**2:
			f = (48/(rmin2**7))-(24/(rmin2**4))

			data["fx"][i,k] += f*xmin
			data["fy"][i,k] += f*ymin
			data["fz"][i,k] += f*zmin

			data["fx"][j,k] += -f*xmin
			data["fy"][j,k] += -f*ymin
			data["fz"][j,k] += -f*zmin

print "Running main loop:..."

for k in range(1,NSTEP):

	if k%(NSTEP/100)==0:
		print "\t"+str(k/(NSTEP/100))+" percent complete..."

	# MOVER (k-1)
	for i in range(0,n):
		data["x"][i,k] = data["x"][i,k-1] + dt*data["vx"][i,k-1] + ((dt**2)/2)*data["fx"][i,k-1]
		data["y"][i,k] = data["y"][i,k-1] + dt*data["vy"][i,k-1] + ((dt**2)/2)*data["fy"][i,k-1]
		data["z"][i,k] = data["z"][i,k-1] + dt*data["vz"][i,k-1] + ((dt**2)/2)*data["fz"][i,k-1]

	# FORCE (k)
	for i in range(0,n-1):
		for j in range(i+1,n):
			xmin = (data["x"][i,k]-data["x"][j,k])-L*round( (data["x"][i,k]-data["x"][j,k]) / L )
			ymin = (data["y"][i,k]-data["y"][j,k])-L*round( (data["y"][i,k]-data["y"][j,k]) / L )
			zmin = (data["z"][i,k]-data["z"][j,k])-L*round( (data["z"][i,k]-data["z"][j,k]) / L )

			rmin2 = xmin**2 + ymin**2 + zmin**2

			if rmin2 < rcut**2:
				f = (48/rmin2**7)-(24/rmin2**4)
				data["fx"][i,k] += f*xmin
				data["fy"][i,k] += f*ymin
				data["fz"][i,k] += f*zmin
				data["fx"][j,k] += -f*xmin
				data["fy"][j,k] += -f*ymin
				data["fz"][j,k] += -f*zmin

	# MOVEV (k-1)
	for i in range(0,n):
		data["vx"][i,k]=data["vx"][i,k-1]+(dt/2)*(data["fx"][i,k]+data["fx"][i,k-1])
		data["vy"][i,k]=data["vy"][i,k-1]+(dt/2)*(data["fy"][i,k]+data["fy"][i,k-1])
		data["vz"][i,k]=data["vz"][i,k-1]+(dt/2)*(data["fz"][i,k]+data["fz"][i,k-1])

# Mean square displacement -- diffusion coefficient estimate 
DC=((data["x"][0:,1000-1]-data["x"][0:,NSTEP-1])**2+(data["y"][0:,1000-1]-data["y"][0:,NSTEP-1])**2+(data["z"][0:,1000-1]-data["z"][0:,NSTEP-1])**2).sum()/(6*dt*NSTEP*n)
print "The Diffusion Coefficient:", DC

end = time.time()
total_time = end - start

# save the data to file to be re-loaded next time
# fname = "cache_n"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+"_t"+str(time.time())+".p"
# pickle.dump(d,open( fname, "wb" ))


print "Test timer in seconds:", end - start    ##," ","Proportional cost of test: "+str((test_time/total_time)*100)+"%"
