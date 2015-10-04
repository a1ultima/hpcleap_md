import numpy as np
import math as m
import random as r
import time 
import pdb
import pickle

start = time.time() ## Timer Start Point

n 		= 27
NSTEP 	= 10001
vmax 	= 2.7
rcut 	= 3	
L 		= 30
dt 		= 0.001

print "Running multiple simulations... "

runTimes = []
diffusionCoefficients = []

chunkSize = 2000

for simulation in range(0,1):

	print "Simulation No. "+str(simulation)

	print "\tInitializing coordinate array..."

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

	# trying to grow the matrix to: chunkSize instead of NSTEP to save memory, and once an error is reached it would just shrink then grow the matrix 
	# data = { 	"x": np.zeros([n,chunkSize]),\
	# 			"y": np.zeros([n,chunkSize]),\
	# 			"z": np.zeros([n,chunkSize]),\
	# 			"vx":np.zeros([n,chunkSize]),\
	# 			"vy":np.zeros([n,chunkSize]),\
	# 			"vz":np.zeros([n,chunkSize]),\
	# 			"fx":np.zeros([n,chunkSize]),\
	# 			"fy":np.zeros([n,chunkSize]),\
	# 			"fz":np.zeros([n,chunkSize]),\
	# 			"params":{	"n":n,\
	# 						"NSTEP":NSTEP,\
	# 						"vmax":vmax,\
	# 						"rcut":rcut,\
	# 						"L":L,\
	# 						"dt":dt}}

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

	print "\tInitializing random velocities..."

	for i in range(0,n):
		data["vx"][i,0]=vmax*(2*r.random()-1)
		data["vy"][i,0]=vmax*(2*r.random()-1)
		data["vz"][i,0]=vmax*(2*r.random()-1)

	# pdb.set_trace()

	# subtract centre of mass velocity
	vxcm = data["vx"][0:,0].mean()
	vycm = data["vy"][0:,0].mean()
	vzcm = data["vz"][0:,0].mean()

	data["vx"] = data["vx"]-vxcm
	data["vy"] = data["vy"]-vycm
	data["vz"] = data["vz"]-vzcm


	print "\tInitializing forces..."

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

	print "\tRunning main loop:..."

	for k in range(1,NSTEP):

		if k%(NSTEP/100)==0:
			print "\t\t"+str(k/(NSTEP/100))+" percent complete..."

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
	print "\tThe Diffusion Coefficient:", DC

	end = time.time()
	total_time = end - start

	# save the data to file to be re-loaded next time
	# fname = "cache_n"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+"_t"+str(time.time())+".p"
	# pickle.dump(d,open( fname, "wb" ))
	file_x 	= open("b_data_x_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_y 	= open("b_data_y_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_z 	= open("b_data_z_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_vx = open("b_data_vx_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_vy = open("b_data_vy_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_vz = open("b_data_vz_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_fx = open("b_data_fx_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_fy = open("b_data_fy_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_fz = open("b_data_fz_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())

	# write to file 
	for row in data["x"]:
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_x.write(row)
	file_x.close()

	for row in data["y"]:
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_y.write(row)
	file_y.close()

	for row in data["z"]: 
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_z.write(row)
	file_z.close()

	for row in data["vx"]:
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_vx.write(row)
	file_vx.close()

	for row in data["vy"]:
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_vy.write(row)
	file_vy.close()

	for row in data["vz"]: 
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_vz.write(row)
	file_vz.close()

	for row in data["fx"]:
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_fx.write(row)
	file_fx.close()

	for row in data["fy"]:
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_fy.write(row)
	file_fy.close()

	for row in data["fz"]: 
		row = row.astype("|S")
		row = "\t".join(row)+"\n"
		file_fz.write(row)
	file_fz.close()



	print "\tTest timer in seconds:", end - start    ##," ","Proportional cost of test: "+str((test_time/total_time)*100)+"%"

	runTimes.append(total_time)
	diffusionCoefficients.append(DC)

	runTimes_arr 				= np.array(runTimes)
	diffusionCoefficients_arr 	= np.array(diffusionCoefficients)

	meanRunTime = runTimes_arr.mean()
	meanDC 		= diffusionCoefficients_arr.mean()

print "multiple simulations complete!! "
print "     ...average runtime: "+str(meanRunTime)+"  "+str(runTimes)