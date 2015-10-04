import numpy as np
import random as r
import time 
import copy

start = time.time() ## Timer Start Point

n 		= 27
NSTEP 	= 1000001
vmax 	= 2.7
rcut 	= 3	
L 		= 30
dt 		= 0.001

print "Running multiple simulations... "

runTimes = []
diffusionCoefficients = []

for simulation in range(0,3):

	print "Simulation No. "+str(simulation)

	print "\tInitializing coordinate array..."

	# Trying to reduce memory: from: md_7_multiRun_dictMode.py
	data = { 	"x":np.zeros([n,2]),\
				"y":np.zeros([n,2]),\
				"z":np.zeros([n,2]),\
				"vx":np.zeros([n,2]),\
				"vy":np.zeros([n,2]),\
				"vz":np.zeros([n,2]),\
				"fx":np.zeros([n,2]),\
				"fy":np.zeros([n,2]),\
				"fz":np.zeros([n,2]),\
				"t":np.array(range(0,NSTEP))*dt,\
				"T":np.zeros([NSTEP]),\
				"params":{	"n":n,\
							"2":NSTEP,\
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

	# file object for writing into...
	file_x 	= open("../data/a_data_s"+str(simulation)+"_x_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_y 	= open("../data/a_data_s"+str(simulation)+"_y_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_z 	= open("../data/a_data_s"+str(simulation)+"_z_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_vx = open("../data/a_data_s"+str(simulation)+"_vx_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_vy = open("../data/a_data_s"+str(simulation)+"_vy_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_vz = open("../data/a_data_s"+str(simulation)+"_vz_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_fx = open("../data/a_data_s"+str(simulation)+"_fx_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_fy = open("../data/a_data_s"+str(simulation)+"_fy_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())
	file_fz = open("../data/a_data_s"+str(simulation)+"_fz_K"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+".txt","w") # "_t"+str(time.time())

	print "\tInitializing forces..."

	test_start = time.time()

	k = 0

	# FORCE (k)
	for i in range(0,n-1):
		for j in range(i+1,n):
			
			xmin = (copy.copy(data["x"][i,0])-copy.copy(data["x"][j,0]))-L*round(copy.copy((data["x"][i,0])-copy.copy(data["x"][j,0]))/L)
			ymin = (copy.copy(data["y"][i,0])-copy.copy(data["y"][j,0]))-L*round(copy.copy((data["y"][i,0])-copy.copy(data["y"][j,0]))/L)
			zmin = (copy.copy(data["z"][i,0])-copy.copy(data["z"][j,0]))-L*round(copy.copy((data["z"][i,0])-copy.copy(data["z"][j,0]))/L)

			rmin2 = xmin**2 + ymin**2 + zmin**2

			if rmin2 < rcut**2:
				f = (48/rmin2**7)-(24/rmin2**4)
				data["fx"][i,0] = f*xmin
				data["fy"][i,0] = f*ymin
				data["fz"][i,0] = f*zmin
				data["fx"][j,0] = -f*xmin
				data["fy"][j,0] = -f*ymin
				data["fz"][j,0] = -f*zmin

	print "\tRunning main loop:..."

	for k in range(1,NSTEP+1):

		if k%(NSTEP/100)==0:
			print "\t"+str(k/(NSTEP/100))+" percent complete..."

		if k==1000-1:
			print "\t\tk: "+str(k)
			# Coordinates: 	UPDATE the time steps, setting the previous time-step's values as the current time-step
			x_start=copy.copy(data["x"][0:,0])
			y_start=copy.copy(data["y"][0:,0])
			z_start=copy.copy(data["z"][0:,0])

			# Forces: 		UPDATE the time steps, setting the previous time-step's values as the current time-step
			fx_start=copy.copy(data["fx"][0:,0])
			fy_start=copy.copy(data["fy"][0:,0])
			fz_start=copy.copy(data["fz"][0:,0])

			# Velocites: 	UPDATE the time steps, setting the previous time-step's values as the current time-step
			vx_start=copy.copy(data["vx"][0:,0])
			vy_start=copy.copy(data["vy"][0:,0])
			vz_start=copy.copy(data["vz"][0:,0])	

		if k==NSTEP-1:
			print "\t\tk: "+str(k)
			# Coordinates: 	UPDATE the time steps, setting the previous time-step's values as the current time-step
			x_end=copy.copy(data["x"][0:,0])
			y_end=copy.copy(data["y"][0:,0])
			z_end=copy.copy(data["z"][0:,0])

			# Forces: 		UPDATE the time steps, setting the previous time-step's values as the current time-step
			fx_end=copy.copy(data["fx"][0:,0])
			fy_end=copy.copy(data["fy"][0:,0])
			fz_end=copy.copy(data["fz"][0:,0])

			# Velocites: 	UPDATE the time steps, setting the previous time-step's values as the current time-step
			vx_end=copy.copy(data["vx"][0:,0])
			vy_end=copy.copy(data["vy"][0:,0])
			vz_end=copy.copy(data["vz"][0:,0])	

		# @TODO: collect string'd data into arrays/lists that build up to len() == 5000, and only then write data into the file outputs
		row = data["x"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_x.write(row)
	
		row = data["y"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_y.write(row)

		row = data["z"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_z.write(row)
		
		row = data["vx"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_vx.write(row)
		
		row = data["vy"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_vy.write(row)
		
		row = data["vz"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_vz.write(row)
		
		row = data["fx"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_fx.write(row)
		
		row = data["fy"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_fy.write(row)
		
		row = data["fz"][0:,0].T.astype("|S")
		row = "\t".join(row)+"\n"
		file_fz.write(row)

		# MOVER (k-1)
		for i in range(0,n):
			# 	Velocity Verlet: 
			# New coordinate 	old coordinate 		velocity*time   			acceleration?
			# data["x"][i,k] = data["x"][i,k-1] + dt*data["vx"][i,k-1] + ((dt**2)/2)*data["fx"][i,k-1]
			# data["y"][i,k] = data["y"][i,k-1] + dt*data["vy"][i,k-1] + ((dt**2)/2)*data["fy"][i,k-1]
			# data["z"][i,k] = data["z"][i,k-1] + dt*data["vz"][i,k-1] + ((dt**2)/2)*data["fz"][i,k-1]
			data["x"][i,1] = copy.copy(data["x"][i,0]) + dt*copy.copy(data["vx"][i,0]) + ((dt**2)/2)*copy.copy(data["fx"][i,0])
			data["y"][i,1] = copy.copy(data["y"][i,0]) + dt*copy.copy(data["vy"][i,0]) + ((dt**2)/2)*copy.copy(data["fy"][i,0])
			data["z"][i,1] = copy.copy(data["z"][i,0]) + dt*copy.copy(data["vz"][i,0]) + ((dt**2)/2)*copy.copy(data["fz"][i,0])

		# FORCE (k)
		for i in range(0,n-1):
			for j in range(i+1,n):
				
				xmin = (copy.copy(data["x"][i,1])-copy.copy(data["x"][j,1]))-L*round(copy.copy((data["x"][i,1])-copy.copy(data["x"][j,1]))/L)
				ymin = (copy.copy(data["y"][i,1])-copy.copy(data["y"][j,1]))-L*round(copy.copy((data["y"][i,1])-copy.copy(data["y"][j,1]))/L)
				zmin = (copy.copy(data["z"][i,1])-copy.copy(data["z"][j,1]))-L*round(copy.copy((data["z"][i,1])-copy.copy(data["z"][j,1]))/L)

				rmin2 = xmin**2 + ymin**2 + zmin**2

				if rmin2 < rcut**2:
					f = (48/rmin2**7)-(24/rmin2**4)
					data["fx"][i,1] = f*xmin
					data["fy"][i,1] = f*ymin
					data["fz"][i,1] = f*zmin
					data["fx"][j,1] = -f*xmin
					data["fy"][j,1] = -f*ymin
					data["fz"][j,1] = -f*zmin

		# MOVEV (k-1)
		for i in range(0,n):
			data["vx"][i,1]=copy.copy(data["vx"][i,0])+(dt/2)*(copy.copy(data["fx"][i,1])+copy.copy(data["fx"][i,0]))
			data["vy"][i,1]=copy.copy(data["vy"][i,0])+(dt/2)*(copy.copy(data["fy"][i,1])+copy.copy(data["fy"][i,0]))
			data["vz"][i,1]=copy.copy(data["vz"][i,0])+(dt/2)*(copy.copy(data["fz"][i,1])+copy.copy(data["fz"][i,0]))

		# UPDATE the arrays
		for z in range(0,n):
			# Coordinates: 	UPDATE the time steps, setting the previous time-step's values as the current time-step
			data["x"][z,0]=copy.copy(data["x"][z,1])
			data["y"][z,0]=copy.copy(data["y"][z,1])
			data["z"][z,0]=copy.copy(data["z"][z,1])

			# Forces: 		UPDATE the time steps, setting the previous time-step's values as the current time-step
			data["fx"][z,0]=copy.copy(data["fx"][z,1])
			data["fy"][z,0]=copy.copy(data["fy"][z,1])
			data["fz"][z,0]=copy.copy(data["fz"][z,1])

			data["fx"][z,1]=0
			data["fy"][z,1]=0
			data["fz"][z,1]=0

			# Velocities: 	UPDATE the time steps, setting the previous time-step's values as the current time-step
			data["vx"][z,0]=copy.copy(data["vx"][z,1])
			data["vy"][z,0]=copy.copy(data["vy"][z,1])
			data["vz"][z,0]=copy.copy(data["vz"][z,1])

	# Mean square displacement -- diffusion coefficient estimate 
	#DC=((data["x"][0:,1000-1]-data["x"][0:,NSTEP-1])**2+(data["y"][0:,1000-1]-data["y"][0:,NSTEP-1])**2+(data["z"][0:,1000-1]-data["z"][0:,NSTEP-1])**2).sum()/(6*dt*NSTEP*n)
	#pdb.set_trace()

	DC=((x_start-x_end)**2+(y_start-y_end)**2+(z_start-z_end)**2).sum()/(6*dt*NSTEP*n)
	print "\tThe Diffusion Coefficient:", DC

	# Close streaming files
	file_x.close()
	file_y.close()
	file_z.close()
	file_vx.close()
	file_vy.close()
	file_vz.close()
	file_fx.close()
	file_fy.close()
	file_fz.close()

	end = time.time()
	total_time = end - start

	# save the data to file to be re-loaded next time
	# fname = "cache_n"+str(NSTEP)+"_dt"+str(dt)+"_L"+str(L)+"_t"+str(time.time())+".p"
	# pickle.dump(d,open( fname, "wb" ))

	print "\tTest timer in seconds:", end - start    ##," ","Proportional cost of test: "+str((test_time/total_time)*100)+"%"

	runTimes.append(total_time)
	diffusionCoefficients.append(DC)

	runTimes_arr 				= np.array(runTimes)
	diffusionCoefficients_arr 	= np.array(diffusionCoefficients)

	meanRunTime = runTimes_arr.mean()
	meanDC 		= diffusionCoefficients_arr.mean()

	print "multiple simulations complete!! "
	print "     ...mean runtime: "+str(meanRunTime)+"  "+str(runTimes)
	print "     ...mean DC: "+str(meanDC)+"  "+str(diffusionCoefficients)
	