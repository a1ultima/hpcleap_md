import numpy as np
import math as m
import random as r
import time 
import pdb
import pickle

print "running multiple simulations..."

runTimes = []
diffusionCoefficients = []


for simulation in range(0,1):

	start = time.time() ## Timer Start Point

	# Parameters (raw)
	n 	= 27 
	NSTEP 	= 10001
	vmax 	= 2.7
	rcut 	= 3 
	L 	= 30
	dt 	= 0.001

	# Parameters

	# print outputs of f
	#fo = open("n"+str(n)+"_N"+str(NSTEP)+"_v"+str(vmax)+"_r"+str(rcut)+"_L"+str(L)+"_d"+str(dt)+".txt","w")
	#fo.close()

	#
	# Initialize lattice
	#

	# x = np.zeros([n,NSTEP])
	# y = np.zeros([n,NSTEP])
	# z = np.zeros([n,NSTEP])
	# x = [[0]*n]*NSTEP
	# y = [[0]*n]*NSTEP
	# z = [[0]*n]*NSTEP

	x = np.zeros([n])
	y = np.zeros([n])
	z = np.zeros([n])

	x1 = np.zeros([n])
	y1 = np.zeros([n])
	z1 = np.zeros([n])

	c = 0

	dist  = int(round(L/((n**(1/3.0)))))
	max_n = int(round(n**(1/3.0)))*dist

	for i in range(0,max_n,dist):
		for j in range(0,max_n,dist):
			for k in range(0,max_n,dist):
				x[c] = i
				y[c] = j
				z[c] = k
				c += 1
	# for i in range(0,max_n,dist):
	# 	for j in range(0,max_n,dist):
	# 		for k in range(0,max_n,dist):
	# 			x[c][0]= i
	# 			y[c][0]= j
	# 			z[c][0]= k
	# 			c += 1

	vx = np.zeros([n])
	vy = np.zeros([n])
	vz = np.zeros([n])

	vx1 = np.zeros([n])
	vy1 = np.zeros([n])
	vz1 = np.zeros([n])
	# vx = [[0]*n]*NSTEP
	# vy = [[0]*n]*NSTEP
	# vz = [[0]*n]*NSTEP

	for i in range(0,n):
		vx[i]=vmax*((2*r.random())-1)
		vy[i]=vmax*((2*r.random())-1)
		vz[i]=vmax*((2*r.random())-1)

	# for i in range(0,n):
	# 	vx[0:]=vmax*((2*0.2)-1)
	# 	vy[0:]=vmax*((2*0.5)-1)
	# 	vz[0:]=vmax*((2*0.7)-1)

	#
	# subtract centre of mass velocity
	#
	vxcm = vx[0:].mean()
	vycm = vy[0:].mean()
	vzcm = vz[0:].mean()

	# vx_0 = [vx[i,0] for i in vx]
	# vy_0 = [vy[i,0] for i in vy]
	# vz_0 = [vz[i,0] for i in vz]
	# mean_initial_vx = float(sum(initial_vx))/len(initial_vx)
	# mean_initial_vy = float(sum(initial_vy))/len(initial_vy)
	# mean_initial_vz = float(sum(initial_vz))/len(initial_vz)

	# vx[0:] -= vxcm	
	# vy[0:] -= vycm
	# vz[0:] -= vzcm

	fx = np.zeros([n])
	fy = np.zeros([n])
	fz = np.zeros([n])

	# t+1
	fx1 = np.zeros([n])
	fy1 = np.zeros([n])
	fz1 = np.zeros([n])

	test_start = time.time()

	k = 0

	store_x_1000 = 0
	store_y_1000 = 0
	store_z_1000 = 0

	store_x_NSTEP = 0
	store_y_NSTEP = 0
	store_z_NSTEP = 0

	# FORCE (k)
	for i in range(0,n-1):
		for j in range(i+1,n):

			xmin = (x[i]-x[j])-L*round( (x[i]-x[j]) / L )
			ymin = (y[i]-y[j])-L*round( (y[i]-y[j]) / L )
			zmin = (z[i]-z[j])-L*round( (z[i]-z[j]) / L )

			rmin2 = xmin**2 + ymin**2 + zmin**2

			if rmin2 < rcut**2:

				f = (48/(rmin2**7))-(24/(rmin2**4)) 

				"""f = (48/(rmin2**7))-(24/(rmin2**4))  // original """

				# fxmin = f*xmin
				# fymin = f*ymin
				# fzmin = f*zmin

				# fx[i,k] += fxmin
				# fy[i,k] += fymin
				# fz[i,k] += fzmin

				# fx[j,k] += -fxmin
				# fy[j,k] += -fymin
				# fz[j,k] += -fzmin

				fx[i] += f*xmin
				fy[i] += f*ymin
				fz[i] += f*zmin

				fx[j] += -f*xmin
				fy[j] += -f*ymin
				fz[j] += -f*zmin



	for k in range(1,NSTEP):

		if k%(NSTEP/100)==0:
			print "\t"+str(k/(NSTEP/100))+" percent complete..."

		# MOVER (k-1)
		for i in range(0,n):
			x1[i] = x[i] + dt*vx[i] + ((dt**2)/2)*fx[i]
			y1[i] = y[i] + dt*vy[i] + ((dt**2)/2)*fy[i]
			z1[i] = z[i] + dt*vz[i] + ((dt**2)/2)*fz[i]

		# FORCE (k)
		for i in range(0,n-1):
			for j in range(i+1,n):

				# OPTI: each x[i,k] needs only be defined once then accessed later

				xmin = (x1[i]-x1[j])-L*round((x1[i]-x1[j])/L)
				ymin = (y1[i]-y1[j])-L*round((y1[i]-y1[j])/L)
				zmin = (z1[i]-z1[j])-L*round((z1[i]-z1[j])/L)

				rmin2 = xmin**2 + ymin**2 + zmin**2

				if rmin2 < rcut**2:
					f = (48/rmin2**7)-(24/rmin2**4)

					# each each i-min * f value needs only be defined once

					# fxmin=fx[i,k]+f*xmin
					# fymin=fy[i,k]+f*ymin
					# fzmin=fz[i,k]+f*zmin

					# fx[i,k]=fxmin
					# fy[i,k]=fymin
					# fz[i,k]=fzmin
					# fx[j,k]=-fxmin
					# fy[j,k]=-fymin
					# fz[j,k]=-fzmin

					fx1[i] = fx[i] + f*xmin
					fy1[i] = fy[i] + f*ymin
					fz1[i] = fz[i] + f*zmin
					fx1[j] = fx[j] + -f*xmin
					fy1[j] = fy[j] + -f*ymin
					fz1[j] = fz[j] + -f*zmin

		# MOVEV (k-1)
		for i in range(0,n):
			vx1[i]=vx[i]+(dt/2)*(fx1[i]+fx[i])
			vy1[i]=vy[i]+(dt/2)*(fy1[i]+fy[i])
			vz1[i]=vz[i]+(dt/2)*(fz1[i]+fz[i])

			# Reset variables
			x[i]=x1[i]
			y[i]=y1[i]
			z[i]=z1[i]

			fx[i]=fx1[i]
			fy[i]=fy1[i]
			fz[i]=fz1[i]

			vx[i]=vx1[i]
			vy[i]=vy1[i]
			vz[i]=vz1[i]

		if k == 1000-1:
			store_x_1000 = x
			store_y_1000 = y
			store_z_1000 = z
		if k == NSTEP-1:
			store_x_NSTEP = x
			store_y_NSTEP = y
			store_z_NSTEP = z

	# Mean square displacement -- diffusion coefficient estimate 

	pdb.set_trace()

	#DC=((x[0:,1000-1]-x[0:,NSTEP-1])**2+(y[0:,1000-1]-y[0:,NSTEP-1])**2+(z[0:,1000-1]-z[0:,NSTEP-1])**2).sum()/(6*dt*NSTEP*n)
	DC=((store_x_1000-store_x_NSTEP)**2+(store_x_NSTEP-store_y_NSTEP)**2+(store_z_1000-store_z_NSTEP)**2).sum()/(6*dt*NSTEP*n)
	print "The Diffusion Coefficient:", DC

	end = time.time()
	total_time = end - start


	#d = {"x":x,"y":y,"z":z,"vx":vx,"vy":vy,"vz":vz,"fx":fx,"fy":fy,"fz":fz}
	#pickle.dump(d,open("n"+str(n)+"_N"+str(NSTEP)+"_v"+str(vmax)+"_r"+str(rcut)+"_L"+str(L)+"_d"+str(dt)+str(end)+".p","r"))

	#print "Test timer in seconds:", end - start    ##," ","Proportional cost of test: "+str((test_time/total_time)*100)+"%"

	runTimes.append(total_time)
	diffusionCoefficients.append(DC)

	runTimes_arr 				= np.array(runTimes)
	diffusionCoefficients_arr 	= np.array(diffusionCoefficients)

	meanRunTime = runTimes_arr.mean()
	meanDC 		= diffusionCoefficients_arr.mean()

print "multiple simulations complete!! "
print "     ...average runtime: "+str(meanRunTime)+"  "+str(runTimes)
#print "     ...average DC     : "+str(meanDC)+"  "+str(diffusionCoefficients)