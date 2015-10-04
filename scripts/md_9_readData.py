
import numpy as np
import pdb

fi_x = open("../data/a_data_x_K250001_dt0.001_L30.txt","r")
fi_y = open("../data/a_data_y_K250001_dt0.001_L30.txt","r")
fi_z = open("../data/a_data_z_K250001_dt0.001_L30.txt","r")

#DC=((x_start-x_end)**2+(y_start-y_end)**2+(z_start-z_end)**2).sum()/(6*dt*NSTEP*n)

k = 0

x_start = 0 
x_end 	= 0

y_start = 0
y_end  	= 0

z_start = 0
z_end  	= 0

# Parameters
n 		= 27
NSTEP 	= 250000
vmax 	= 2.7
rcut 	= 3	
L 		= 30
dt 		= 0.001

# Main
while True:

	line_x = fi_x.readline()
	line_y = fi_y.readline()
	line_z = fi_z.readline()

	if line_x=="":
		print "ended line x"
		if line_y=="":
			print "ended line y"
			if line_z=="":
				print "ended line z"
				break
	if "#" in line_x:
		continue

	k += 1

	if k == 1000:
		x_start = np.array(line_x.rstrip().split("\t")).astype("|f")
		y_start = np.array(line_y.rstrip().split("\t")).astype("|f")
		z_start = np.array(line_z.rstrip().split("\t")).astype("|f")
	elif k == NSTEP:
		x_end 	= np.array(line_x.rstrip().split("\t")).astype("|f")
		y_end 	= np.array(line_y.rstrip().split("\t")).astype("|f")
		z_end 	= np.array(line_z.rstrip().split("\t")).astype("|f")

# Calculate Diffusion Coefficient 
DC=((x_start-x_end)**2+(y_start-y_end)**2+(z_start-z_end)**2).sum()/(6*dt*NSTEP*n)

print DC

fi_x.close()
fi_y.close()
fi_z.close()




