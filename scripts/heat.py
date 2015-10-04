
# Dependencies: sudo apt-get install python-matplotlib

# Instantaneous temperature vs. time
import pickle
#import pdb
import numpy as np

import matplotlib.pyplot as plt 

f=open("./cache_n100001_t1443608397.21.p","r")

pick = pickle.load(f)

L 		= pick['params']['L']
rcut 	= pick['params']['rcut']
n 		= pick['params']['n']
vmax 	= pick['params']['vmax']
dt 		= pick['params']['dt']
NSTEP 	= pick['params']['NSTEP']

vx = pick['vx']
vy = pick['vy']
vz = pick['vz']

f.close()

# temp_vs_time = np.array([])
times = np.array(range(0,NSTEP))*dt
heats = []

for i,k in enumerate(times):
	sumv = (vx[0:,i]**2 + vy[0:,i]**2 + vz[0:,i]**2).sum()
	heat = sumv/(3*n)
	heats.append(heat)

# red dashes, blue squares and green triangles
plt.plot(times, heats, 'r')

fig = plt.gcf()

#plot_url = py.plot_mpl(fig, filename='mpl-line-style')

plt.show()
