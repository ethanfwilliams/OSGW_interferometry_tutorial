'''
09_visualize_results.py

Author: Ethan Williams (efwillia@caltech.edu)

- Load the differential and mean current, sum, and plot them together 

'''

import numpy as np
import matplotlib.pyplot as plt


## Flag to save plots or not
save_plots = True

nt = 2 # number of hours

## Load mean current data
tmp = np.load('./disp/U.npy')
x_m = tmp[:,0]
U_m = tmp[:,1]
H = np.load('./disp/U.npy')[:,1]

## Load differential current data
x_d = np.load('./disp/midpts.npy')
U_d = np.load('./disp/U_diff.npy')

## Interpolate onto the same grid
N = 100
x = np.linspace(0,10000,N)
Um = np.interp(x,x_m,U_m)
U = np.zeros((nt,N))
for i in range(nt):
  U[i,:] = np.interp(x,x_d,U_d[:,i]) + Um

## Plot
fig,ax = plt.subplots(1,1,figsize=(8,5))
ax.plot(x,U[0,:],'r',label='hour 0')
ax.plot(x,U[1,:],'b',label='hour 1')
ax.set_xlim([0,10000])
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Current velocity (m/s)')
ax.legend()
ax.set_title('Absolute current')

## Save
if save_plots:
  fig.savefig('./figs/09_abs_current.png')



plt.show()



