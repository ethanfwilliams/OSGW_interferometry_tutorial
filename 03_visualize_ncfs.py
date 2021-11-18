'''
03_visualize_ncfs.py

Author: Ethan Williams (efwillia@caltech.edu)

- Load and plot cross-correlation gathers

'''

import os
import numpy as np
import matplotlib.pyplot as plt


## Flag to save plots or not
save_plots = True

## Set path and source gather to plot
ncfdir = './ncfs'
src = 500

## Set parameters
nx = 200          # number of channels in source gather
ns = int(2**9)    # number of samples in source gather
fs = 1            # sampling rate (Hz)
dx = 10           # channel spacing (m)

x = np.arange(nx)*dx # offsets
t = np.arange(-ns//2,ns//2)/fs # time lags

ncf_neg = np.zeros((nx,ns))
ncf_pos = np.zeros((nx,ns))

## Load and stack data
hrs = [0,1]
for hr in hrs:
  fpath1 = os.path.join(ncfdir,'neg_ncf_src%04d_hr%01d.npy' % (src,hr))
  fpath2 = os.path.join(ncfdir,'pos_ncf_src%04d_hr%01d.npy' % (src,hr))
  ncf_neg += np.load(fpath1)
  ncf_pos += np.load(fpath2)
ncf_neg /= len(hrs)
ncf_pos /= len(hrs)

## Reorganize time lags
ncf_neg = np.concatenate((ncf_neg[:,ns//2:],ncf_neg[:,:ns//2]),axis=1)
ncf_pos = np.concatenate((ncf_pos[:,ns//2:],ncf_pos[:,:ns//2]),axis=1)
ncf_com = np.concatenate((ncf_neg[:,:ns//2],ncf_pos[:,ns//2:]),axis=1)

## Plot data
v = 1
fig,ax = plt.subplots(1,3,figsize=(15,6),sharex=True,sharey=True)
ax[0].pcolormesh(x,t,ncf_neg.T,cmap='RdBu',vmin=-v,vmax=v)
ax[0].set_xlabel('Offset (m)')
ax[0].set_ylabel('Time lag (s)')
ax[0].set_title('Anti-causal (negative)')
ax[0].set_xlim([0,1000])
ax[0].set_ylim([-100,100])
ax[1].pcolormesh(x,t,ncf_pos.T,cmap='RdBu',vmin=-v,vmax=v)
ax[1].set_xlabel('Offset (m)')
ax[1].set_title('Causal (positive)')
ax[2].pcolormesh(x,t,ncf_com.T,cmap='RdBu',vmin=-v,vmax=v)
ax[2].set_xlabel('Offset (m)')
ax[2].set_title('Combined')

## Save
if save_plots:
  fig.savefig('./figs/03_ncf_stack.png')


plt.show()

