'''
07_visualize_mean.py

Author: Ethan Williams (efwillia@caltech.edu)

- Load dispersion picks and inversion results
- Plot dispersion picks and dispersion model from results
- Plot the results for (U,h) and compare with bathymetry

'''

import os 
import numpy as np
import matplotlib.pyplot as plt

## Flag to save plots or not
save_plots = True

## Set path
dispdir = './disp'

## Load picks, model, and solutions
H = np.load(os.path.join(dispdir,'H.npy'))
U = np.load(os.path.join(dispdir,'U.npy'))
picks = np.load(os.path.join(dispdir,'dispersion_picks.npy'))
model = np.load(os.path.join(dispdir,'dispersion_model.npy'))
freqs = np.load(os.path.join(dispdir,'dispersion_freqs.npy'))
mpts = H[:,0]

## Load bathymetry profile
bathy = np.genfromtxt('./data/bathymetry.txt')

## Plot dispersion picks and model
fig1,ax = plt.subplots(3,1,figsize=(8,10),sharex=True,sharey=True)
im1 = ax[0].pcolormesh(mpts,freqs,picks.T,cmap='jet',vmin=5,vmax=20)
ax[0].set_ylabel('Frequency (Hz)')
ax[0].set_title('Dispersion picks')
fig1.colorbar(im1,ax=ax[0])
im2 = ax[1].pcolormesh(mpts,freqs,model.T,cmap='jet',vmin=5,vmax=20)
ax[1].set_ylabel('Frequency (Hz)')
ax[1].set_title('Dispersion model')
fig1.colorbar(im2,ax=ax[1])
im3 = ax[2].pcolormesh(mpts,freqs,(picks-model).T,cmap='RdBu',vmin=-5,vmax=5)
ax[2].set_ylabel('Frequency (Hz)')
ax[2].set_xlabel('Subarray midpoint (m)')
ax[2].set_title('Residual')
ax[2].set_ylim([-0.2,0.2])
fig1.colorbar(im3,ax=ax[2])

## Plot inversion results with 1-sigma uncertainty
fig2,ax = plt.subplots(2,1,figsize=(8,7),sharex=True)
ax[0].fill_between(H[:,0],H[:,1]+H[:,2],H[:,1]-H[:,2],color='gray')
ax[0].plot(H[:,0],H[:,1],'k',label='inversion result')
ax[0].plot(bathy[:,0],bathy[:,1],'r',label='EMODNET bathymetry')
ax[0].set_ylabel('Water depth (m)')
ax[0].set_title('Inverted depth profile')
ax[0].invert_yaxis()
ax[0].legend()
ax[1].fill_between(U[:,0],U[:,1]+U[:,2],U[:,1]-U[:,2],color='gray')
ax[1].plot(U[:,0],U[:,1],'k')
ax[1].set_ylabel('Current velocity (m/s)')
ax[1].set_xlabel('Subarray midpoint (m)')
ax[1].set_title('Inverted current velocity')
ax[1].set_xlim([min(H[:,0]),max(H[:,0])])


## Save
if save_plots:
  fig1.savefig('./figs/07a_dispersion.png')
  fig2.savefig('./figs/07b_mean_state.png')

plt.show()


