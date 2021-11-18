'''
04_calc_dispersion_example.py

Author: Ethan Williams (efwillia@caltech.edu)

For one example....
- Load and stack cross-correlation gathers
- Compute phase-shift beamforming to get dispersion image
- Pick dispersion curves

'''

import os
import math
import numpy as np
import matplotlib.pyplot as plt


def calcDispersion2(tr_,fs,dx,v_min,v_max,v_step,f_min,f_max):
    ''' calculate dispersion curve by phase-shift beamforming '''
    vel = np.arange(v_min,v_max,v_step)
    Nv = len(vel)
    Nx = tr_.shape[0]
    Ns = tr_.shape[1]
    sp_ = np.fft.rfft(tr_,axis=1)
    frq = np.fft.rfftfreq(Ns,d=1./fs)
    sp_shift = np.zeros((Nx,Ns//2+1,Nv),dtype=np.complex_)
    disp = np.zeros((Nv,Ns//2+1))
    for iv in range(Nv):
        for ix in range(tr_.shape[0]):
            sp_shift[ix,:,iv] = sp_[ix,:] * np.exp(2j*np.pi*frq*ix*dx/vel[iv])
        disp[iv,:] = np.mean(np.real(sp_shift[:,:,iv]),axis=0)
    return frq,vel,disp
   
#######################################################################################


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

## Parameters for dispersion image
vmin = 5 # velocity range
vmax = 20
vstep = 0.1 
fmin = 0.04 # frequency range
fmax = 0.2 
xmax = 400 # offsets to calculate dispersion (m)
tmax = 128 # time to calculate dispersion (s)
# the time also sets the frequency resolution

## Cut positive lags
idx = (x<xmax)
idt = np.logical_and(t>=0,t<=tmax)
trxc1 = ncf_com[np.ix_(idx,idt)]

## Cut negative lags
idx = (x<xmax)
idt = np.logical_and(t<=0,t>=-tmax)
trxc2 = ncf_com[np.ix_(idx,idt)]
trxc2 = trxc2[:,::-1] # flip to positive time axis

## Compute dispersion images
frq,vel,disp1 = calcDispersion2(trxc1,fs,dx,vmin,vmax,vstep,fmin,fmax)
frq,vel,disp2 = calcDispersion2(trxc2,fs,dx,vmin,vmax,vstep,fmin,fmax)

## Pick dispersion curves
N = len(frq)
pick1 = np.zeros(N)
pick2 = np.zeros(N)
vals1 = np.zeros(N)
vals2 = np.zeros(N)
for i in range(N):
  pick1[i] = vel[np.argmax(disp1[:,i])]
  pick2[i] = vel[np.argmax(disp2[:,i])]
  vals1[i] = max(disp1[:,i])
  vals2[i] = max(disp2[:,i])

## Remove picks below a threshold
threshold = 2 # manually set
bounds = np.logical_and(frq<fmax,frq>fmin) # inside frequency range of interest
# positive side
index1 = np.logical_and(vals1>threshold,bounds)
good_pick1 = pick1[index1]
good_freq1 = frq[index1]
bad_pick1 = pick1[~index1]
bad_freq1 = frq[~index1]
# negative side
index2 = np.logical_and(vals2>threshold,bounds)
good_pick2 = pick2[index2]
good_freq2 = frq[index2]
bad_pick2 = pick2[~index2]
bad_freq2 = frq[~index2]



## Plot
v = 1
fig,ax = plt.subplots(2,1,figsize=(6,10))
ax[0].pcolormesh(t[idt],x[idx],trxc1[:,::-1],vmin=-v,vmax=v,cmap='RdBu')
ax[0].pcolormesh(-t[idt],x[idx],trxc2[:,::-1],vmin=-v,vmax=v,cmap='RdBu')
ax[0].set_xlim([-60,60])
ax[0].set_xlabel('Time lag (s)')
ax[0].set_ylabel('Distance (m)')

ax[1].pcolormesh(frq,vel,disp1,cmap='jet')
ax[1].pcolormesh(-frq,vel,disp2,cmap='jet')
ax[1].scatter(good_freq1,good_pick1,c='k',label='accepted')
ax[1].scatter(-good_freq2,good_pick2,c='k')
ax[1].scatter(bad_freq1,bad_pick1,c='w',label='rejected')
ax[1].scatter(-bad_freq2,bad_pick2,c='w')
ax[1].set_xlim([-0.3,0.3])
ax[1].set_ylim([vmin,vmax])
ax[1].set_xlabel('Frequency (Hz)')
ax[1].set_ylabel('Phase speed (m/s)')
plt.legend()


plt.show()



