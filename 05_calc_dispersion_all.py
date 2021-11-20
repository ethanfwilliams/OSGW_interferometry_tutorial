'''
05_calc_dispersion_all.py

Author: Ethan Williams (efwillia@caltech.edu)

For all subarrays...
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

## Set paths
ncfdir = './ncfs'
outdir = './disp'

## Set parameters
nx = 200          # number of channels in source gather
ns = int(2**9)    # number of samples in source gather
fs = 1            # sampling rate (Hz)
dx = 10           # channel spacing (m)

x = np.arange(nx)*dx # offsets
t = np.arange(-ns//2,ns//2)/fs # time lags

## Get array of sources (from step 02)
src_list = np.arange(0,1000-200,10,dtype=int)
nsrc = len(src_list)

## Parameters for dispersion image
vmin = 5 # velocity range
vmax = 20
vstep = 0.1
fmin = 0.04 # frequency range
fmax = 0.2
xmax = 400 # offsets to calculate dispersion (m)
tmax = 128 # time to calculate dispersion (s)
# the time also sets the frequency resolution
npck = int(tmax//2 + 1) # number of frequency bins where picks are made
threshold = 2 # manually set



## Iterate over each source and get dispersion picks
picks = np.zeros((nsrc,2*npck))
for j, src in enumerate(src_list):
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
    bounds = np.logical_and(frq<fmax,frq>fmin) # inside frequency range of interest
    index1 = np.logical_and(vals1>threshold,bounds)
    pick1[~index1] = np.nan
    index2 = np.logical_and(vals2>threshold,bounds)
    pick2[~index2] = np.nan

    ## Add to pre-allocated array
    picks[j,:] = np.concatenate((pick2[::-1],pick1),axis=0)

freqs = np.concatenate((-frq[::-1],frq),axis=0)

np.save(os.path.join(outdir,'dispersion_picks.npy'),picks)
np.save(os.path.join(outdir,'dispersion_freqs.npy'),freqs)
