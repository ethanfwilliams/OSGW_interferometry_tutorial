'''
02_calc_ncfs.py

Author: Ethan Williams (efwillia@caltech.edu)

- Load raw data
- FK filter data to isolate positive and negative velocity OSGW
- Whiten traces
- Calculate cross-correlations (for certain subarrays only)

'''

import os
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.io as spio
from obspy.signal.filter import bandpass as bp
from scipy.ndimage import gaussian_filter as gf

def detrend(y):
    ''' Remove a linear trend '''
    x = np.arange(len(y))
    m = (np.sum(x*y)-np.sum(x)*np.sum(y)/len(x))/(np.sum(x**2)-np.sum(x)**2/len(x))
    b = np.mean(y) - m*np.mean(x)
    return y - (m*x + b)

def whiten2(spc,fs,low,high):
    ''' Whiten the spectrum between low and high '''
    N = (len(spc)-1)*2
    i1 = int(math.ceil(low/(fs/float(N))))
    i2 = int(math.ceil(high/(fs/float(N))))
    spc[:i1] = np.cos(np.linspace(np.pi/2,np.pi,i1))**2 * np.exp(1j*np.angle(spc[:i1]))
    spc[i1:i2] = np.exp(1j*np.angle(spc[i1:i2]))
    spc[i2:] = np.cos(np.linspace(np.pi,np.pi/2,len(spc)-i2))**2 * np.exp(1j*np.angle(spc[i2:]))
    return spc

def fk_filt(tr_in,fs,dx,sgn='pos',cmin=5,cmax=50):
    ''' F-K filtering between cmin and cmax contours '''
    Nx = tr_in.shape[0]
    Ns = tr_in.shape[1]
    f0 = np.fft.fftshift(np.fft.fftfreq(Ns,d=1./fs))
    k0 = np.fft.fftshift(np.fft.fftfreq(Nx,d=dx))
    ft2 = np.fft.fftshift(np.fft.fft2(tr_in))
    F,K = np.meshgrid(f0,k0)
    C = F/K
    filt = np.zeros(ft2.shape)
    if sgn=='pos': # which direction to choose
        filt[np.logical_and(C>cmin,C<cmax)] = 1.
    else:
        filt[np.logical_and(C<-cmin,C>-cmax)] = 1.
    filt = gf(filt,3) # blur the filter a little to reduce Gibbs ringing
    ft2f = ft2*filt
    tr_out = np.fft.ifft2(np.fft.fftshift(ft2f)).astype(float)
    return tr_out

def get_src_gather_fk(data,src,rec_arr,nnx,nns,nnw,nwn,nov,fmin,fmax,fs,dx,sgn):
    ''' 
    Get a virtual source gather of cross-correlations from src to rec_arr 

    Same processing as above, but only for a single virtual source
    '''
    src_idx = np.argmin(np.abs(rec_arr-src))
    nc = len(rec_arr)
    trxc = np.zeros((nc,nns))
    spxc = np.zeros((nc,nnw),dtype=np.complex_)
    for n in range(nwn):
        tr = data[rec_arr,n*nov:(n*nov + nns)]
        for ic in range(nc):
            tr[ic,:] = detrend(tr[ic,:])
            tr[ic,:] = bp(tr[ic,:],df=fs,freqmin=fmin,freqmax=fmax,zerophase=True)
        tr = fk_filt(tr,fs,dx,sgn)
        sp = np.zeros((nc,nnw),dtype=np.complex_)
        for ic in range(nc):
            sp[ic,:] = np.fft.rfft(tr[ic,:])
            sp[ic,:] = whiten2(sp[ic,:],fs,fmin,fmax)
        for ic in range(nc):
            spxc[ic,:] += np.conj(sp[ic,:]) * sp[src_idx,:]
    for ic in range(nc):
        trxc[ic,:] = np.fft.irfft(spxc[ic,:])
    return trxc


#######################################################################################

## Set up paths
datafile = './data/Fs1Hz_dx10m_Belgium_data.npy'
outdir = './ncfs'

## Set parameters
fs = 1
dx = 10
fmin = 0.03 # min and max frequencies to whiten
fmax = 0.3
nns = int(2**9) # number of samples in a subwindow
nnw = int(nns/2+1) # number of samples in RFFT
nov = nns//2 # overlap between windows

## Set geometry
nx = 1000
ddx = 10 # spatial interval for sources
nr = 200 # offset to compute cross-correlations over
src_list = np.arange(0,nx-nr,ddx,dtype=int)

## Read data
all_data = np.load(datafile)

## Split data into two partially overlapping windows
ns = 2800 # number of samples in each window
datasets = [all_data[:,:ns],all_data[:,(all_data.shape[1]-ns):]]
nf = len(datasets)

## Set directions to filter
sgns = ['neg','pos']

## For each hour of data
for i in range(nf):
  print('Hour %d/%d' % (i+1,nf))
  ## Get data and dimensions
  data = datasets[i]
  nx = data.shape[0]
  ns = data.shape[1]
  nwn = int(np.floor(ns/nov))-1

  ## For each wavefield direction
  for sgn in sgns:
    print('\t calculating NCFs for '+sgn+' direction')
    ## Calculate cross-correlations for each souce
    for src in src_list:
      rec_arr = np.arange(src,src+nr)
      xc = get_src_gather_fk(data,src,rec_arr,nx,nns,nnw,nwn,nov,fmin,fmax,fs,dx,sgn)
      # Save
      np.save(os.path.join(outdir,sgn+'_ncf_src%04d_hr%01d.npy'%(src,i)),xc)



