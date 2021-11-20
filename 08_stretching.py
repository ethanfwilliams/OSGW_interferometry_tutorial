'''
08_stretching.py

Author: Ethan Williams (efwillia@caltech.edu)

- Load a time-lapse gather for a each pair
- Calculate dictionary of stretched traces
- Grid search to find the best fit

'''

import os
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.filter import bandpass as bp
import scipy.io as spio

class TL_gather:
  ncfdir = './ncfs'
  ns = int(2**9)
  fs = 1
  dx = 10
  dt = 2800/3600
  g = 9.8

  def __init__(self,src,off):
    self.src = src # source channel number
    self.off = off # offset in channels
    self.dr = self.dx*self.off # offset in meters
    self.hrs = [0,1] # hours where NCFs are available
    self.nt = len(self.hrs)

  def load_gather(self):
    ''' Load a time-lapse gather for both pos/neg sides and filter it '''
    # Load data
    srcdir = './ncfs' 
    xc1 = np.zeros((self.nt,self.ns))
    xc2 = np.zeros((self.nt,self.ns))
    for it, hr in enumerate(self.hrs): 
      fpath1 = os.path.join(self.ncfdir,'pos_ncf_src%04d_hr%01d.npy' % (self.src,hr))
      xc1[it,:] = np.load(fpath1)[off,:] # get only channel pair with given offset
      fpath2 = os.path.join(self.ncfdir,'neg_ncf_src%04d_hr%01d.npy' % (self.src,hr))
      xc2[it,:] = np.load(fpath2)[off,:]
    # Assign values
    self.pos = xc2 # I accidentally swapped pos/neg in the code... will fix at some point
    self.neg = xc1
    self.lags = np.arange(-self.ns//2,self.ns//2) / self.fs
    self.time = np.arange(self.nt) * self.dt
    # Filter
    for it in range(self.nt):
      self.pos[it,:] = bp(self.pos[it,:],freqmin=0.03,freqmax=0.2,df=self.fs,zerophase=True)
      self.neg[it,:] = bp(self.neg[it,:],freqmin=0.03,freqmax=0.2,df=self.fs,zerophase=True)

  def load_mean(self):
    ''' Load the mean state parameters from step 06 '''
    x_mean = np.load('./disp/H.npy')[:,0]
    h_mean = np.load('./disp/H.npy')[:,1]
    u_mean = np.load('./disp/U.npy')[:,1]
    mpt = np.mean((self.src,self.src+self.off))*self.dx
    self.h = np.interp(mpt,x_mean,h_mean)
    self.u = np.interp(mpt,x_mean,u_mean)

  def get_dispersion(self):
    ''' Evaluate the dispersion relation for our frequencies and water depth '''
    w = 2.*np.pi*self.f
    k_ = np.linspace(0,2,200)
    w_ = np.sqrt(self.g*k_*np.tanh(self.h*k_))
    k = np.interp(w,np.concatenate((-w_[::-1],w_)),np.concatenate((k_[::-1],k_)))
    k[k<=1e-12] = 1e-12
    self.c0 = w/k + self.u
    self.c0[self.c0<=1e-12] = 1e-12

  def pre_stretch(self):
    # Cut time
    self.tt = self.lags[self.lags>=0]
    self.nns = len(self.tt)
    self.f = np.fft.rfftfreq(self.nns,d=1./self.fs)
    self.nf = len(self.f)
    # Calculate dispersion
    self.get_dispersion()
    # Set stretching parameters
    self.nu = 100
    self.us = np.linspace(-1,1,self.nu)
    self.idx = np.logical_and(self.tt<=self.dr/4,self.tt>=self.dr/30)
    # Flip and calc reference
    self.cur_pos = self.pos[:,self.lags>=0]
    self.cur_pos = self.cur_pos[:,::-1]
    self.ref_pos = np.mean(self.cur_pos,axis=0)
    self.cur_neg = self.neg[:,self.lags<0]
    self.ref_neg = np.mean(self.cur_neg,axis=0)

  def stretch(self):
    ''' Create the dictionary of shifted reference traces '''
    ref_pos_f = np.fft.rfft(self.ref_pos)
    ref_neg_f = np.fft.rfft(self.ref_neg)
    sft_pos_f = np.zeros((self.nu,self.nf),dtype=np.complex_)
    sft_neg_f = np.zeros((self.nu,self.nf),dtype=np.complex_)
    for i,u in enumerate(self.us):
      ds = 1./(self.c0 + u) - 1./self.c0
      sft_pos_f[i,:] = ref_pos_f * np.exp(2j*np.pi*self.f*ds*self.dr)
      sft_neg_f[i,:] = ref_neg_f * np.exp(-2j*np.pi*self.f*ds*self.dr)
    self.sft_pos = np.fft.irfft(sft_pos_f,axis=1)
    self.sft_neg = np.fft.irfft(sft_neg_f,axis=1)
      
  def get_corrcoef(self):
    ''' Cross-correlate the trace at each hour with the reference traces '''
    self.cc_pos = np.zeros((self.nt,self.nu))
    self.cc_neg = np.zeros((self.nt,self.nu))
    for i in range(self.nt):
      for j in range(self.nu):
        self.cc_pos[i,j] = np.corrcoef(self.sft_pos[j,self.idx],self.cur_pos[i,self.idx])[1,0]
        self.cc_neg[i,j] = np.corrcoef(self.sft_neg[j,self.idx],self.cur_neg[i,self.idx])[1,0]
    self.cc = (self.cc_pos + self.cc_neg)/2 # Doppler-averaged correlation coefficient

  def get_current(self):
    self.current = np.zeros(self.nt)
    for it in range(self.nt):
      self.current[it] = self.us[self.cc[it,:].argmax()]

#######################################################################################



## Flag to save plots or not
save_plots = True

## Set geometry
off = 20 # 20-channel/200-m offset
srcs = np.arange(0,1000-200,10,dtype=int) # source channels for each subarray
U = []
M = []

## Load data, stretch, and find optimum CC
for m,src in enumerate(srcs):
  print('source %d/%d' % (m+1,len(srcs)))
  trxc = TL_gather(src,off)
  trxc.load_gather()
  trxc.load_mean()
  trxc.pre_stretch()
  trxc.stretch()
  trxc.get_corrcoef()
  trxc.get_current()

  U.append(trxc.current)
  M.append(np.mean((src,src+off))*trxc.dx) 
U = np.array(U)
M = np.array(M)

## Save
np.save('./disp/U_diff.npy',U)
np.save('./disp/midpts.npy',M)


## Plot 
fig1,ax = plt.subplots(2,1,figsize=(8,8),sharex=True,sharey=True)
ax[0].plot(trxc.tt,trxc.ref_pos,'k',label='reference')
ax[0].plot(trxc.tt,trxc.cur_pos[0,:],'r',label='hour 0')
ax[0].plot(trxc.tt,trxc.cur_pos[1,:],'b',label='hour 1')
ax[0].set_title('Anti-causal')
ax[0].legend()
ax[1].plot(trxc.tt,trxc.ref_neg,'k')
ax[1].plot(trxc.tt,trxc.cur_neg[0,:],'r')
ax[1].plot(trxc.tt,trxc.cur_neg[1,:],'b')
ax[1].set_xlim([trxc.dr/30,trxc.dr/4])
ax[1].set_xlabel('Time lag (s)')
ax[1].set_title('Causal')

fig2,ax = plt.subplots(1,3,figsize=(10,5))
ax[0].plot(trxc.us,trxc.cc_pos[0,:],'r',label='hour 0')
ax[0].plot(trxc.us,trxc.cc_pos[1,:],'b',label='hour 1')
ax[0].axvline(trxc.us[np.argmax(trxc.cc_pos[0,:])],c='r',linestyle='--')
ax[0].axvline(trxc.us[np.argmax(trxc.cc_pos[1,:])],c='b',linestyle='--')
ax[0].axvline(0,c='k')
ax[0].set_xlabel('Current velocity (m/s)')
ax[0].set_ylabel('Corr. coeff.')
ax[0].set_title('Anti-causal')
ax[0].legend()
ax[1].plot(trxc.us,trxc.cc_neg[0,:],'r')
ax[1].plot(trxc.us,trxc.cc_neg[1,:],'b')
ax[1].axvline(trxc.us[np.argmax(trxc.cc_neg[0,:])],c='r',linestyle='--')
ax[1].axvline(trxc.us[np.argmax(trxc.cc_neg[1,:])],c='b',linestyle='--')
ax[1].axvline(0,c='k')
ax[1].set_xlabel('Current velocity (m/s)')
ax[1].set_title('Causal')
ax[2].plot(trxc.us,trxc.cc[0,:],'r')
ax[2].plot(trxc.us,trxc.cc[1,:],'b')
ax[2].axvline(trxc.us[np.argmax(trxc.cc[0,:])],c='r',linestyle='--')
ax[2].axvline(trxc.us[np.argmax(trxc.cc[1,:])],c='b',linestyle='--')
ax[2].axvline(0,c='k')
ax[2].set_xlabel('Current velocity (m/s)')
ax[2].set_title('Doppler-averaged')


## Save plots
if save_plots:
  fig1.savefig('./figs/08a_traces.png')
  fig2.savefig('./figs/08b_corrcoeff.png')



plt.show()


