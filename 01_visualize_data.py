'''
01_visualize_data.py

Author: Ethan Williams (efwillia@caltech.edu)

- Load data file (.npy)
- Plot raw data
- Plot frequency-wavenumber spectrum
- Plot frequency spectrum along cable

'''

import numpy as np
import matplotlib.pyplot as plt


# Location of the data file
datafile = './data/Fs1Hz_dx10m_Belgium_data.npy'

# Flag to save plots or not
save_plots = True

# Parameters
nx = 1000                # number of channels
ns = 4200                # number of samples
fs = 1                   # sampling rate (Hz)
dx = 10                  # channel spacing (m)



## (A) Plot raw data ######################################################

# Load data
data = np.load(datafile)

# Create distance and time axis arrays
x = np.arange(nx) * dx
t = np.arange(ns) / fs

# Select a slice of the dataset
xmin = 0
xmax = 500
idx = np.logical_and(x>=xmin,x<=xmax)

tmin = 1000
tmax = 1200
idt = np.logical_and(t>=tmin,t<=tmax)

data_slice = data[np.ix_(idx,idt)]
xx = x[idx]
tt = t[idt]

# Plot the data
v = 1e2
fig1,ax = plt.subplots(1,2,figsize=(12,8))
ax[0].pcolormesh(xx,tt,data_slice.T,cmap='RdBu',vmin=-v,vmax=v)
ax[0].set_xlabel('Distance (m)')
ax[0].set_ylabel('Time (s)')
ax[0].set_title('Raw DAS data')



## (B) Plot frequency-wavenumber spectrum  ################################

# Take a larger slice of the dataset
xmin = 0
xmax = 10000
idx = np.logical_and(x>=xmin,x<=xmax)

tmin = 0
tmax = 600
idt = np.logical_and(t>=tmin,t<=tmax)

data_slice = data[np.ix_(idx,idt)]
xx = x[idx]
tt = t[idt]

# Apply two FFTs to get the FK spectrum
fk = np.fft.fft2(data_slice)
fk = np.fft.fftshift(fk)
fk = 20*np.log10(abs(fk))

# Get FK axes
f = np.fft.fftshift(np.fft.fftfreq(len(tt),d=1./fs))
k = np.fft.fftshift(np.fft.fftfreq(len(xx),d=dx))

# Get theoretical dispersion relation
h = 20 # water depth (m)
f_OSGW = np.sqrt(9.8*2*np.pi*k * np.tanh(2*np.pi*k*h)) * (0.5/np.pi)

# Plot FK spectrum 
ax[1].pcolormesh(f,k,fk,cmap='jet',vmin=60,vmax=120)
ax[1].plot(f_OSGW,k,'k--')
ax[1].set_xlabel('Frequency (1/s)')
ax[1].set_ylabel('Wavenumber (1/m)')
ax[1].set_title('FK spectrum')
ax[1].set_xlim([0,max(f)])



## (C) Plot spectrum along cable ##########################################

# For each channel, compute the spectrum
nns = ns//2 + 1
spec = np.zeros((nx,nns))
for ix in range(nx):
  spec[ix,:] = 20*np.log10(abs(np.fft.rfft(data[ix,:])))

# Get frequency axis
f = np.fft.rfftfreq(ns,d=1./fs)

# Plot spectrum along cable
fig2,ax = plt.subplots(1,1,figsize=(10,5))
ax.pcolormesh(x,f,spec.T,cmap='jet',vmin=10,vmax=90)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Frequency (Hz)')
ax.set_title('Spectrum vs. distance')



## Save
if save_plots:
  fig1.savefig('./figs/01a_raw_data.png')
  fig2.savefig('./figs/01b_spectrum.png')

## Show plots
plt.show()
