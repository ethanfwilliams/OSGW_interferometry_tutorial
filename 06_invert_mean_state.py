'''
06_mean_state.py

Author: Ethan Williams (efwillia@caltech.edu)

- Load dispersion picks for each subarray
- Invert for (U,h)

'''

import os
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


def dispersion(f,h,U):
    k_ = np.linspace(1e-3,1e0,200)
    w_ = U*k_ + np.sqrt(9.8*k_ * np.tanh(k_*h))
    w_2 = U*k_ - np.sqrt(9.8*k_ * np.tanh(k_*h))
    w_ = np.concatenate((w_2[::-1],w_),axis=0)
    k_ = np.concatenate((k_[::-1],k_),axis=0)
    k = np.interp(2*np.pi*f,w_,k_)
    return abs((2*np.pi*f)/k)

def residual(x,f,picks):
    return dispersion(f,x[0],x[1]) - picks


#######################################################################################


## Set paths
dispdir = './disp'
outdir = './disp'

## Load picks
freqs = np.load(os.path.join(dispdir,'dispersion_freqs.npy'))
picks = np.load(os.path.join(dispdir,'dispersion_picks.npy'))

## Set geometry
srcs = np.arange(0,1000-200,10,dtype=int)
mpts = (srcs + srcs + 200)/2
nx = picks.shape[0]
nf = picks.shape[1]

## Preallocate results
H = np.zeros(nx)
Hs = np.zeros(nx)
U = np.zeros(nx)
Us = np.zeros(nx)
model = np.zeros(picks.shape)

## Initial guess
x0 = (10,0.1) # 10 m depth, 0.1 m/s current

## For each subarray location
for ix in range(nx):
    # get rid of nans
    idf = ~np.isnan(picks[ix,:])
    f = freqs[idf]
    p = picks[ix,idf]
    # do Levenberg-Marquardt 
    result = least_squares(residual,x0,args=(f,p),method='lm')
    H[ix] = result.x[0]
    U[ix] = result.x[1]
    model[ix,idf] = dispersion(f,result.x[0],result.x[1])
    model[ix,~idf] = np.nan
    # calculate uncertainty
    J = result.jac
    cov = np.linalg.inv(J.T @ J)
    RSS = np.sum(result.fun**2)
    m = len(f); n = len(x0)
    std = np.sqrt(cov * RSS/(m-n))
    Hs[ix] = std[0,0]
    Us[ix] = std[1,1]

## Save
np.save(os.path.join(outdir,'H.npy'),np.column_stack((mpts,H,Hs)))
np.save(os.path.join(outdir,'U.npy'),np.column_stack((mpts,U,Us)))
np.save(os.path.join(outdir,'dispersion_model.npy'),model)

