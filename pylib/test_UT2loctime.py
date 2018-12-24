#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script tests the transformation of UTC to localtime in UT2locTime.
The test is based on a series of functions that are defined as functions of local time.
Cases 2 and 4 show the effect of interpolation on discontinuous functions

Created on Sun Jun 17 11:08:03 2018

@author: Bernard Legras
"""

from UT2locTime import UT2locTime
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import sawtooth

# Synthetic signal
N = 360
M = 8
MU = np.zeros(shape=(5,N,M))
lambda0 = 0
for m in range(M):
    UT = 24/M * m
    LT = UT + 24/360 * lambda0 + 24*np.arange(N)/N  
    MU[0,:,m] = np.cos(2*np.pi*LT/24)
    MU[1,:,m] = np.cos(4*np.pi*LT/24)+np.sin(2*np.pi*(LT/24))
    MU[2,:,m] = sawtooth(2*np.pi*LT/24,1)
    MU[3,:,m] = sawtooth(2*np.pi*LT/24,0.5)
    MU[4,:,m] = 2*np.heaviside(np.cos(2*np.pi*LT/24),0.5)-1
ML = UT2locTime(MU,0)

# case 0
plt.subplot(1,2,1)
for m in range(M):
    plt.plot(ML[0,:,m])
plt.plot(MU[0,:,0])
plt.subplot(1,2,2)
plt.plot(np.mean(ML[0,:,:],axis=1))
plt.plot(np.mean(MU[0,:,:],axis=1))
plt.show()

# case 1
plt.subplot(1,2,1)
for m in range(M):
    plt.plot(ML[1,:,m])
plt.plot(MU[1,:,0])
plt.subplot(1,2,2)
plt.plot(np.mean(ML[1,:,:],axis=1))
plt.plot(np.mean(MU[1,:,:],axis=1))
plt.show()

# case 2
plt.subplot(1,2,1)
for m in range(M):
    plt.plot(ML[2,:,m])
plt.plot(MU[2,:,0])
plt.subplot(1,2,2)
plt.plot(np.mean(ML[2,:,:],axis=1))
plt.plot(np.mean(MU[2,:,:],axis=1))
plt.show()

# case 3
plt.subplot(1,2,1)
for m in range(M):
    plt.plot(ML[3,:,m])
plt.plot(MU[3,:,0])
plt.subplot(1,2,2)
plt.plot(np.mean(ML[3,:,:],axis=1))
plt.plot(np.mean(MU[3,:,:],axis=1))
plt.show()

# case 4
plt.subplot(1,2,1)
for m in range(M):
    plt.plot(ML[4,:,m])
plt.plot(MU[4,:,0])
plt.subplot(1,2,2)
plt.plot(np.mean(ML[4,:,:],axis=1))
plt.plot(np.mean(MU[4,:,:],axis=1))
plt.show()