#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function that convert an array calculated with an index in UTC to an array in
local time
input: 
    MU: array in UTC, the last two indexes must be in order longitude and the
    UTC time
    lambda0: origin of the longitude grid (degree)
    N: optional parameter giving the total number of longitude points around a 
    latitude circle if the array does not cover the full circle
    It is assumed that the longitudes of the domain run from lambda0 with a step 360/N
    
output:
    ML: array with same dimension as MU but with last index in local time
    The local time is discretized in the same way as the UTC

    A necessary condition is that N is an integer multiple of the number of 
    time slots that discretize both UTC and local time. The routine fails and
    returns error code otherwise.

    The linear interpolatio is such that the sum over all times and the varaiance 
    are preserved at any longitude

    MU can have any number of dimensions before longitude and time      
    
    see notes in h2/manus/Methods/LocalTimeInterp-v3.pdf
    
Created on Sun Jun 17 02:30:48 2018

@author: Bernard Legras
"""
import numpy as np

def UT2locTime(MU,lambda0,N=None):
    # Find M and N
    dims = MU.shape
    M = dims[-1]
    Nl = dims[-2]
    if N==None:
        N = Nl
    if Nl>N:
        print ('Error on the longitude index')
        return -1
    # Test whether N/M is an integer
    if int(N/M) != N/M:
        print('N not a multiple of M')
        # should better raise an exception
        return -1
    J = int(N/M)
    # Set m0 (with offset to avoid truncation errors)
    m0 = int(lambda0*N/360+0.0001)
    # Define nstar
    nstar = - m0 - J*np.arange(M) %N
    # Initialize ML
    ML = np.empty(shape=dims)
    # Loops on p and s=q-p
    for p in range(M):
        for s in range(M):
            # Define range of longitude indexes for this segment
            nn = nstar[s] + range(J) %N
            # Select part of the segment within the domain
            sel = nn<Nl
            if len(sel)>0:
                nns = nn[sel]
                js = np.arange(J)[sel]
                ML[...,nns,p] = (1-js/J)*MU[...,nns,p+s %M] + js/J*MU[...,nns,p+s+1 %M]
    return ML