#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to make a plot from the 1D model using the 
hdf5 output expo.hdf5
The 1D model producing the expo files is in diffusive-advective-loss-v2.nb

Created on Mon Nov  4 17:49:57 2019

@author: bernard Legras
"""
import deepdish as dd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os

listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']            
mymap=colors.ListedColormap(listcolors)

figpdf=True
aEA = dd.io.load('expoEAcore.hdf5')['Dataset1'][...,0]
aEI = dd.io.load('expoEIcore.hdf5')['Dataset1'][...,0]
times = np.arange(0,70.1,0.25)
theta = np.arange(310,450.25,0.5)

degree = 6372 * 2 * np.pi / 360
ff_s = (1/24) * degree**2
fst = 16
fs  =20

#%%
# first figure showing the non normalized and the normalized impact on the same figure
fig = plt.figure(figsize=(9,9))
plt.subplot(2,2,1)
im1 = plt.imshow(aEA,extent=(0.,70,310,450),origin='lower',aspect='auto',cmap=mymap,
                   norm = LogNorm(vmin=10,vmax=5e7))
plt.ylim(320,420)
plt.xlim(0,62)
plt.ylabel(r'Potential temperature (K)',fontsize=fs)
plt.title(r'ERA5',fontsize=fs)
plt.tick_params(labelsize=fst)
plt.subplot(2,2,3)
im3 = plt.imshow(aEI,extent=(0.,70,310,450),origin='lower',aspect='auto',cmap=mymap,
                   norm = LogNorm(vmin=10,vmax=5e7))
plt.ylim(320,420)
plt.xlim(0,62)
plt.ylabel(r'Potential temperature (K)',fontsize=fs)
plt.xlabel(r'Age (day)',fontsize=fs)
plt.title('ERA-I',fontsize=fs)
plt.tick_params(labelsize=fst)

# Normalisation accounting for the 6h interval of the time axis
nEA = 4 * aEA / np.sum(aEA,axis=1)[:,np.newaxis]
nEI = 4 * aEI / np.sum(aEA,axis=1)[:,np.newaxis]

plt.subplot(2,2,2)
im2 = plt.imshow(nEA,extent=(0,70,310,450),origin='lower',aspect='auto',cmap=mymap,
                   norm = LogNorm(vmin=1.e-6,vmax=0.1))
plt.ylim(320,420)
plt.xlim(0,62)
#plt.ylabel(r'Potential temperature (K)',fontsize=fs)
plt.title(r'ERA5',fontsize=fs)
plt.tick_params(labelsize=fst)
plt.subplot(2,2,4)
im4 = plt.imshow(nEI,extent=(0,70,310,450),origin='lower',aspect='auto',cmap=mymap,
                   norm = LogNorm(vmin=1.e-6,vmax=0.1))
plt.ylim(320,420)
plt.xlim(0,62)
#plt.ylabel(r'Potential temperature (K)',fontsize=fs)
plt.xlabel(r'Age (day)',fontsize=fs)
plt.title('ERA-I',fontsize=fs)
plt.tick_params(labelsize=fst)

cax1 = fig.add_axes([0.15,-0.04,0.30,0.05])
cbar1 = fig.colorbar(im1,cax1,orientation='horizontal')
cbar1.set_label(r'Impact (day km$^2$ K$^{-1}$)',labelpad=-1,size=fs)
cbar1.ax.tick_params(labelsize=fs)
cax2 = fig.add_axes([0.57,-0.04,0.30,0.05])
cbar2 = fig.colorbar(im2,cax2,orientation='horizontal')
cbar2.set_label(r'Age spectrum  (day$^{-1}$)',labelpad=-1,size=fs)
cbar2.ax.tick_params(labelsize=fs)

if figpdf:
    plt.savefig(os.path.join('figs','1Dmodel_composit-ACP.png'),bbox_inches='tight',dpi=300)
    plt.savefig(os.path.join('figs','1Dmodel-composit-ACP.pdf'),bbox_inches='tight',dpi=300)
plt.show()

#%%
# second figure showing the non normalized impact only
fig = plt.figure(figsize=(9,4.5))
plt.subplot(1,2,1)
im1 = plt.imshow(aEA,extent=(0.,70,310,450),origin='lower',aspect='auto',cmap=mymap,
                   norm = LogNorm(vmin=10,vmax=5e7))
plt.ylim(320,420)
plt.xlim(0,62)
plt.ylabel(r'Potential temperature (K)',fontsize=fs)
plt.xlabel(r'Age (day)',fontsize=fs)
plt.title(r'ERA5',fontsize=fs)
plt.tick_params(labelsize=fst)
plt.subplot(1,2,2)
im3 = plt.imshow(aEI,extent=(0.,70,310,450),origin='lower',aspect='auto',cmap=mymap,
                   norm = LogNorm(vmin=10,vmax=5e7))
plt.ylim(320,420)
plt.xlim(0,62)
plt.xlabel(r'Age (day)',fontsize=fs)
plt.title('ERA-I',fontsize=fs)
plt.tick_params(labelsize=fst)

cax1 = fig.add_axes([0.15,-0.08,0.70,0.05])
cbar1 = fig.colorbar(im1,cax1,orientation='horizontal')
cbar1.set_label(r'Cumulated impact par age (day km$^2$ K$^{-1}$)',labelpad=-1,size=fs)
cbar1.ax.tick_params(labelsize=fs)

if figpdf:
    plt.savefig(os.path.join('figs','1Dmodel_composit-ACP-v2.png'),bbox_inches='tight',dpi=300)
    plt.savefig(os.path.join('figs','1Dmodel-composit-ACP-v2.pdf'),bbox_inches='tight',dpi=300)
plt.show()
