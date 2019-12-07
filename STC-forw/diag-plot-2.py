#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the forward runs to determine the proportion of parcels going up 
with respect to that going down

This script exploits the results of diag-make-all

In this version the grouping is redone from histdisptheta and a 
latitude-altitude filter is applied to remove high latitude 
spurious sources.

The altitude where crossing from down to up displacement occurs is calculated
and printed for each region and group of regions.

Plots displacements and crossover properties for a set of regions
['Ocean','Land','Tibet','All','Asia','AfricaArabia']

1) display the cumul histogram of equal, above, below in the global domain for EID & EIZ

2) the same for FullAMA and EAD, EAZ

3) display the ratios to total in the global domain for EID & EIZ

4) the same for FullAMA and EAD, EAZ

5) cumul histograms and ratios for All region and ERA-I global (EIZ,EID), 
ERA-I FullAMA and ERA5 FullAMA (EAZ, EAD)

6) Total sources with respect to main Asia and AfricaArabia contributions

7) Proportion of total sources per region 

4Spe) Special version of figure 4 for Asia

Created on Tue Apr 30 22:52:48 2019

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
from os.path import join
import collections
from scipy.interpolate import Akima1DInterpolator
from scipy.optimize import brentq
import constants as cst
from group import group
#from numba import jit

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif') 

#supertype = 'EID-FULL'
#target = 'global'
#target = 'FullAMA'
vert = 'theta'
age_max = 62
age_max_inter = 30

# derived parameters
#run_type = supertype+'-Box-'+vert+'-'+target
hmax = 24*(age_max + 10)
hinter = 24*(age_max_inter + 10)
# read mask an initialize edges
print ('open mask')
with gzip.open(join('..','mkSTCmask','MaskCartopy2-STCfine.pkl'),'rb') as f: 
    mm =pickle.load(f) 
xedge = np.arange(-10,160+0.5*mm['icx'],mm['icx'])
yedge = np.arange(0,50+0.5*mm['icy'],mm['icx'])
minv = 322.5
maxv = 422.5
binv = 20
binvl = 5
vcent = np.arange(325,421,(maxv-minv)/binv)
savefig = False

# %%
diag = {}
diag['global'] = {}
diag['FullAMA'] = {}

#hdisptheta = {}
#hdisptheta['global'] = {}
#hdisptheta['FullAMA'] = 0

supertypes = {'global':['EID-FULL','EIZ-FULL'],'FullAMA':['EID-FULL','EIZ-FULL','EAD','EAZ','EAT']}
hightypes = ['mh','sh']

for target in ['global','FullAMA']:
    for supertype in supertypes[target]:
        # read the histogram for theta mean
        with gzip.open('diag-'+supertype+'-'+target+'-'+vert+'.pkl','rb') as f:
            [[_,hdisptheta],_] = pickle.load(f)
        # filtering of hdisptheta
        # set histog to zero above 360K in regions where spurious convection occurs
        for reg in ['Europe','NorthAsia','Atlantic','Mediterranea']:
            hdisptheta['sh'][: ,8:, mm['regcode'][reg]] = 0
            hdisptheta['mh'][: ,8:, mm['regcode'][reg]] = 0
        # set histog to zero above 370K for North China
        hdisptheta['mh'][:, 10:, mm['regcode']['NorthChina']] = 0 
        hdisptheta['sh'][:, 10:, mm['regcode']['NorthChina']] = 0
        # calculate the regional diagnostics
        diag[target][supertype] = {}    
        for hightype in hightypes:
            diag[target][supertype][hightype] = {}
            for var in ['equal','above','below','total']:
                diag[target][supertype][hightype][var] = np.zeros(shape=(binv,len(mm['regcode'])+1),dtype=np.float32)
            for ireg in range(1,len(mm['regcode'])+1):
                for lev in range(binv):
                    diag[target][supertype][hightype]['equal'][lev,ireg] = hdisptheta[hightype][lev,lev,ireg]
                    diag[target][supertype][hightype]['above'][lev,ireg] = np.sum(hdisptheta[hightype][lev+1:,lev,ireg])
                    diag[target][supertype][hightype]['below'][lev,ireg] = np.sum(hdisptheta[hightype][:lev,lev,ireg])
                    diag[target][supertype][hightype]['total'][lev,ireg] = np.sum(hdisptheta[hightype][:,lev,ireg])
            ddt = diag[target][supertype][hightype]['total'].copy()
            ddt[ddt==0] = 1
            diag[target][supertype][hightype]['equal_prop'] = diag[target][supertype][hightype]['equal']/ddt
            diag[target][supertype][hightype]['above_prop'] = diag[target][supertype][hightype]['above']/ddt
            diag[target][supertype][hightype]['below_prop'] = diag[target][supertype][hightype]['below']/ddt
    
            #%% Grouping of the regions into blocs
            for gr in group.keys():
                diag[target][supertype][hightype][gr] = {}
                for var in ['equal','above','below','total']:
                    diag[target][supertype][hightype][gr][var] = np.zeros(binv,dtype=np.float32)
                    for reg in group[gr]:
                        diag[target][supertype][hightype][gr][var] += diag[target][supertype][hightype][var][:,mm['regcode'][reg]]
                ddt = diag[target][supertype][hightype][gr]['total'].copy()
                ddt[ddt==0] = 1
                diag[target][supertype][hightype][gr]['equal_prop'] = diag[target][supertype][hightype][gr]['equal']/ddt
                diag[target][supertype][hightype][gr]['above_prop'] = diag[target][supertype][hightype][gr]['above']/ddt
                diag[target][supertype][hightype][gr]['below_prop'] = diag[target][supertype][hightype][gr]['below']/ddt
                del ddt

#%%
# find the crossings of the above and below curves for the different cases
# save the results in a pkl file
cross = {}
for target in ['global','FullAMA']:
    cross[target] = {}
    for supertype in supertypes[target]:
        cross[target][supertype] = {}
        for hightype in hightypes:
            cross[target][supertype][hightype] = {}
            for gr in group.keys():
                amb = diag[target][supertype][hightype][gr]['above_prop'] - diag[target][supertype][hightype][gr]['below_prop']
                fz = Akima1DInterpolator(vcent,amb)
                t1 = 355
                t2 = 370
                if gr=='Tibet':
                    if supertype == 'EAZ': t1 = 350
                    elif supertype == 'EIZ-FULL': t1 = 365
                    elif supertype == 'EAT': t1 = 350
                if gr=='AfricaArabia':
                    if (supertype == 'EIZ-FULL') & (target == 'FullAMA'):
                        t1 = 385
                        t2 = 400
                try:       
                    root = brentq(fz,t1,t2)
                    print(target,'\t',supertype,'\t','\t{:.2f} K\t'.format(root),hightype,'\t',gr)
                except:
                    print('ERROR',target,supertype,hightype,gr,fz(350),fz(375))
                    root = np.nan
                cross[target][supertype][hightype][gr] = root
with gzip.open('crossover.pkl','wb') as f:
    pickle.dump(cross,f)
#%% diag 
#im=plt.imshow(np.log10(hdisptheta_raw[...,6]),extent=(minv,maxv,minv,maxv),origin='lower',cmap='jet')
#plt.colorbar(im)
#plt.xlabel('source level')
#plt.ylabel('mean level')
#plt.show()
#    
                
#%% plot of the above, equal, below per levels
ff = 1/24
fs = 16
fig = plt.figure(figsize=(12,8))
fig.suptitle(r'1] source $\Delta \theta \pm$ global  EID and EIZ',fontsize=fs)
groups = ['Ocean','Land','Tibet','All','Asia','AfricaArabia']
splt = 0
target = 'global'
for gr in groups:
    splt += 1
    plt.subplot(2,3,splt)
    im = plt.semilogx(ff*diag[target]['EIZ-FULL']['sh'][gr]['equal'],vcent,'.c',
                           ff*diag[target]['EIZ-FULL']['sh'][gr]['below'],vcent,'--c',
                           ff*diag[target]['EIZ-FULL']['sh'][gr]['above'],vcent,'c',
                           ff*diag[target]['EID-FULL']['sh'][gr]['equal'],vcent,'.b',
                           ff*diag[target]['EID-FULL']['sh'][gr]['below'],vcent,'--b',
                           ff*diag[target]['EID-FULL']['sh'][gr]['above'],vcent,'b',
                           linewidth=4,alpha=0.9 )
    plt.ylim(330,400)
    if splt ==3: plt.xlim(0.01,300)
    else: plt.xlim(0.1,3000)
    plt.tick_params(labelsize=fs)
    plt.title(gr,fontsize=fs)
plt.show()
    
fig = plt.figure(figsize=(12,8))
fig.suptitle(r'2] source $\Delta \theta \pm$ FullAMA  EID and EAD',fontsize=fs)
groups = ['Ocean','Land','Tibet','All','Asia','AfricaArabia']
splt = 0    
target = 'FullAMA'
for gr in groups:
    splt += 1
    plt.subplot(2,3,splt)
    im = plt.semilogx(
                           #ff*diag[target]['EIZ-FULL']['sh'][gr]['equal'],vcent,'.c',
                           #ff*diag[target]['EIZ-FULL']['sh'][gr]['below'],vcent,'--c',
                           #ff*diag[target]['EIZ-FULL']['sh'][gr]['above'],vcent,'c',
                           ff*diag[target]['EID-FULL']['sh'][gr]['equal'],vcent,'.b',
                           ff*diag[target]['EID-FULL']['sh'][gr]['below'],vcent,'--b',
                           ff*diag[target]['EID-FULL']['sh'][gr]['above'],vcent,'b',
                           #ff*diag[target]['EAZ']['sh'][gr]['equal'],vcent,'.m',
                           #ff*diag[target]['EAZ']['sh'][gr]['below'],vcent,'--m',
                           #ff*diag[target]['EAZ']['sh'][gr]['above'],vcent,'m',
                           ff*diag[target]['EAD']['sh'][gr]['equal'],vcent,'.r',
                           ff*diag[target]['EAD']['sh'][gr]['below'],vcent,'--r',
                           ff*diag[target]['EAD']['sh'][gr]['above'],vcent,'r',
                           linewidth=4,alpha=0.9 )
    plt.ylim(330,400)
    if splt ==3: plt.xlim(0.01,300)
    else: plt.xlim(0.1,3000)
    plt.tick_params(labelsize=fs)
    plt.title(gr,fontsize=fs)   
plt.show()


#%% plot of the above, equal, below ratios per levels
ff = 1/24
fs = 16
fig = plt.figure(figsize=(12,8))
fig.suptitle(r'3] source $\Delta \theta \pm$ ratios global  EID and EIZ',fontsize=fs)
groups = ['Ocean','Land','Tibet','All','Asia','AfricaArabia']
splt = 0
target = 'global'
for gr in groups:
    splt += 1
    plt.subplot(2,3,splt)
    im = plt.plot(ff*diag[target]['EIZ-FULL']['sh'][gr]['equal_prop'],vcent,'.c',
                           diag[target]['EIZ-FULL']['sh'][gr]['below_prop'],vcent,'--c',
                           diag[target]['EIZ-FULL']['sh'][gr]['above_prop'],vcent,'c',
                           diag[target]['EID-FULL']['sh'][gr]['equal_prop'],vcent,'.b',
                           diag[target]['EID-FULL']['sh'][gr]['below_prop'],vcent,'--b',
                           diag[target]['EID-FULL']['sh'][gr]['above_prop'],vcent,'b',
                           linewidth=4,alpha=0.9 )
    plt.ylim(330,400)
    plt.tick_params(labelsize=fs)
    plt.title(gr,fontsize=fs)
plt.show()

#%%    
fig = plt.figure(figsize=(12,8))
fig.suptitle(r'4] source $\Delta \theta \pm$ ratios FullAMA  EAZ and EAD',fontsize=fs)
groups = ['Ocean','Land','Tibet','All','Asia','AfricaArabia']
splt = 0    
target = 'FullAMA'
for gr in groups:
    splt += 1
    plt.subplot(2,3,splt)
    im = plt.plot(
                           #diag[target]['EIZ-FULL']['sh'][gr]['equal'],vcent,'.c',
                           #diag[target]['EIZ-FULL']['sh'][gr]['below'],vcent,'--c',
                           #diag[target]['EIZ-FULL']['sh'][gr]['above'],vcent,'c',
                           diag[target]['EID-FULL']['sh'][gr]['equal_prop'],vcent,'.b',
                           diag[target]['EID-FULL']['sh'][gr]['below_prop'],vcent,'--b',
                           diag[target]['EID-FULL']['sh'][gr]['above_prop'],vcent,'b',
                           diag[target]['EAZ']['sh'][gr]['equal_prop'],vcent,'.m',
                           diag[target]['EAZ']['sh'][gr]['below_prop'],vcent,'--m',
                           diag[target]['EAZ']['sh'][gr]['above_prop'],vcent,'m',
                           diag[target]['EAD']['sh'][gr]['equal_prop'],vcent,'.r',
                           diag[target]['EAD']['sh'][gr]['below_prop'],vcent,'--r',
                           diag[target]['EAD']['sh'][gr]['above_prop'],vcent,'r',
                           linewidth=4,alpha=0.9 )
    plt.ylim(330,400)
    plt.tick_params(labelsize=fs)
    plt.title(gr,fontsize=fs)   
plt.show()

#%% plots for All   
# plot of the cumul and of the ratios
ff = 1/24
fs = 16
fig = plt.figure(figsize=(12,8))
fig.suptitle(r'5] Cumul and ratios for All',fontsize=fs)
splt = 0
for target in ['global','FullAMA']:
    splt += 1
    plt.subplot(2,3,splt)            
    im = plt.semilogx(ff*diag[target]['EIZ-FULL']['sh']['All']['equal'],vcent,'.b',
                           ff*diag[target]['EIZ-FULL']['sh']['All']['below'],vcent,'--b',
                           ff*diag[target]['EIZ-FULL']['sh']['All']['above'],vcent,'b',
                           ff*diag[target]['EID-FULL']['sh']['All']['equal'],vcent,'.r',
                           ff*diag[target]['EID-FULL']['sh']['All']['below'],vcent,'--r',
                           ff*diag[target]['EID-FULL']['sh']['All']['above'],vcent,'r',
                           linewidth=4,alpha=0.9)
    plt.ylim(330,400)
    plt.xlim(1,10000)
    plt.tick_params(labelsize=fs)
    plt.title('ERA-I '+target,fontsize=fs)
    plt.subplot(2,3,splt+3)            
    im = plt.plot(diag[target]['EIZ-FULL']['sh']['All']['equal_prop'],vcent,'.b',
                           diag[target]['EIZ-FULL']['sh']['All']['below_prop'],vcent,'--b',
                           diag[target]['EIZ-FULL']['sh']['All']['above_prop'],vcent,'b',
                           diag[target]['EID-FULL']['sh']['All']['equal_prop'],vcent,'.r',
                           diag[target]['EID-FULL']['sh']['All']['below_prop'],vcent,'--r',
                           diag[target]['EID-FULL']['sh']['All']['above_prop'],vcent,'r',
                           linewidth=4,alpha=0.9)
    plt.ylim(330,400)
    plt.tick_params(labelsize=fs)
    #plt.title(gr,fontsize=fs)
target = 'FullAMA'
plt.subplot(2,3,3)            
im = plt.semilogx(ff*diag[target]['EAZ']['sh']['All']['equal'],vcent,'.b',
                       ff*diag[target]['EAZ']['sh']['All']['below'],vcent,'--b',
                       ff*diag[target]['EAZ']['sh']['All']['above'],vcent,'b',
                       ff*diag[target]['EAD']['sh']['All']['equal'],vcent,'.r',
                       ff*diag[target]['EAD']['sh']['All']['below'],vcent,'--r',
                       ff*diag[target]['EAD']['sh']['All']['above'],vcent,'r',
                       linewidth=4,alpha=0.9)
plt.ylim(330,400)
plt.xlim(1,10000)
plt.tick_params(labelsize=fs)
plt.title('ERA5 FullAMA',fontsize=fs)
plt.subplot(2,3,6)            
im = plt.plot(diag[target]['EAZ']['sh']['All']['equal_prop'],vcent,'.b',
                       diag[target]['EAZ']['sh']['All']['below_prop'],vcent,'--b',
                       diag[target]['EAZ']['sh']['All']['above_prop'],vcent,'b',
                       diag[target]['EAD']['sh']['All']['equal_prop'],vcent,'.r',
                       diag[target]['EAD']['sh']['All']['below_prop'],vcent,'--r',
                       diag[target]['EAD']['sh']['All']['above_prop'],vcent,'r',
                       linewidth=4,alpha=0.9)

plt.tick_params(labelsize=fs)
plt.show()
#plt.title(gr,fontsize=fs)
#%% Test
# Total sources with respect to the contribution of Asia an Africa
fig.suptitle(r'6] total sources All, Asia, AficaArabia',fontsize=fs)
plt.semilogx(ff*diag[target]['EAZ']['sh']['All']['total'],vcent,'r',
             ff*diag[target]['EAZ']['sh']['Asia']['total'],vcent,'b',
             ff*diag[target]['EAZ']['sh']['AfricaArabia']['total'],vcent,'k',)
plt.ylim(330,400)
#%%
# Total source per region
fig = plt.figure(figsize=(7,21))
fig.suptitle(r'7] total sources proportion per region',fontsize=fs)
splt = 0
for reg in mm['regcode'].keys():
    splt += 1
    plt.subplot(9,3,splt)
    plt.plot(diag[target]['EAZ']['sh']['total'][:,mm['regcode'][reg]]/diag[target]['EAZ']['sh']['All']['total'],vcent,'r')
    plt.ylim(330,410)
    plt.title('total '+reg)
    plt.xlim(0,0.3)
plt.show()
#%% special view Asia
#%%    
fig = plt.figure(figsize=(6,6))
#fig.suptitle(r'4spe] source $\Delta \theta \pm$ ratios FullAMA EAD and global EID',fontsize=fs) 
fs=22
gr = 'Asia'
im = plt.plot(
       #diag[target]['EIZ-FULL']['sh'][gr]['equal'],vcent,'.c',
       #diag[target]['EIZ-FULL']['sh'][gr]['below'],vcent,'--c',
       #diag[target]['EIZ-FULL']['sh'][gr]['above'],vcent,'c',
       diag['global']['EID-FULL']['sh'][gr]['equal_prop'],vcent,'Db',
       diag['global']['EID-FULL']['sh'][gr]['below_prop'],vcent,'--b',
       diag['global']['EID-FULL']['sh'][gr]['above_prop'],vcent,'b',
       diag[target]['EAD']['sh'][gr]['equal_prop'],vcent,'Dk',
       diag[target]['EAD']['sh'][gr]['below_prop'],vcent,'--k',
       diag[target]['EAD']['sh'][gr]['above_prop'],vcent,'k',
       diag['FullAMA']['EID-FULL']['sh'][gr]['equal_prop'],vcent,'Dr',
       diag['FullAMA']['EID-FULL']['sh'][gr]['below_prop'],vcent,'--r',
       diag['FullAMA']['EID-FULL']['sh'][gr]['above_prop'],vcent,'r',
       linewidth=5,markersize=7,alpha=0.8)
plt.ylim(330,400)
plt.ylabel('Potential temperature (K)',fontsize=fs)
plt.xlabel('Proportion',fontsize=fs)
plt.tick_params(labelsize=fs)
plt.title('Asia crossover',fontsize=fs)
if savefig:
    plt.savefig(join('figs-Box','Asia-crossover.pdf'),dpi=300,bbox_inches='tight')
    plt.savefig(join('figs-Box','Asia-crossover.png'),dpi=300,bbox_inches='tight')
plt.show()
