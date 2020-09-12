#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script provides statistics of the forward runs of the FORWBox-meanhigh family

It is based on the output of ageStat.py and is a variant of plotAgeStat
In plotAgeStat, the diagnostics are shown for all streams for a given supertype

This version is a special version of plotAgeStatComposit for the 2019 STC final meeting
It differs by removing the plots for the Return runs, rescaling due to change of 
time unit from hour to day and accounting for the output step and beautifying 
the figures, and improving the printed version.
It differs from plotAgeStatComposit4ACP by changing spatial units from degree to km

Here we show the results for the whole summer 2017 and as a comparison 
between the supertypes in the FullAMA domain (EAZ, EAD, EIZ, EID, EIZ-Return and EID-Return)

It generates for both sh and mh hightypes 4  figures which are
(actually mh not plotted)

1) The 2d histogram of parcels as a function of age and potential temperature

2) The same histogram normalized per level to better see the vertical propagation

3) The mean age and the modal age for each level

4) The number of active parcels as a function of age

5) Plot of the vertical distribution of sources at (age 0 parcels)

6) Plot for each of the supertype of the normalized age histogram at 370 K, 380 K and 400 K ,
that is positions 94, 104, 124
Show mean and modal peak

7) Comparison of the histograms for EAZ, EAD, EID, EIZ, EID-Return, EIZ-Return
superimposed to the source curve

7bis) Plot of the proportion of parcels above a given level relative to the same 
proportion in the sources

This is made for each of the 9 decades and can be called for all the runs

EAZ, EAD, EIZ, EID are showing statistics of parcels in the FullAMA domain with the rule that a particle that exits
once is discarded.
EIZ-FULL, EID-FULL show the statistics of parcels in the global domain
EIZ-Return, EID-Return show the statistics of parcels in the FullAMA domain by applying a simple mask to the FULL runs, 
that is parcels that leave the domain are counted if they return.  

Modifications to show EAT diagnostics: 27 August 2019

Created on Sat 23 Feb 2019

@author: Bernard Legras
"""
import numpy as np
#from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
import argparse
import socket
import os

# Color list with 20 colors      
listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']            
mymap=colors.ListedColormap(listcolors)
dpi = 300

parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",type=str,choices=["EAD","EAZ","EIZ","EID","EIZ-Return","EID-Return","EIZ-FULL","EID-FULL"],help="type")
#parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')

#supertypes = ["EAZ","EIZ","EIZ-Return","EIZ-FULL","EAD","EID","EID-Return","EID-FULL","EAT"]
#hightypes = ['sh','mh']
hightypes = ['sh',]

step = 6
hmax = 1728
# 62 days
age_max = 1488
nstep  =  int(hmax/step)
figsave = False
figpdf = False
nbins = 425
theta = np.arange(275.5,700,1)
ages = np.arange(0,age_max+1,6)/24
ageaxis = np.arange(0.,62.25,0.25) 
thetaxis = np.arange(275.5,700)

fs = 16
trea = {'EAZ':'ERA5 kinematic','EAD':'ERA5 diabatic','EAT':'ERA5 total diabatic',
        'EIZ':'ERA-I kinematic','EID':'ERA-I diabatic','EID-Return':'ERA-I diabatic return'}
trea_short = {'EAZ':'ERA5 kin','EAD':'ERA5 dia','EAT':'ERA5 tot',
        'EIZ':'ERA-I kin','EID':'ERA-I dia','EID-Return':'ERA-I dia ret'}
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 300

supertypes = ['EAT','EAD','EID','EID-Return']

args = parser.parse_args()
if args.type is not None: supertype = args.type

dates = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]
#dates = ["Jun-01","Jun-11","Jun-21"]
#dates = ['Aug-21']

# Find the main work dir according to the computer
if socket.gethostname() == 'gort':
    forw_dir =  '/home/legras/data/STC/STC-forw'
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif 'ciclad' in socket.gethostname():
    forw_dir =  '/home/legras/STC/STC-forw'
    out_dir = '/data/legras/STC/STC-FORWBox-meanhigh-OUT'
elif 'satie' in socket.gethostname():
    forw_dir =  '/home/legras/data/STC/STC-forw'
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif 'Graphium' in socket.gethostname():
    forw_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-forw"
    out_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-FORWBox-meanhigh-OUT" 

result = {}
result['sh'] = {}
result['mh'] = {}
for supertype in supertypes:
    for hightype in hightypes:       
        result[hightype][supertype] = {}
        result[hightype][supertype]['histog'] = np.zeros(shape=(nstep,nbins))
    for date in dates:
        print('')
        print('Processing '+date+' for '+supertype)         
        # get the histograms and mean data calculated by checkmean.py from the part files
        file_out = os.path.join(out_dir,'ageStat-'+supertype+'-2017-'+date)
        #if supertype == 'EAD':
        #    file_out = os.path.join(out_dir,'ageStat-'+supertype+'-Core-2017-'+date)
        print('reading',file_out )
        with gzip.open(file_out,'rb') as f:
           [_,histog,_,_] = pickle.load(f)
        for hightype in hightypes:
            result[hightype][supertype]['histog'] += histog[hightype]
    for hightype in hightypes:
        result[hightype][supertype]['nactiv'] = np.sum(result[hightype][supertype]['histog'],axis=1)

#%%        
# need to load the true source distribution
with gzip.open(os.path.join(forw_dir,'source_dist_updated.pkl'),'rb') as f:
        source_dist = pickle.load(f)

# insert the source distribution       
for hightype in hightypes:
    for supertype in supertypes:
        result[hightype][supertype]['histog'] = np.insert(result[hightype][supertype]['histog'],0,source_dist['sh'],axis=0)

#%% 1) Plot of the age/thet histogram
#for hightype in hightypes:
#    fig = plt.figure(figsize=(14,7))
#    fig.suptitle(hightype+' age theta histogram' )
#    n = 1
#    for supertype in supertypes:
#        plt.subplot(2,4,n)
#        plt.imshow(np.log10(result[hightype][supertype]['histog'][0:248,:]).T,
#                   extent=(0.25,62,275,700),origin='lower',aspect='auto',cmap='jet',clim=(-1,5))
#        plt.ylim(310,450)
#        plt.colorbar( in [1,5]: plt.ylabel('potential temperature (K)')
#        if n in [5,6,7,8]: plt.xlabel('age (day)')
#        n += 1
#    if figsave:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethethist-composit'))
#    plt.show()

#%% 1ACP) Plot of the age/thet histogram
# normalization : conversion of time into days and degree in km
# The factor Delta theta is 1 here (1K layer)
degree = 6372 * 2 * np.pi / 360
ff_s = (1/24) * degree**2
fst = 16
fs  =20

for hightype in hightypes:
    fig = plt.figure(figsize=(9,9))
    #fig.suptitle(hightype+' age theta histogram',fontsize=fs )
    n = 1
    for supertype in supertypes:
        plt.subplot(2,2,n)
        #im = plt.imshow(np.log10(ff_s*result[hightype][supertype]['histog'][0:248,:]).T,
        #           extent=(0.25,62,275,700),origin='lower',aspect='auto',cmap=mymap,clim=(-3,3))
        im = plt.imshow(ff_s*result[hightype][supertype]['histog'][0:249,:].T,
                   extent=(0.,62,275,700),origin='lower',aspect='auto',cmap=mymap,
                   norm = LogNorm(vmin=10,vmax=5e7))
        plt.ylim(320,440)
        plt.xlim(0,62)
        plt.tick_params(labelsize=fst)
        #plt.colorbar()
        plt.title(trea[supertype],fontsize=fs)      
        if n in [1,3]: plt.ylabel('potential temperature (K)',fontsize=fs)
        if n in [3,4]: plt.xlabel('age (day)',fontsize=fs)
        n += 1
    cax = fig.add_axes([0.17,-0.04,0.67,0.05])
    cbar = fig.colorbar(im,cax,orientation='horizontal')
    cbar.set_label(r'cumulative impact per age (day km$^2$ K$^{-1}$)',labelpad=-1,size=fs)
    cbar.ax.tick_params(labelsize=fs)
    if figpdf:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethethist-composit-withEAT-ACP.png'),bbox_inches='tight',dpi=dpi)
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethethist-composit-withEAT-ACP.pdf'),bbox_inches='tight',dpi=dpi)
    plt.show()
    
#%% 1 bis ACP) diagnostics on this histogram
    # integral of source dist over the potential temperature
    d_theta = 5
    SS = ff_s*np.sum(source_dist['sh']) * d_theta
    tc = 140
    for hightype in hightypes:
        impact = {}
        impact_age = {}
        impact_theta = {}
        for supertype in supertypes:            
            impact[supertype] = ff_s*result[hightype][supertype]['histog'][0:249,:].T
            # integral of impact over theta as a function of age
            impact_age[supertype] = np.sum(impact[supertype],axis=0) * d_theta
            # proba to each the level theta
            impact_theta[supertype] = np.sum(impact[supertype],axis=1)*d_theta/SS
        fig = plt.figure(figsize=(12,6))    
        plt.subplot(1,2,1)
        ages = np.arange(0.,62.1,0.25)
        plt.semilogy(ages,impact_age['EID-Return'],'--r',
                            ages,impact_age['EID'],'r',
                            ages,impact_age['EAT'],'--b',
                            ages,impact_age['EAD'],'b',linewidth=6)
        #plt.legend(('EIZ 17 days','EID 15 days','EAZ 13.3 days','EAD 13.3 days'),fontsize=fs,loc='upper right')
        # fit slope for ages > 35 days
        for supertype in supertypes:
            [slope,i0] = np.polyfit(ages[tc:],np.log(impact_age[supertype][tc:]),1)
            print ('decay time and anchor',supertype,1/slope,np.exp(i0))
            plt.semilogy(ages[tc:],np.exp(i0 + slope*ages[tc:]),'k',linewidth=12,alpha=0.3)
        plt.xlabel('ages (day)',fontsize=fs)
        plt.ylabel('integrated impact',fontsize=fs)
        plt.tick_params(labelsize=fst)
        plt.subplot(1,2,2)
        thetas = np.arange(275.5,700,1)
        plt.semilogx(impact_theta['EID-Return'],thetas,'--r',
                            impact_theta['EID'],thetas,'r',
                            impact_theta['EAT'],thetas,'--b',
                            impact_theta['EAD'],thetas,'b',linewidth=6)
        #plt.legend(('EIZ 10.6 K','EID 20.1 K','EAZ 14.7 K','EAD 15.4 K'),fontsize=fs,loc='lower left')
        # fit slope for lveles between 370K and 410K
        for supertype in supertypes:
            [slope,pp] = np.polyfit(thetas[95:136],np.log(impact_theta[supertype][95:136]),1)
            print('proba theta time and anchor',supertype,1/slope,np.exp(pp))
            plt.semilogx(np.exp(pp + slope*thetas[95:136]),thetas[95:136],'k',linewidth=12,alpha=0.3)
        plt.xlabel('probability of hiting the level',fontsize=fs)
        plt.ylabel('potential temperature (K)',fontsize=fs)
        plt.tick_params(labelsize=fs)
        plt.ylim(320,440)
        plt.xlim(1.e-3,2)
        plt.show()
    
#%% 2) Plot of the age/thet histogram, normalized for each level
#for hightype in hightypes
#    fig = plt.figure(figsize=(14,7))
#    fig.suptitle(hightype+' age theta normalized histogram' )
#    n = 1
#    for supertype in supertypes:
#        plt.subplot(2,4,n)
#        hh = result[hightype][supertype]['histog'][0:248,:]
#        ss = np.sum(hh,axis=0)
#        hh = hh / ss[np.newaxis,:]
#        plt.imshow(np.log10(hh.T),
#                   extent=(0.25,62,275,700),origin='lower',aspect='auto',cmap='jet',clim=(-6,0))
#        plt.ylim(310,450)
#        plt.colorbar()
#        plt.title(supertype)      
#        if n in [1,5]: plt.ylabel('potential temperature (K)')
#        if n in [5,6,7,8]: plt.xlabel('age (day)')
#        n += 1
#    if figsave:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethetnormhist-composit'))
#    plt.show()
    
#%% 2ACP) Plot of the age/thet histogram, normalized for each level by integrating in time
lt1 = 94   # index for theta = 369.5 K
lt2 = 144  # index for theta = 419.5 K
for hightype in hightypes:
    fig = plt.figure(figsize=(9,9))
    #fig.suptitle(hightype+' age theta normalized histogram' )
    n = 1
    for supertype in supertypes:
        plt.subplot(2,2,n)
        hh = result[hightype][supertype]['histog'][0:249,:].copy()
        # factor 0.25 because the sampling step is 6h
        ss = 0.25*np.sum(hh,axis=0)
        hh = hh / ss[np.newaxis,:]
        im = plt.imshow(hh.T,norm = LogNorm(vmin=1e-6,vmax=0.1),
                   extent=(0.,62,275,700),origin='lower',aspect='auto',cmap=mymap)
        # 0.25 is 6h
        # The offset 30 is here to avoid capturing a wrong max at high altitude in EAD and EIZ
        [slope,org] = np.polyfit(ages[30+np.argmax(hh[30:,lt1:lt2+1],axis=0)],theta[lt1:lt2+1],1)
        age1 = (theta[lt1] - org)/slope
        age2 = (theta[lt2] - org)/slope
        plt.plot([age1,age2],[theta[lt1],theta[lt2]])
        plt.ylim(320,420)
        plt.xlim(0.,62)
        plt.title(trea[supertype],fontsize=fs)
        print(supertype,slope)
        
        if n in [1,3]: plt.ylabel('potential temperature (K)',fontsize=fs)
        if n in [3,4]: plt.xlabel('age (day)',fontsize=fs)
        n += 1
    cax = fig.add_axes([0.17,-0.04,0.67,0.05])
    cbar = fig.colorbar(im,cax,orientation='horizontal')
    cbar.set_label(r'normalized age spectrum per level (day$^{-1}$)',labelpad=-1,size=fs)
    cbar.ax.tick_params(labelsize=fst)
    if figpdf:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethetnormhist-composit-withEAT-ACP.png'),bbox_inches='tight',dpi=dpi)
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethetnormhist-composit-withEAT-ACP.pdf'),bbox_inches='tight',dpi=dpi)
    plt.show()
    
#%% 3) Mean age
#ageaxis = np.arange(0.25,62.25,0.25) 
#thetaxis = np.arange(275.5,700)
#for hightype in hightypes:
#    fig = plt.figure(figsize=(14,7))
#    fig.suptitle(hightype+'  mean and modal age ')
#    n = 1
#    for supertype in supertypes:
#        plt.subplot(2,4,n)
#        hh = result[hightype][supertype]['histog'][0:248,:]
#        ss = np.sum(hh,axis=0)
#        ss[ss==0]=1
#        hh = hh / ss[np.newaxis,:]
#        agemean = np.sum(hh*ageaxis[:,np.newaxis],axis=0)
#        agemode = ageaxis[np.argmax(hh,axis=0)]
#        plt.plot(agemean,thetaxis,'k',agemode,thetaxis,'r')
#        plt.ylim(325,425)
#        plt.title(supertype)
#        if n in [1,5]: plt.ylabel('potential temperature (K)')
#        if n in [5,6,7,8]: plt.xlabel('age (day)')
#        n +=1
#    if figsave:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agemean-composit'))
#    plt.show()

#%% Plot of the number of active parcels as a function of age
#for hightype in hightypes:
#    fig = plt.figure(figsize=(14,7))
#    fig.suptitle(hightype+' nactiv as a function of age')
#    n = 1
#    for supertype in supertypes:
#        plt.subplot(2,4,n)
#        plt.semilogy(ageaxis,result[hightype][supertype]['nactiv'][0:248])
#        plt.title(supertype)
#        if n in [1,5]: plt.ylabel('nactiv')
#        if n in [5,6,7,8]: plt.xlabel('age (day)')
#        n += 1
#    if figsave:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agenactiv-composit'))
#    plt.show()
#    
#%% 4) Plot of the number of active parcels as a function of age in three separate 
# vertical bands of theta (and total)
#for hightype in hightypes:
#    fig = plt.figure(figsize=(14,7))
#    fig.suptitle(hightype+' nactiv as a function of age  per range tot(k) top(r) bot(b) mid(g)')
#    n = 1
#    for supertype in supertypes:
#        plt.subplot(2,4,n)
#        # 65 is first thetaxis above 340 K and 95 is first thetaxis above 370 K
#        total = np.sum(result[hightype][supertype]['histog'][0:248,:],axis=1)
#        toprange = np.sum(result[hightype][supertype]['histog'][0:248,95:],axis=1)
#        botrange = np.sum(result[hightype][supertype]['histog'][0:248,0:65],axis=1)
#        midrange = np.sum(result[hightype][supertype]['histog'][0:248,65:95],axis=1)
#        plt.semilogy(ageaxis,total,'k',ageaxis,toprange,'r',
#                     ageaxis,botrange,'b',ageaxis,midrange,'g',linewidth=3)
#        
#        if n in [1,5]: plt.ylabel('nactiv')
#        if n in [5,6,7,8]: plt.xlabel('age (day)')
#        # Calculate the slope of the curves during the second half and fit a curve
#        tc = 188
#        [total_slope,total_y0] = np.polyfit(ageaxis[tc:],np.log(total[tc:]),1)
#        [toprange_slope,toprange_y0] = np.polyfit(ageaxis[tc:],np.log(toprange[tc:]),1)
#        [midrange_slope,midrange_y0] = np.polyfit(ageaxis[tc:],np.log(midrange[tc:]),1)
#        [botrange_slope,botrange_y0] = np.polyfit(ageaxis[tc:],np.log(botrange[tc:]),1)
#        print('decay',supertype)
#        print('total   ',[total_slope,total_y0])
#        print('toprange',[toprange_slope,toprange_y0])
#        print('midrange',[midrange_slope,midrange_y0])
#        print('botrange',[botrange_slope,botrange_y0])
#        tc2 = 124
#        plt.semilogy(ageaxis[tc2:],np.exp(total_y0)*np.exp(total_slope*ageaxis[tc2:]),'k',
#                     ageaxis[tc2:],np.exp(midrange_y0)*np.exp(midrange_slope*ageaxis[tc2:]),'g',
#                     ageaxis[tc2:],np.exp(botrange_y0)*np.exp(botrange_slope*ageaxis[tc2:]),'b',
#                     linewidth=8,alpha=0.3)
#        if n in [1,2,5,6]: plt.ylim(1,2*10**6)
#        else:
#            plt.ylim(5*10**3,2*10**6)
#        plt.title('{} t:{:3.1f} m:{:3.1f} b:{:3.1f}'.format(supertype,-1/total_slope,-1/midrange_slope,-1/botrange_slope))
#        n += 1
#        
#    if figsave:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agenactiv-composit'))
#        
#    plt.show()

#%% 4ACP) Plot of the number of active parcels as a function of age in three separate 
# vertical bands of theta (and total)
# The rescaling factor is here due to converting from hour to day in tau, multiplying 
# by the integration factor Delta t = 6h and converting dedree in km
ff = (1/24) * (1/4) * degree**2
for hightype in hightypes:
    fig = plt.figure(figsize=(9,9))
    #fig.suptitle(hightype+' nactiv as a function of age  per range tot(k) top(r) bot(b) mid(g)')
    n = 1
    for supertype in supertypes:
        plt.subplot(2,2,n)
        # 65 is first thetaxis above 340 K and 95 is first thetaxis above 370 K
        total = np.sum(ff*result[hightype][supertype]['histog'][0:249,:],axis=1)
        toprange = np.sum(ff*result[hightype][supertype]['histog'][0:249,95:],axis=1)
        botrange = np.sum(ff*result[hightype][supertype]['histog'][0:249,0:65],axis=1)
        midrange = np.sum(ff*result[hightype][supertype]['histog'][0:249,65:95],axis=1)
        plt.semilogy(ageaxis,total,'k',ageaxis,toprange,'r',
                     ageaxis,botrange,'b',ageaxis,midrange,'g',linewidth=4)
        
        if n in [1,3]: plt.ylabel(r'cumulative impact (day$^2$ km$^2$)',fontsize=fs)
        if n in [3,4]: plt.xlabel(r'age (day)',fontsize=fs)
        plt.tick_params(labelsize=fst)
        plt.ylim(1000,3e8)
        # Calculate the slope of the curves during the second half and fit a curve
        tc = 140
        [total_slope,total_y0] = np.polyfit(ageaxis[tc:],np.log(total[tc:]),1)
        [toprange_slope,toprange_y0] = np.polyfit(ageaxis[tc:],np.log(toprange[tc:]),1)
        [midrange_slope,midrange_y0] = np.polyfit(ageaxis[tc:],np.log(midrange[tc:]),1)
        [botrange_slope,botrange_y0] = np.polyfit(ageaxis[tc:],np.log(botrange[tc:]),1)
        print('decay',supertype)
        print('total   ',[total_slope,total_y0])
        print('toprange',[toprange_slope,toprange_y0])
        print('midrange',[midrange_slope,midrange_y0])
        print('botrange',[botrange_slope,botrange_y0])
        tc2 = 124
        plt.semilogy(ageaxis[tc2:],np.exp(total_y0)*np.exp(total_slope*ageaxis[tc2:]),'k',
                     ageaxis[tc2:],np.exp(midrange_y0)*np.exp(midrange_slope*ageaxis[tc2:]),'g',
                     ageaxis[tc2:],np.exp(botrange_y0)*np.exp(botrange_slope*ageaxis[tc2:]),'b',
                     linewidth=12,alpha=0.3)
        #if n in [1,3]: plt.ylim(1,2*10**6)
        #else:
        #    plt.ylim(5*10**3,2*10**6)
        
        plt.title('{} t:{:3.1f} m:{:3.1f} b:{:3.1f}'.format(trea_short[supertype],-1/total_slope,-1/midrange_slope,-1/botrange_slope),fontsize=fs)
        n += 1
        
    if figpdf:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agenactiv-composit-withEAT-ACP.png'),bbox_inches='tight',dpi=dpi)
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agenactiv-composit-withEAT-ACP.pdf'),bbox_inches='tight',dpi=dpi)
    plt.show()
    
#%% 5) Plot of the vertical distribution of sources at (age 0 parcels)
# We only have data at 6h but that should not be too different in terms of vertical distribution
# By construction, the distribution of sources is the same for all runs
    # We use now the source distribution calculated directly from part_000 data in sources-stat.py

#with gzip.open(os.path.join(forw_dir,'source_dist.pkl'),'rb') as f:
#    source_dist = pickle.load(f)

#for hightype in hightypes:
#    fig = plt.figure(figsize=(14,7))
#    fig.suptitle(hightype+' vertical distribution of sources')
#    n = 1
#    for supertype in supertypes:
#        plt.subplot(2,4,n)
#        plt.semilogx(result[hightype][supertype]['histog'][0,:],thetaxis,linewidth=3)
#        plt.ylim(315,425)
#        if n in [1,5]: plt.ylabel('potential temperature (K)')
#        if n in [5,6,7,8]: plt.xlabel('souce density')
#        plt.title(supertype)
#        n += 1
    # define source_dist which is the same for all supertypes
    # notice that the source distribution is actually the distribution at 6h
    # to do: calculate the true source distribution from part_000 file
#    source_dist[hightype] = result[hightype][supertype]['histog'][0,:]
#    if figsave:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-vertsources-composit'))
#    plt.show()

#%% 6) Plot for each of the supertype of the normalized age histogram at 370 K, 380 K and 400 K ,
# that is positions 94, 104, 124
# show mean and modal peak
#for hightype in hightypes:
#    fig = plt.figure(figsize=(14,7))
#    fig.suptitle(hightype+' normalized age histogram at 370 k (b), 380 K (r) and 400 K (k)')
#    n = 1
#    for supertype in supertypes:
#        plt.subplot(2,4,n)
#        hh = result[hightype][supertype]['histog'][0:248,:]
#        ss = np.sum(hh,axis=0)
#        hh = hh / ss[np.newaxis,:]
#        agemean = np.sum(hh*ageaxis[:,np.newaxis],axis=0)
#        agemode = ageaxis[np.argmax(hh,axis=0)]
#        plt.plot(ageaxis,hh[:,94],'b',ageaxis,hh[:,104],'r',ageaxis,hh[:,124],'k',linewidth=3)
#        plt.scatter(agemean[94],0.015,c='b',marker='v',s=64)
#        plt.scatter(agemean[104],0.015,c='r',marker='v',s=64)
#        plt.scatter(agemean[124],0.015,c='k',marker='v',s=64)
#        plt.scatter(agemode[94],0.013,c='b',marker='D',s=64)
#        plt.scatter(agemode[104],0.013,c='r',marker='D',s=64)
#        plt.scatter(agemode[124],0.013,c='k',marker='D',s=64)
#        #plt.plot([agemean[94],],[0.015,],'bv',[agemean[104],],[0.015,],'rv',[agemean[124],],[0.015,],'kv',linewidth=12)
#        plt.title(supertype)
#        plt.ylim(0,0.016)
#        if n in [5,6,7,8]: plt.xlabel('age (day)')
#        n += 1
#    if figsave:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agehistog3levels-composit'))
#    plt.show() 

#%% 6ACP) Plot the supertypes EAZ, EAD, EIZ, EID of the normalized age histogram at 370 K, 380 K and 400 K ,
# that is positions 94, 104, 124
# show mean and modal peak
fs = 16
for hightype in hightypes:
    fig = plt.figure(figsize=(9,9))
    #fig.suptitle(hightype+' normalized age histogram at 370 k (b), 380 K (r) and 400 K (k)')
    n = 1
    for supertype in supertypes:
        plt.subplot(2,2,n)
        hh = result[hightype][supertype]['histog'][0:249,:].copy()
        ss = 0.25*np.sum(hh,axis=0)
        hh = hh / ss[np.newaxis,:]
        agemean = 0.25*np.sum(hh*ageaxis[:,np.newaxis],axis=0)
        agemode = ageaxis[np.argmax(hh,axis=0)]
        plt.plot(ageaxis,hh[:,94],'b',ageaxis,hh[:,104],'r',ageaxis,hh[:,124],'k',linewidth=5)
        plt.scatter(agemean[94],0.065,c='b',marker='v',s=128)
        plt.scatter(agemean[104],0.065,c='r',marker='v',s=128)
        plt.scatter(agemean[124],0.065,c='k',marker='v',s=128)
        plt.scatter(agemode[94],0.053,c='b',marker='D',s=96)
        plt.scatter(agemode[104],0.053,c='r',marker='D',s=96)
        plt.scatter(agemode[124],0.053,c='k',marker='D',s=96)
        #plt.plot([agemean[94],],[0.015,],'bv',[agemean[104],],[0.015,],'rv',[agemean[124],],[0.015,],'kv',linewidth=12)
        plt.title(trea[supertype],fontsize=fs)
        plt.tick_params(labelsize=fst)
        plt.ylim(0,0.07)
        if n in [1,3]: plt.ylabel('normalized age spectrum (day$^{-1}$)',fontsize=fs)
        if n in [3,4]: plt.xlabel('age (day)',fontsize=fs)
        n += 1
    if figpdf:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agehistog3levels-composit-withEAT-ACP.png'),bbox_inches='tight',dpi=dpi)
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agehistog3levels-composit-withEAT-ACP.pdf'),bbox_inches='tight',dpi=dpi)
    plt.show()  

#%% 7ACP) Comparison of the histograms for EAZ, EAD, EID, EIZ, EID-Return, EIZ-Return
# superimposed to the source curve    
# Set TeX to write the labels and annotations
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# normalization factor 1/24 * 1/4 (to convert units of previous time factor in
# day and to multiply per the time integration interval 6h and the inverse of the layer
# thickness (no factor for this later since layer thickness is 1 K)
ff = (1/24) * (1/4) * degree**2
# for the source distribution, the Delta t factor is not used
ff_s = (1/24)  * degree**2
fs = 22
for hightype in hightypes:
    fig,ax = plt.subplots(figsize=(9,9))
    ax.semilogx(ff*np.sum(result[hightype]['EID-Return']['histog'][0:249,:],axis=0),thetaxis,'--r',
                ff*np.sum(result[hightype]['EID']['histog'][0:249,:],axis=0),thetaxis,'r',
                ff*np.sum(result[hightype]['EAT']['histog'][0:249,:],axis=0),thetaxis,'--b',
                ff*np.sum(result[hightype]['EAD']['histog'][0:249,:],axis=0),thetaxis,'b',
                linewidth=6)
    ax.set_ylabel(r'potential temperature $\theta$ (K)',fontsize=fs)
    ax.set_xlabel(r'target cumulative impact distribution (day$^2$ km$^2$ K$^{-1}$)',fontsize=fs)
    ax.set_ylim(320,440)
    ax.set_xlim(1e5,5e8)
    ax.tick_params(labelsize=fs) 
    # superimpose the source curve in separate axis (same log range)
    ax2 = ax.twiny()
    ax2.set_xlabel(r'source distribution (day km$^2$ K$^{-1}$)',fontsize=fs,color='g')
    ax2.semilogx(ff_s*source_dist[hightype],thetaxis,'g',linewidth=8,alpha=0.5)
    ax2.tick_params(axis='x',labelcolor='g',labelsize=fs)
    ax2.set_xlim(1e4,5e8)
    ax.legend([r'\textit{ERA-I dia ret}',r'\textbf{ERA-I dia}',
               r'\textit{ERA5 }',r'\textbf{ERA5 dia}'],
               fontsize=fs,loc='upper right')
    # plt.title('Cumulative impact in the AMA region as a function of altitude',fontsize=fs)
    # plot horizontal line at modal max of sources
    ax2.plot([1e4,5e8],[thetaxis[np.argmax(source_dist[hightype])],
                    thetaxis[np.argmax(source_dist[hightype])]],'g')
    ax2.annotate(r'$\theta$ = 349.5 K',(1.50e4,351),fontsize=fs)
    print('theta for max source',thetaxis[np.argmax(source_dist[hightype])])
    if figpdf:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-age-compar-impact-profile-composit-withEAT-ACP.png'),bbox_inches='tight',dpi=dpi)
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-age-compar-impact-profile-composit-withEAT-ACP.pdf'),bbox_inches='tight',dpi=dpi)
    plt.show()      
        
##%% 7bisACP) Plot of the proportion of parcels above a given level relative to the same proportion in the sources
## Calculation of the proportion for the sources
#fs = 22  
#
#for hightype in hightypes:
#    source_above = np.flip(np.cumsum(np.flip(source_dist[hightype]))/np.sum(source_dist[hightype]))
#    source_below = np.cumsum(source_dist[hightype])/np.sum(source_dist[hightype])
#    target_above = {}
#    target_below = {}
#    fig,ax = plt.subplots(figsize=(9,9))
#    for supertype in supertypes:
#        target_sum = np.sum(result[hightype][supertype]['histog'][0:249,:])
#        target_above[supertype] = np.flip(np.cumsum(np.flip(np.sum(result[hightype][supertype]['histog'][0:249,:],axis=0))))\
#            / (target_sum * source_above)
#        target_below[supertype] = np.cumsum(np.sum(result[hightype][supertype]['histog'][0:249,:],axis=0))\
#            / (target_sum * source_below)
#    ax.plot(target_above['EIZ']*target_below['EAZ'],thetaxis,'--r',
#            target_above['EID']*target_below['EAD'],thetaxis,'r',
#            target_above['EAZ']*target_below['EIZ'],thetaxis,'--b',
#            target_above['EAD']*target_below['EID'],thetaxis,'b',
#            linewidth=6)
#    # plot a vertical line for the unit ratio
#    ax.plot(np.ones(len(thetaxis)),thetaxis,'k')
#    ax.set_ylabel('target potential temperature (K)',fontsize=fs)
#    ax.set_xlabel('target/source cumulative ratio distribution',fontsize=fs)
#    ax.set_ylim(320,440)
#    ax.set_xlim(0,14)
#    ax.tick_params(labelsize=fs)
#    # superimpose the source curve
#    ax2 = ax.twiny()
#    ax2.semilogx(ff_s*source_dist[hightype],thetaxis,'g',linewidth=8,alpha=0.5)
#    ax2.tick_params(axis='x',labelcolor='g',labelsize=fs)
#    ax2.set_xlabel('source distribution (day km$^2$ K$^{-1}$)',fontsize=fs,color='g')
#    ax2.set_xlim(1e4,5e8)
#    ax2.plot([1e4,5e8],[thetaxis[np.argmax(source_dist[hightype])],
#                    thetaxis[np.argmax(source_dist[hightype])]],'g')
#    ax2.annotate(r'$\theta$ = 349.5 K',(15e4,351),fontsize=fs)
#    ax.legend([r'\textbf{ERA5 kin}',r'\textit{ERA5 dia}',
#                r'\textbf{ERA-I kin}',r'\textit{ERA-I dia}'],
#                fontsize=fs,loc='upper right')
#    #plt.title(hightype+' cumulative impact ratio as a function of age',fontsize=fs)
#    # plot horizontal line at modal max of sources
#    #ax.plot([0,14],[thetaxis[np.argmax(source_dist[hightype])],
#    #                thetaxis[np.argmax(source_dist[hightype])]],'g')
#    if figpdf:
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-age-compar-impact-ratio-composit-ACP.png'),bbox_inches='tight',dpi=dpi)
#        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-age-compar-impact-ratio-composit-ACP.pdf'),bbox_inches='tight',dpi=dpi)  
#    plt.show()
#    
##%% analysis of diffusion
#alpha = 1/13.3
#AA = {'EAD':1.07,'EAZ':1.11,'EIZ':0.98,'EID':1.35}
#for supertype in supertypes:
#    # analysis as a function of the vertical coordinate theta by summing over time
#    hht = result['sh'][supertype]['histog'][0:249,74:140].copy()
#    thets = thetaxis[74:140]
#    # sum over time
#    sst = np.sum(hht,axis=0)
#    # max along the time axis 
#    hhtmax = np.max(hht,axis=0)
#    # normalizing each histogram for a given theta
#    hht /= sst[np.newaxis,:]
#    fig = plt.figure(figsize=(12,6))
#    plt.subplot(2,4,1)
#    # mean age for theta values
#    agemean = np.sum(hht*ageaxis[:,np.newaxis],axis=0)
#    agemode = ageaxis[np.argmax(hht,axis=0)]
#    plt.plot(agemean,thets,agemode,thets)
#    plt.ylim(370,400)
#    plt.title(supertype+' mean and mode (t)')
#    plt.subplot(2,4,2)
#    var2age = np.sum(hht*ageaxis[:,np.newaxis]**2,axis=0) - agemean**2
#    [b,a] = np.polyfit(agemean[10:30],var2age[10:30],1)
#    dcor = (0.5*AA[supertype]**2 * b)/(1 - 2*alpha*b)
#    print('diff1 '+supertype,0.5*b,dcor)
#    plt.plot(agemean,var2age,agemean,a + b*agemean)
#    plt.title(supertype+' age: var vs mean' )
#    plt.subplot(2,4,3)
#    tags=[5,10,15,20,25,30,35,40]
#    for tag in tags:
#        plt.plot(ageaxis,hht[:,tag])
#    plt.subplot(2,4,4)
#    plt.plot(hhtmax[20:50]*np.exp(alpha*agemode[20:50])*np.sqrt(agemode[20:50]),thets[20:50])
#   
#    # analysis as a function of time by summing over theta
#    # this method does not work possibly for lack of separation in kinematic
#    hhp = result['sh'][supertype]['histog'][0:249,74:140].copy()
#    thets = thetaxis[74:140]
#    ssp = np.sum(hhp,axis=1)
#    hhpmax = np.max(hhp,axis=1)
#    hhp /= ssp[:,np.newaxis]
#    plt.subplot(2,4,5)
#    thetmean = np.sum(hhp*thets[np.newaxis,:],axis=1)
#    thetmode = thets[np.argmax(hh,axis=1)]
#    plt.plot(ageaxis,thetmean,ageaxis,thetmode)
#    plt.title(supertype+' mean and mode (thet)')
#    plt.subplot(2,4,6)
#    var2thet = np.sum(hhp*thets[np.newaxis,:]**2,axis=1) - thetmean**2
#    [b,a] = np.polyfit(ageaxis[10:60],var2thet[10:60],1)
#    print('diff2 '+supertype,0.5*b)
#    plt.plot(ageaxis,var2thet,ageaxis,a + b*ageaxis)
#    plt.subplot(2,4,7)
#    tags=[5,10,20,30,40,50,60,70,80,100,140,160,180,200,220,240]
#    for tag in tags:
#        plt.plot(hhp[tag,:],thets)
#    plt.subplot(2,4,8)
#    print('hhpmax0 '+supertype,hhpmax[0])
#    plt.plot(ageaxis,hhpmax*np.exp(alpha*ageaxis)*np.sqrt(ageaxis))
#    plt.show()    