#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script provides statistics of the forward runs of the FORWBox-meanhigh family

It is based on the output of ageStat.py and is a variant of plotAgeStat
In plotAgeStat, the diagnostics are shown for all streams for a given supertype

Here we show the results for the whole summer 2017 (summing the streams) 
and as a comparison between the supertypes in the FullAMA domain 
(EAZ, EAD, EIZ, EID, EIZ-Return and EID-Return)

It generates for both sh and mh hightypes 4  figures which are
(actually mh not plotted)

This script is the basis for the more specialized versions used later
(ACP and CONF2019)

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

EAZ, EAD, EIZ, EID are showing statistics of parcels in the FullAMA domain with the rule that a particle that exits
once is discarded.
EIZ-FULL, EID-FULL show the statistics of parcels in the global domain
EIZ-Return, EID-Return show the statistics of parcels in the FullAMA domain by applying a simple mask to the FULL runs, 
that is parcels that leave the domain are counted if they return.  

ACHTUNG: the normalizations need to be corrected according to the ACP version

Created on Sat 23 Feb 2019

@author: Bernard Legras
"""
import numpy as np
#from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
import argparse
import socket
import os

parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",type=str,choices=["EAD","EAZ","EIZ","EID","EIZ-Return","EID-Return","EIZ-FULL","EID-FULL"],help="type")
#parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')

supertypes = ["EAZ","EIZ","EIZ-Return","EIZ-FULL","EAD","EID","EID-Return","EID-FULL"]
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

args = parser.parse_args()
if args.type is not None: supertype = args.type

dates = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]

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
        print('reading',file_out )
        with gzip.open(file_out,'rb') as f:
           [_,histog,_,_] = pickle.load(f)
        for hightype in hightypes:
            result[hightype][supertype]['histog'] += histog[hightype]
    for hightype in hightypes:
        result[hightype][supertype]['nactiv'] = np.sum(result[hightype][supertype]['histog'],axis=1)

#%% 1) Plot of the age/thet histogram
for hightype in hightypes:
    fig = plt.figure(figsize=(14,7))
    fig.suptitle(hightype+' age theta histogram' )
    n = 1
    for supertype in supertypes:
        plt.subplot(2,4,n)
        plt.imshow(np.log10(result[hightype][supertype]['histog'][0:248,:]).T,
                   extent=(0.25,62,275,700),origin='lower',aspect='auto',cmap='jet',clim=(-1,5))
        plt.ylim(310,450)
        plt.colorbar()
        plt.title(supertype)      
        if n in [1,5]: plt.ylabel('potential temperature (K)')
        if n in [5,6,7,8]: plt.xlabel('age (day)')
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethethist-composit'))
    plt.show()
    
#%% 2) Plot of the age/thet histogram, normalized for each level
for hightype in hightypes:
    fig = plt.figure(figsize=(14,7))
    fig.suptitle(hightype+' age theta normalized histogram' )
    n = 1
    for supertype in supertypes:
        plt.subplot(2,4,n)
        hh = result[hightype][supertype]['histog'][0:248,:]
        ss = np.sum(hh,axis=0)
        hh = hh / ss[np.newaxis,:]
        plt.imshow(np.log10(hh.T),
                   extent=(0.25,62,275,700),origin='lower',aspect='auto',cmap='jet',clim=(-6,0))
        plt.ylim(310,450)
        plt.colorbar()
        plt.title(supertype)      
        if n in [1,5]: plt.ylabel('potential temperature (K)')
        if n in [5,6,7,8]: plt.xlabel('age (day)')
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethetnormhist-composit'))
    plt.show()
    
#%% 3) Mean age
ageaxis = np.arange(0.25,62.25,0.25) 
thetaxis = np.arange(275.5,700)
for hightype in hightypes:
    fig = plt.figure(figsize=(14,7))
    fig.suptitle(hightype+'  mean and modal age ')
    n = 1
    for supertype in supertypes:
        plt.subplot(2,4,n)
        hh = result[hightype][supertype]['histog'][0:248,:]
        ss = np.sum(hh,axis=0)
        ss[ss==0]=1
        hh = hh / ss[np.newaxis,:]
        agemean = np.sum(hh*ageaxis[:,np.newaxis],axis=0)
        agemode = ageaxis[np.argmax(hh,axis=0)]
        plt.plot(agemean,thetaxis,'k',agemode,thetaxis,'r')
        plt.ylim(325,425)
        plt.title(supertype)
        if n in [1,5]: plt.ylabel('potential temperature (K)')
        if n in [5,6,7,8]: plt.xlabel('age (day)')
        n +=1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agemean-composit'))
    plt.show()

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
for hightype in hightypes:
    fig = plt.figure(figsize=(14,7))
    fig.suptitle(hightype+' nactiv as a function of age  per range tot(k) top(r) bot(b) mid(g)')
    n = 1
    for supertype in supertypes:
        plt.subplot(2,4,n)
        # 65 is first thetaxis above 340 K and 95 is first thetaxis above 370 K
        total = np.sum(result[hightype][supertype]['histog'][0:248,:],axis=1)
        toprange = np.sum(result[hightype][supertype]['histog'][0:248,95:],axis=1)
        botrange = np.sum(result[hightype][supertype]['histog'][0:248,0:65],axis=1)
        midrange = np.sum(result[hightype][supertype]['histog'][0:248,65:95],axis=1)
        plt.semilogy(ageaxis,total,'k',ageaxis,toprange,'r',
                     ageaxis,botrange,'b',ageaxis,midrange,'g',linewidth=3)
        
        if n in [1,5]: plt.ylabel('nactiv')
        if n in [5,6,7,8]: plt.xlabel('age (day)')
        # Calculate the slope of the curves during the second half and fit a curve
        tc = 188
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
                     linewidth=8,alpha=0.3)
        if n in [1,2,5,6]: plt.ylim(1,2*10**6)
        else:
            plt.ylim(5*10**3,2*10**6)
        plt.title('{} t:{:3.1f} m:{:3.1f} b:{:3.1f}'.format(supertype,-1/total_slope,-1/midrange_slope,-1/botrange_slope))
        n += 1
        
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agenactiv-composit'))
        
    plt.show()
    
#%% 5) Plot of the vertical distribution of sources at (age 0 parcels)
# We only have data at 6h but that should not be too different in terms of vertical distribution
# By construction, the distribution of sources is the same for all runs
source_dist = {}
for hightype in hightypes:
    fig = plt.figure(figsize=(14,7))
    fig.suptitle(hightype+' vertical distribution of sources')
    n = 1
    for supertype in supertypes:
        plt.subplot(2,4,n)
        plt.semilogx(result[hightype][supertype]['histog'][0,:],thetaxis,linewidth=3)
        plt.ylim(315,425)
        if n in [1,5]: plt.ylabel('potential temperature (K)')
        if n in [5,6,7,8]: plt.xlabel('souce density')
        plt.title(supertype)
        n += 1
    # define source_dist which is the same for all supertypes
    source_dist[hightype] = result[hightype][supertype]['histog'][0,:]
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-vertsources-composit'))
    plt.show()

#%% 6) Plot for each of the supertype of the normalized age histogram at 370 K, 380 K and 400 K ,
# that is positions 94, 104, 124
# show mean and modal peak
for hightype in hightypes:
    fig = plt.figure(figsize=(14,7))
    fig.suptitle(hightype+' normalized age histogram at 370 k (b), 380 K (r) and 400 K (k)')
    n = 1
    for supertype in supertypes:
        plt.subplot(2,4,n)
        hh = result[hightype][supertype]['histog'][0:248,:]
        ss = np.sum(hh,axis=0)
        hh = hh / ss[np.newaxis,:]
        agemean = np.sum(hh*ageaxis[:,np.newaxis],axis=0)
        agemode = ageaxis[np.argmax(hh,axis=0)]
        plt.plot(ageaxis,hh[:,94],'b',ageaxis,hh[:,104],'r',ageaxis,hh[:,124],'k',linewidth=3)
        plt.scatter(agemean[94],0.015,c='b',marker='v',s=64)
        plt.scatter(agemean[104],0.015,c='r',marker='v',s=64)
        plt.scatter(agemean[124],0.015,c='k',marker='v',s=64)
        plt.scatter(agemode[94],0.013,c='b',marker='D',s=64)
        plt.scatter(agemode[104],0.013,c='r',marker='D',s=64)
        plt.scatter(agemode[124],0.013,c='k',marker='D',s=64)
        #plt.plot([agemean[94],],[0.015,],'bv',[agemean[104],],[0.015,],'rv',[agemean[124],],[0.015,],'kv',linewidth=12)
        plt.title(supertype)
        plt.ylim(0,0.016)
        if n in [5,6,7,8]: plt.xlabel('age (day)')
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agehistog3levels-composit'))
    plt.show()  

#%% 7) Comparison of the histograms for EAZ, EAD, EID, EIZ, EID-Return, EIZ-Return
# superimposed to the source curve    
# Set TeX to write the labels and annotations
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# normalization factor 1/24 * 1/4 (to convert units of previous time factor in
# day and to multiply per the time integration interval 6h)
ff = (1/24) * (1/4)

for hightype in hightypes:
    fig,ax = plt.subplots(figsize=(8,8))
    fs = 18
    ax.semilogx(ff*np.sum(result[hightype]['EAZ']['histog'][0:248,:],axis=0),thetaxis,'b',
                ff*np.sum(result[hightype]['EAD']['histog'][0:248,:],axis=0),thetaxis,'-.b',
                ff*np.sum(result[hightype]['EIZ']['histog'][0:248,:],axis=0),thetaxis,'r',
                ff*np.sum(result[hightype]['EID']['histog'][0:248,:],axis=0),thetaxis,'-.r',
                ff* np.sum(result[hightype]['EIZ-Return']['histog'][0:248,:],axis=0),thetaxis,'k',
                ff* np.sum(result[hightype]['EID-Return']['histog'][0:248,:],axis=0),thetaxis,'-.k',
                 linewidth=4)
    ax.set_ylabel('target potential temperature (K)',fontsize=fs)
    ax.set_xlabel('target impact distribution (day$^2$ degree$^2$ K$^{-1}$)',fontsize=fs)
    ax.set_ylim(320,420)
    ax.set_xlim(5*10**1,5*10**4)
    ax.tick_params(labelsize=fs) 
    # superimpose the source curve in separate axis (same log range)
    ax2 = ax.twiny()
    ax2.set_xlabel('source distribution (day degree$^2$ K$^{-1}$)',fontsize=fs,color='g')
    ax2.semilogx(source_dist[hightype],thetaxis,'g',linewidth=8,alpha=0.5)
    ax2.tick_params(axis='x',labelcolor='g',labelsize=fs)
    ax2.set_xlim(10**2,10**6)
    ax.legend([r'\textbf{ERA5 kinematic AMA}',r'\textit{ERA5 diabatic AMA}',
                r'\textbf{ERA-I kinematic AMA}',r'\textit{ERA-I diabatic AMA}',
                r'\textbf{ERA-I kinematic Return}',r'\textit{ERA-I diabatic Return}'],
                fontsize=fs,loc='center left')
    plt.title('Total '+hightype+' impact in the FULLAMA region as a function of altitude',fontsize=fs)
    # plot horizontal line at modal max of sources
    ax2.plot([10**2,10**6],[thetaxis[np.argmax(source_dist[hightype])],
                    thetaxis[np.argmax(source_dist[hightype])]],'g')
    print('theta for max source',thetaxis[np.argmax(source_dist[hightype])])
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-age-compar-impact-profile-composit'))
    plt.show()      
        
#%% 7bis) Plot of the proportion of parcels above a given level relative to the same proportion in the sources
# Calculation of the proportion for the sources
    
source_above = np.flip(np.cumsum(np.flip(source_dist['sh']))/np.sum(source_dist['sh']))
source_below = np.cumsum(source_dist['sh'])/np.sum(source_dist['sh'])
target_above = {}
target_below = {}
for hightype in hightypes:
    fig,ax = plt.subplots(figsize=(8,8))
    fs = 18
    for supertype in supertypes:
        target_sum = np.sum(result[hightype][supertype]['histog'][0:248,:])
        target_above[supertype] = np.flip(np.cumsum(np.flip(np.sum(result[hightype][supertype]['histog'][0:248,:],axis=0))))\
            / (target_sum * source_above)
        target_below[supertype] = np.cumsum(np.sum(result[hightype][supertype]['histog'][0:248,:],axis=0))\
            / (target_sum * source_below)
    ax.plot(target_above['EAZ']*target_below['EAZ'],thetaxis,'b',
            target_above['EAD']*target_below['EAD'],thetaxis,'-.b',
            target_above['EIZ']*target_below['EIZ'],thetaxis,'r',
            target_above['EID']*target_below['EID'],thetaxis,'-.r',
            target_above['EIZ-Return']*target_below['EIZ-Return'],thetaxis,'k',
            target_above['EID-Return']*target_below['EID-Return'],thetaxis,'-.k',
            linewidth=4)
    # plot a vertical line for the unit ratio
    ax.plot(np.ones(len(thetaxis)),thetaxis,'k')
    ax.set_ylabel('target potential temperature (K)',fontsize=fs)
    ax.set_xlabel('target/source cumul ratio distribution',fontsize=fs)
    ax.set_ylim(320,420)
    ax.set_xlim(0,14)
    ax.tick_params(labelsize=fs)
    # superimpose the source curve
    ax2 = ax.twiny()
    ax2.semilogx(source_dist[hightype],thetaxis,'g',linewidth=8,alpha=0.5)
    ax2.tick_params(axis='x',labelcolor='g',labelsize=fs)
    ax2.set_xlabel('source distribution',fontsize=fs,color='g')
    ax2.set_xlim(10**2,5*10**5)
    ax.legend([r'\textbf{ERA5 kinematic AMA}',r'\textit{ERA5 diabatic AMA}',
                r'\textbf{ERA-I kinematic AMA}',r'\textit{ERA-I diabatic AMA}',
                r'\textbf{ERA-I kinematic Return}',r'\textit{ERA-I diabatic Return}'],
                fontsize=fs,loc='upper right')
    plt.title(hightype+' cumulative impact ratio as a function of age',fontsize=fs)
    # plot horizontal line at modal max of sources
    ax.plot([0,14],[thetaxis[np.argmax(source_dist[hightype])],
                    thetaxis[np.argmax(source_dist[hightype])]],'g')
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-age-compar-impact-ratio-composit'))
    plt.show() 