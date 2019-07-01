# -*- coding: utf-8 -*-
"""
Analysis of the forward runs to determine the proportion of parcels going up with 
respect to that going down.
A normalization factor is required to reduce this to impact as in ageStat
This code filters the wrong image of 30 August at 11h and all the sources above 
360K at latitudes >= 40N

The code establishes regional statistics in hdisptheta
It also collects the raw statistics in hdisptheta_raw but cannot blacklist the 
spurious data here because the regional histograms are made directly in 
history-forw-uniq (to be corrected for that).

hdisptheta is a 3D histogram in (lev,lev0,reg) where lev is the mean theta level
of the parcel during its trajectory, lev0 is the source level and reg is the region
of the group of regions to which belong the parcel. This is calculated after 
the properties determined in history.
The histogram is both established for sh and mh hightype

The histogram is then used to diagnose for each level and region in the source
domain, the proportion of parcels remaining at the same level (5K vertical bin),
the proportion of parcels going above, and those going below

In this diagnostic, each parcel has the same weight, whatever its life time.
 
A second histogram is established from the raw regional histogram calculated 
in history (only for mh so far). This histogram uses all the positions occupied 
by a parcel along its trajectory and is therefore weighted by the duration of 
each trajectory.

The analysis produces an output file diag*** that contains the calculated 
histograms, and the diagnostics of vertical displacement. This is done for each
valid couple of supertype + target, that is : EAD, EAZ, EID, EIZ in the FullAMA
domain and EID, EIZ in the global domain.

The results are plotted by diag-plot-2.py

Created on Tue Apr 30 22:52:48 2019

@author: Bernard Legras
"""
import numpy as np
#from datetime import datetime,timedelta
import pickle,gzip
#import matplotlib.pyplot as plt
from os.path import join
import socket
import constants as cst
from numba import njit

#%%
@njit
def accumul(hist,lev,lev0,regs,factor):
    for i in range(len(lev)):
        hist[lev[i],lev0[i],regs[i]] += factor[i]

#%%
# parameters
dates = ['Jun-01','Jun-11','Jun-21','Jul-01','Jul-11','Jul-21','Aug-01','Aug-11','Aug-21']
vert = 'theta'
age_max = 62
age_max_inter = 30
# family: association of supertype with target parameter
family = {1:['EID-FULL','global'],2:['EIZ-FULL','global'],3:['EAD','FullAMA'],4:['EAZ','FullAMA'], \
          5:['EID-FULL','FullAMA'],6:['EIZ-FULL','FullAMA']}

# read mask an initialize edges at 0.25° resolution
print ('open mask')
with gzip.open(join('..','mkSTCmask','MaskCartopy2-STCforwfine.pkl'),'rb') as f: 
    mm =pickle.load(f) 
xedge = np.arange(-10,160+0.5*mm['icx'],mm['icx'])
yedge = np.arange(0,50+0.5*mm['icy'],mm['icx'])
# 20 bins àf 5K in the vertical direction, centered on [325, 330, ..., 420]
minv = 322.5
maxv = 422.5
binv = 20
binvl = 5
vcent = np.arange(325,421,(maxv-minv)/binv)

# fixed derived parameters
hmax = 24*(age_max + 10)
hinter = 24*(age_max_inter + 10)

# normalizing factor (1h *  (0.1degree)^2  / 5K)
ff0 = 0.01 / 5
 
# input directory
if 'ciclad' in socket.gethostname():
    in_dir = '/data/legras/STC/STC-FORWBox-meanhigh-OUT'
elif  socket.gethostname() in ['gort','satie']:
    in_dir = "/home/legras/data/STC/STC-FORWBox-meanhigh-OUT"

hightypes = ['mh','sh']

for ff in family:
    supertype = family[ff][0]
    target = family[ff][1]
    print(supertype,target)
    # derived parameters
    run_type = supertype+'-Box-'+vert+'-'+target
   
    # initializes histograms
    # hdisptheta is based on the mean theta
    # hdisptheta_raw is based on the raw histogram that has been constructed
    # in history using all the location of the parcel along its trajectory
    # it has only be calculated for mh
    hdisptheta = {} 
    hdisptheta['mh'] = np.zeros(shape=(binv,binv,len(mm['regcode'])+1),dtype=np.float32)
    hdisptheta['sh'] = np.zeros(shape=(binv,binv,len(mm['regcode'])+1),dtype=np.float32)
    # to be updated when sh is processed in history
    hdisptheta_raw = {} 
    hdisptheta_raw['mh'] = np.zeros(shape=(binv,binv,len(mm['regcode'])+1),dtype=np.float32)
    #hdisptheta_raw['sh'] = np.zeros(shape=(binv,binv,len(mm['regcode'])+1),dtype=np.int)
    
    # loop on dates
    for date in dates:
        # data file is read from the product of history
        history_stream = join(in_dir,'history-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')
        with gzip.open(history_stream,'rb') as f: 
            sources = pickle.load(f)
        print('file has been read',date)
        # accumulate the raw histogram
        # to be updated when sh is processed in history
        hdisptheta_raw['mh'] += sources['hdisptheta'] * ff0
        # statistics for mh
        selec = sources['age']>0
        # filtering out 30 August at 11h
        if dates == 'Aug-21':
            selec = selec & (sources['ir_start'] != 227)
        # filtering out all souces above 360K at latitudes larger than 360K
        thet = sources['t']*(cst.p0/sources['p'])**cst.kappa
        highhigh = (thet>360) & (sources['y']>=40)
        selec = selec &  ~highhigh
        del highhigh
        del thet
        # statistics for mh 
        # selecting the sources in the regonal grid
        ix0 = np.floor((sources['x'][selec]-xedge[0])/mm['icx']).astype(int)
        jy0 = np.floor((sources['y'][selec]-yedge[0])/mm['icy']).astype(int)
        ix0 = np.clip(ix0,0,mm['nlons']-1)
        jy0 = np.clip(jy0,0,mm['nlats']-1)
        # full scaling factor including the cosine of the source latitude
        factor = ff0 * np.cos(np.deg2rad(sources['y'][selec]))
        # attribution of the regional code from the mask
        regs =( mm['mask'][jy0,ix0]).astype(int)
        theta = sources['t'][selec]*(cst.p0/ sources['p'][selec])**cst.kappa
        theta_mean = sources['theta_mean'][selec]/sources['count'][selec]
        # selecting the level of the source and of the theta_mean
        lev0 = np.floor((theta-minv)/binvl).astype(int)
        lev0 = np.clip(lev0,0,binv-1)
        lev =  np.floor((theta_mean-minv)/binvl).astype(int)
        lev = np.clip(lev,0,binv-1)
        # accumulate 
        accumul(hdisptheta['mh'],lev,lev0,regs,factor)
        # statistics for sh repeating the previous sequence until accumulation
        selec = selec & sources['silviahigh']
        ix0 = np.floor((sources['x'][selec]-xedge[0])/mm['icx']).astype(int)
        jy0 = np.floor((sources['y'][selec]-yedge[0])/mm['icy']).astype(int)
        ix0 = np.clip(ix0,0,mm['nlons']-1)
        jy0 = np.clip(jy0,0,mm['nlats']-1)
        factor = ff0 * np.cos(np.deg2rad(sources['y'][selec]))
        regs =( mm['mask'][jy0,ix0]).astype(int)
        theta = sources['t'][selec]*(cst.p0/ sources['p'][selec])**cst.kappa
        theta_mean = sources['theta_mean'][selec]/sources['count'][selec]
        lev0 = np.floor((theta-minv)/binvl).astype(int)
        lev0 = np.clip(lev0,0,binv-1)
        lev =  np.floor((theta_mean-minv)/binvl).astype(int)
        lev = np.clip(lev,0,binv-1)
        accumul(hdisptheta['sh'],lev,lev0,regs,factor)
        
    #%% diagnostic of the proportions for histograms of the mean theta
    # here each trajectory counts for one unit
    # diagnose for each level and region, using impact metric, the trajectories 
    # - that have a theta mean in the same bin
    # - the trajectories that have a theta mean above
    # - the trajectories that have a theta mean below
    # and calculate ratios to total
    diag = {}
    for hightype in hightypes:
        diag[hightype] = {}
        for var in ['equal','above','below','total']:
            diag[hightype][var] = np.zeros(shape=(binv,len(mm['regcode'])+1),dtype=np.float32)
        for ireg in range(1,len(mm['regcode'])+1):
            for lev in range(binv):
                diag[hightype]['equal'][lev,ireg] = hdisptheta[hightype][lev,lev,ireg]
                diag[hightype]['above'][lev,ireg] = np.sum(hdisptheta[hightype][lev+1:,lev,ireg])
                diag[hightype]['below'][lev,ireg] = np.sum(hdisptheta[hightype][:lev,lev,ireg])
                diag[hightype]['total'][lev,ireg] = np.sum(hdisptheta[hightype][:,lev,ireg])
        ddt = diag[hightype]['total'].copy()
        ddt[ddt==0] = 1
        diag[hightype]['equal_prop'] = diag[hightype]['equal']/ddt
        diag[hightype]['above_prop'] = diag[hightype]['above']/ddt
        diag[hightype]['below_prop'] = diag[hightype]['below']/ddt        
    
    #%% Grouping of the regions into blocks and cumulating the previous dignostics
    # into these blocks
    group = {}
    group['Land']  = ['IndianSub','SouthChina','Pen','Pakistan','Bangladesh']
    group['Seas'] = ['BoB','SCSPhi']
    group['Ocean'] = group['Seas'] + ['IndianOcean','Indonesia','WestPacific']
    group['Tibet'] = ['TibetanPlateau',]
    for hightype in hightypes:
        for gr in group.keys():
            diag[hightype][gr] = {}
            # cumulating the diagnostics
            for var in ['equal','above','below','total']:
                diag[hightype][gr][var] = np.zeros(binv,dtype=np.float32)
                for reg in group[gr]:
                    diag[hightype][gr][var] += diag[hightype][var][:,mm['regcode'][reg]]
            ddt = diag[hightype][gr]['total'].copy()
            ddt[ddt==0] = 1
            diag[hightype][gr]['equal_prop'] = diag[hightype][gr]['equal']/ddt
            diag[hightype][gr]['above_prop'] = diag[hightype][gr]['above']/ddt
            diag[hightype][gr]['below_prop'] = diag[hightype][gr]['below']/ddt
            del ddt
        
    #%% diagnostic of the proportions for raw histograms
    # no weighting according to the latitude here and no processing of silviahigh for the moment (therefore quite useless)
    # provision for that processing is made for that in diag_raw, however 
    # Here each trajectory is weighted according to its residence time
    diag_raw = {}
    for hightype in ['mh',]:
        diag_raw[hightype] = {}
        for var in ['equal','above','below','total']:
            diag_raw[hightype][var] = np.zeros(shape=(binv,len(mm['regcode'])+1),dtype=np.float32)
        for ireg in range(1,len(mm['regcode'])+1):
            for lev in range(binv):
                diag_raw[hightype]['equal'][lev,ireg] = hdisptheta_raw[hightype][lev,lev,ireg]
                diag_raw[hightype]['above'][lev,ireg] = np.sum(hdisptheta_raw[hightype][lev+1:,lev,ireg])
                diag_raw[hightype]['below'][lev,ireg] = np.sum(hdisptheta_raw[hightype][:lev,lev,ireg])
                diag_raw[hightype]['total'][lev,ireg] = np.sum(hdisptheta_raw[hightype][:,lev,ireg])
        ddt = diag_raw[hightype]['total'].copy()
        ddt[ddt==0] = 1
        diag_raw[hightype]['equal_prop'] = diag_raw[hightype]['equal']/ddt
        diag_raw[hightype]['above_prop'] = diag_raw[hightype]['above']/ddt
        diag_raw[hightype]['below_prop'] = diag_raw[hightype]['below']/ddt        
    
    #%% Grouping of the regions into blocs 
        for gr in group.keys():
            diag_raw[hightype][gr] = {}
            for var in ['equal','above','below','total']:
                diag_raw[hightype][gr][var] = np.zeros(binv,dtype=np.float32)
                for reg in group[gr]:
                    diag_raw[hightype][gr][var] += diag_raw[hightype][var][:,mm['regcode'][reg]]
            ddt = diag_raw[hightype][gr]['total'].copy()
            ddt[ddt==0] = 1
            diag_raw[hightype][gr]['equal_prop'] = diag_raw[hightype][gr]['equal']/ddt
            diag_raw[hightype][gr]['above_prop'] = diag_raw[hightype][gr]['above']/ddt
            diag_raw[hightype][gr]['below_prop'] = diag_raw[hightype][gr]['below']/ddt
            del ddt
            
    # %% output the result for each family (supertype + target)
    with gzip.open('diag-'+supertype+'-'+target+'-'+vert+'.pkl','wb') as f:
        pickle.dump([[diag,hdisptheta],[diag_raw,hdisptheta_raw]],f)