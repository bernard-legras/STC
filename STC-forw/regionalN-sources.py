# -*- coding: utf-8 -*-
"""
Created on 30 March 2018

This code processes the forward runs to provide transit statistics

It calculate transit properties by accumulating for all the parcels during their 
lifetime under some filtering conditions given below.
The two major filters are based on the age and on the cloud type associated to 
the parcel.

Two values of the age max are managed simultaneously generating two transit 
tables. They are given in age_max and in age_max_inter

The veryhigh and silviahigh sources are flagged before calling transit.
The forward runs are therefore expected to have been initialized with sources
that include these two types of sources, like meanhigh.

This version should work for all cases

Generate printed diagnostics and two output pickle files.
The two output pickle files are generated for the two values of the max age. They
are distinguished by the suffix hxxx where xxx is the max age in hours. 
Take input into runs archived in /data/legras/flexout/STC/FORWBox-meanhigh    
and generate outputs in /data/legras/STC/STC-FORWBox-meanhigh-OUT

Parameters:
--type: supertype of wind data used and domain used in this analysis
"EAD","EAZ","EIZ","EID","EIZ-FULL","EID-FULL"
--vert: choice of the type of vertical levels used to determine the transit properties
"baro" or "theta"
--quiet: determines whether printed outputs go to a file or to the screen
--full: determines whether the domain is global, can only be used with EID-FULL and EIZ-FULL
--date: determines the run stream following the segmentation of the 
runs
["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]
--age_max: age max (in days)
--age_max_inter: second age max, should be less than the first
--inc:  step inc for which the part file are processed    

Parameterized for Box-meanhigh on 15 Jan 2019
Forked from stat-forw to process a unique run
Needs a subsequent script stat-forw-gather to gather the data

The parcels launched from the erroneous image on 30 August at 11 are eliminated.

@author: Bernard Legras
"""

from os.path import join
import numpy as np
#import matplotlib.pyplot as plt
#from datetime import datetime, timedelta
import gzip, pickle
import sys
import argparse
import socket
#from subprocess import call
import io107
import constants as cst

#%%
parser = argparse.ArgumentParser()
#parser.add_argument("-y","--year",type=int,help="year")
#parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
parser.add_argument("-t","--type",choices=["EAT","EAD","EAZ","EIZ","EID","EIZ-FULL","EID-FULL"],help="type")
#parser.add_argument("-s","--saf",choices=["O","N"],help="SAF version")
parser.add_argument("-v","--vert",choices=["theta","baro"],help="vertical discretization")
parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
parser.add_argument("-f","--full",type=str,choices=["y","n"],help="full (y) or not (n)")
parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')
parser.add_argument("-am","--age_max",type=int,help="max age to be processed")
parser.add_argument("-i","--inc",type=int,help='step inc for processing')

# default values for the date
year = 2017
#month = 7
#day = 1
date = 'Jun-01'
supertype = 'EAD'
quiet = False
vert = 'theta'
target = 'FullAMA'

# ages max in days
age_max = 62
step_inc = 6
#step_start=240
#step_end=240
step_start  = step_inc

area_pix = 0.1*0.1

args = parser.parse_args()
#if args.year is not None: year = args.year
#if args.month is not None: month = args.month
if args.type is not None: supertype = args.type
#if args.saf is not None:
#    if args.saf == 'O':
#        saf = ''
if args.vert is not None: vert = args.vert
if args.quiet is not None:
        if args.quiet=='y': quiet=True
        else: quiet=False
if args.full is not None:
        if args.full=='y':
            target = 'global'
if args.date is not None: date = args.date
if 'FULL' not in supertype: target = 'FullAMA'
if args.age_max is not None: age_max = args.age_max
if args.inc is not None: step_inc = args.inc
    
hmax = 24*(age_max + 10)
step_start = step_inc

print('target',target)

# other parameters 
#traj_dir = "/data/legras/flexout/STC/FORW"+saf    
#out_dir = "/data/legras/STC/STC-FORW"+saf+"-OUT"
traj_dir = "/data/legras/flexout/STC/FORWBox-meanhigh"
if 'ciclad' in socket.gethostname():
    out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"
elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
    out_dir = "/homedata/legras/STC/STC-FORWBox-meanhigh-OUT"

#%% Central section
run_type = supertype+'-Box-'+vert+'-'+target
    
#pile_sav_name = os.path.join(out_dir,'pile-save-'+run_type+'.pkl')
regional_sav_stream = join(out_dir,'regionalN-save-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')

# Manage the file that receives the print output
if quiet:
    # Output file
    print_file = join(out_dir,'out','regionalN-'+run_type+'-'+date+'.out')
    fsock = open(print_file,'w')
    sys.stdout=fsock
    
print("regional-sources> process "+run_type+' '+date)

run_dir = join(traj_dir,'FORW-'+supertype+'-'+str(year)+'-'+date)
sources = io107.readpart107(0,run_dir)
IDX_ORGN = sources['idx_orgn']
# barometric altitude of the source (not assumed to be above 55 hPa)
#id1 = sources['p']>22632.
#sources['alt'] = np.empty(sources['numpart'])
#sources['alt'][~id1] = z2(0.01*sources['p'][~id1])
#sources['alt'][id1] = z1(0.01*sources['p'][id1])
# potential temperature of the source
sources['thet'] = sources['t'] * (cst.p0/sources['p'])**cst.kappa
# initialized live record
sources['live'] = np.empty(sources['numpart'],dtype=bool)
sources['live'].fill(True)
# veryhigh filter generating a boolean slice among initial points
ct = sources['flag'] >> 24
sources['veryhigh'] = (ct == 9) | (ct == 13)
#sources['silviahigh'] = ((ct == 9) & (sources['x']>90.75)) | \
#        (((ct==8)|(ct==9)|(ct==13)) & (sources['x']<=90.75))
sources['silviahigh'] = (ct == 9) | (ct == 13) | (ct == 8)
print("regional-sources > date "+ date,np.sum(sources['veryhigh']),np.sum(sources['silviahigh']))

# weight
ww = area_pix * np.cos(np.deg2rad(sources['y']))

# kill the sources of 30 August 2017 at 11h (227h)
if date == 'Aug-21':
    sources['live'][sources['ir_start'] == 3600*227] = False
#  filtering high lat / high alt (spurious sources above 360K at lat>40N)
sources['live'][(sources['thet']>360) & (sources['y']>=40)] = False
# additional filtering in the European region 
sources['live'][(sources['thet']>360) & (sources['y']>=35) & (sources['x']<40)] = False

with gzip.open(join('..','mkSTCmask','MaskCartopy2-STCfine.pkl'),'rb') as f: 
    mm =pickle.load(f)

# Attribution of a region code to all the sources    
regcode = mm['regcode']
mask = mm['mask']
icx=icy=0.25
xlon0 = -10
ylat0 = 0
nlon = mm['nlons']
nlat = mm['nlats']
idy = np.clip(np.floor((sources['y'] - ylat0)/icy).astype(np.int),0,nlat-1)
idx = np.clip(np.floor((sources['x'] - xlon0)/icx).astype(np.int),0,nlon-1)
sources['regn'] = mask[idy,idx]

# target discretization (borrowed from transit.py)
minv = 322.5
maxv = 422.5
binv = 20
deltav = (maxv-minv)/binv
vcent = np.arange(minv+0.5*deltav,maxv,deltav)
vedge = np.arange(minv,maxv+0.01,deltav)

# source distribution in the same interval but with 200 bins
bins = 200 
deltas = (maxv-minv)/bins
scent = np.arange(minv+0.5*deltav,maxv,deltas)
sedge = np.arange(minv,maxv+0.01,deltas)
ids = np.floor((sources['thet']-sedge[0])/deltas).astype(np.int)
# eliminate possible sources outside the (minv,maxv) interval
sources['live'][ids<0] = False
sources['live'][ids>=bins] = False

# filtering of sources with mask = 0 
sources['live'][sources['regn']==0] = False

# Generation of the output dictionary
reg_src = {}
reg_src['supertype'] = supertype
reg_src['target'] = target
reg_src['date'] = date
reg_src['count'] = {}
reg_src['age'] = {}
for reg in regcode:
    reg_src['count'][reg] = np.zeros(shape=(binv,bins))
    reg_src['age'][reg] = np.zeros(shape=(binv,bins))

# Loop on steps
for step in range(step_start,hmax + step_inc ,step_inc):
    print("> step "+str(step))
    # Read the nactive parcels at current step
    data = io107.readpart107(step,run_dir,quiet=True)
    # Get the list of indexes for the active parcels after removing the offset
    idxsel = data['idx_back']-IDX_ORGN
    # Generate the list of ages of active parcels from their launch (in days)
    age = step/24 - sources['ir_start'][idxsel]/86400
    thet = data['t'] * (cst.p0/data['p'])**cst.kappa
    idv = np.floor((thet-vedge[0])/deltav).astype(np.int)
    # Killing exited parcels if FULL in supertype and if target=FullAMA
    # The killed parcels remain killed afterwards, even if they re-enter the domain
    if ('FULL' in supertype) & (target == 'FullAMA'):
        selkill = np.any([data['x']<=-10,data['x']>=160,data['y']<=0,data['y']>=50],axis=0)
        idxkill = idxsel[selkill]
        sources['live'][idxkill] = False
        del selkill
        del idxkill
    # Selecting parcels which are both of age less than age_max (in days) and have never left the FullAMA box
    # selage is a boolean array of dimension nactive
    selage = np.all([age <= age_max, sources['live'][idxsel] == True,idv>=0,idv<20,sources['silviahigh'][idxsel] == True],axis=0)
    idx2 = idxsel[selage]
    idv2 = idv[selage]
    age2 = age[selage]
    for i in range(len(idx2)):
        ii = idx2[i]
        #print(sources['regn'][ii],idv2[i],ids[ii])
        reg_src['count'][mm['regcode_inv'][sources['regn'][ii]]][idv2[i]][ids[ii]] += ww[ii]
        reg_src['age'][mm['regcode_inv'][sources['regn'][ii]]][idv2[i]][ids[ii]] += age2[i]*ww[ii]
    sys.stdout.flush()
    
    
sources.clear()
sys.stdout.flush()
pickle.dump(reg_src,gzip.open(regional_sav_stream,'wb',pickle.HIGHEST_PROTOCOL))
    
# %%
# Final processing and save 
#pile_name = os.path.join(out_dir,'pile-'+run_type+'.pkl')
#pile.complete()        
#pickle.dump(pile,gzip.open(pile_name,'wb',pickle.HIGHEST_PROTOCOL))
#transit.transit_show(pile)
if quiet: fsock.close()
