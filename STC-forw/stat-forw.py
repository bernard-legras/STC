# -*- coding: utf-8 -*-
"""
Created on 30 March 2018

This code process the forward runs to provide transit statistics

This version should work for all cases

Parameterized for Box-meanhigh on 15 Jan 2019

@author: Bernard Legras
"""

import os
import numpy as np
#import matplotlib.pyplot as plt
#from datetime import datetime, timedelta
import gzip, pickle
import sys
import resource
import argparse
from constants import kappa
#from subprocess import call
import transit as tt
import io107
from satratio import satratio
from zISA import z1,z2
p00 = 100000.
#%%
parser = argparse.ArgumentParser()
parser.add_argument("-y","--year",type=int,help="year")
parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
parser.add_argument("-t","--type",choices=["EAD","EAZ","EIZ","EID","EIZ-FULL","EID-FULL"],help="type")
#parser.add_argument("-s","--saf",choices=["O","N"],help="SAF version")
parser.add_argument("-v","--vert",choices=["theta","baro"],help="vertical discretization")
parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
parser.add_argument("-f","--full",type=str,choices=["y","n"],help="full (y) or not (n)")
# default values for the first start day
year = 2017
month = 7
day = 1
#saf = 'N'
supertype = 'EAZ'
water_path = True
quiet = False
vert = 'baro'
target = 'FullAMA'

# number of start dates to be processed
ndate = 6
dates = ["Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]

# for each start date, define the first time read, the last time read and 
# the increment, all in hour
step_start = 6
step_end  = 978
#step_start=240
#step_end=240
step_inc  = 6

args = parser.parse_args()
if args.year is not None:
    year = args.year
if args.month is not None:
    month = args.month
if args.type is not None:
    supertype = args.type
#if args.saf is not None:
#    if args.saf == 'O':
#        saf = ''
if args.vert is not None:
   vert = args.vert
if args.quiet is not None:
        if args.quiet=='y':
            quiet=True
        else:
            quiet=False
if args.full is not None:
        if args.full=='y':
            target = 'global'

if 'FULL' not in supertype:
    target = 'FullAMA'

print('target',target)

# other parameters 
#traj_dir = "/data/legras/flexout/STC/FORW"+saf    
#out_dir = "/data/legras/STC/STC-FORW"+saf+"-OUT"
traj_dir = "/data/legras/flexout/STC/FORWBox-meanhigh"    
out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"

# initialize transit dictionary
pile=tt.transit(water_path=water_path,vert=vert,target=target)
#### Restart from sav of the previous run if interrupted 
#pile=pickle.load(gzip.open('pile2-sav-R.pkl','rb'))
#%% Central section

#if saf == 'N':
run_type = supertype+'-Box-'+vert+'-'+target
#else:
#    run_type = supertype+'-O-'+vert+'-'+target
    
#pile_sav_name = os.path.join(out_dir,'pile-save-'+run_type+'.pkl')
pile_sav_stream = os.path.join(out_dir,'pile-save-stream-'+run_type+'.pkl')

# Manage the file that receives the print output
if quiet:
    # Output file
    print_file = os.path.join(out_dir,'out','pile-'+run_type+'.out')
    fsock = open(print_file,'w')
    sys.stdout=fsock
    
print("stat-forw> process "+run_type)

# Loop on ndate
for i in range(ndate):
    run_dir = os.path.join(traj_dir,'FORW-'+supertype+'-'+str(year)+'-'+dates[i])
    sources = io107.readpart107(0,run_dir)
    IDX_ORGN = sources['idx_orgn']
    if water_path :
        sources['rv_s'] = satratio(sources['p'],sources['t'])
    # barometric altitude of the source (not assumed to be above 55 hPa)
    id1 = sources['p']>22632.
    sources['alt'] = np.empty(sources['numpart'])
    sources['alt'][~id1] = z2(0.01*sources['p'][~id1])
    sources['alt'][id1] = z1(0.01*sources['p'][id1])
    # potential temperature of the source
    sources['thet'] = sources['t'] * (p00/sources['p'])**kappa
    # initialized live record
    sources['live'] = np.empty(sources['numpart'],dtype=bool)
    sources['live'].fill(True)
    # veryhigh filter generating a boolean slice among initial points
    ct = sources['flag'] >> 24
    sources['veryhigh'] = (ct == 9) | (ct == 13)
    sources['silviahigh'] = ((ct == 9) & (sources['x']>90.75)) | \
        (((ct==8)|(ct==9)|(ct==13)) & (sources['x']<=90.75))
    print("stat-forw> date "+ dates[i])
    
    # Loop on steps
    for step in range(step_start,step_end+step_inc ,step_inc):
        print("stat-forw> step "+str(step))
        data = io107.readpart107(step,run_dir,quiet=True)
        idxsel = data['idx_back']-IDX_ORGN      
        data['age'] = step/24 - sources['ir_start'][idxsel]/86400
        # killing exited parcels if FULL in supertype and if target=FullAMA
        if ('FULL' in supertype) & (target == 'FullAMA'):
            selkill = np.any([data['x']<=-10,data['x']>=160,data['y']<=0,data['y']>=50],axis=0)
            idxkill = idxsel[selkill]
            sources['live'][idxkill] = False
        # Selecting parcels which are both of age less than a month and have never left the FullAMA box
        selage = np.all([data['age'] <= 30, sources['live'][idxsel] == True],axis=0)
        # Reduction if the filter applies (only for the variables used in transit)
        if np.all(selage) != True:
            for var in ('x','y','p','t','age'):
                data[var] = data[var][selage]
            idxsel = idxsel[selage]
        data['x0'] = sources['x'][idxsel]
        data['y0'] = sources['y'][idxsel]
        data['p0'] = sources['p'][idxsel]
        data['t0'] = sources['t'][idxsel]
        data['thet0'] = sources['thet'][idxsel]
        data['alt0'] = sources['alt'][idxsel]
        # generate a veryhigh boolean slice among active parcels 
        data['veryhigh'] = sources['veryhigh'][idxsel]
        data['silviahigh'] = sources['silviahigh'][idxsel]
        if water_path :
            sources['rv_s'][idxsel] = np.minimum(sources['rv_s'][idxsel],satratio(data['p'],data['t']))
            data['rv_t'] = sources['rv_s'][idxsel]
        pile.update(data)
        data.clear()
        pile.transit['count'] +=1
        #pickle.dump(pile,gzip.open(pile_sav_name,'wb',pickle.HIGHEST_PROTOCOL))
        print('Memory used: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)+' (kb)')
    sources.clear()
    sys.stdout.flush()
    # save the pile to restart the run if needed
    pickle.dump(pile,gzip.open(pile_sav_stream,'wb',pickle.HIGHEST_PROTOCOL))
    
# %%
# Final processing and save 
pile_name = os.path.join(out_dir,'pile-'+run_type+'.pkl')
pile.complete()        
pickle.dump(pile,gzip.open(pile_name,'wb',pickle.HIGHEST_PROTOCOL))
#transit.transit_show(pile)
if quiet: fsock.close()