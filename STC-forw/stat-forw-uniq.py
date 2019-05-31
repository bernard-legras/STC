# -*- coding: utf-8 -*-
"""
Created on 30 March 2018

This code process the forward runs to provide transit statistics

This version should work for all cases

Parameterized for Box-meanhigh on 15 Jan 2019
Forked from stat-forw to process a unique run
Needs a subsequent script stat-forw-gather to gather the data

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
#from subprocess import call
import transit as tt
import io107
from satratio import satratio
from zISA import z1,z2
import constants as cst

#%%
parser = argparse.ArgumentParser()
#parser.add_argument("-y","--year",type=int,help="year")
#parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
parser.add_argument("-t","--type",choices=["EAD","EAZ","EIZ","EID","EIZ-FULL","EID-FULL"],help="type")
#parser.add_argument("-s","--saf",choices=["O","N"],help="SAF version")
parser.add_argument("-v","--vert",choices=["theta","baro"],help="vertical discretization")
parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
parser.add_argument("-f","--full",type=str,choices=["y","n"],help="full (y) or not (n)")
parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')
parser.add_argument("-am","--age_max",type=int,help="max age to be processed")
parser.add_argument("-ami","--age_max_inter",type=int,help="max intermediate age to be processed")
parser.add_argument("-i","--inc",type=int,help='step inc for processing')

# default values for the date
year = 2017
#month = 7
#day = 1
date = 'Jun-01'
supertype = 'EAZ'
water_path = True
quiet = False
vert = 'baro'
target = 'FullAMA'

# ages max in days
age_max = 62
age_max_inter = 30
step_inc = 6
#step_start=240
#step_end=240
step_start  = step_inc

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
if args.age_max_inter is not None: age_max_inter = args.age_max_inter
if args.inc is not None: step_inc = args.inc
    
# derived times
#if date == 'Jul-21':
    #hmax = 24*(age_max + 11)
    #hinter = 24*(age_max_inter + 11)
#else:
hmax = 24*(age_max + 10)
hinter = 24*(age_max_inter + 10)
step_start = step_inc

print('target',target)

# other parameters 
#traj_dir = "/data/legras/flexout/STC/FORW"+saf    
#out_dir = "/data/legras/STC/STC-FORW"+saf+"-OUT"
traj_dir = "/data/legras/flexout/STC/FORWBox-meanhigh"    
out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"

# initialize transit dictionary
pile = tt.transit(water_path=water_path,vert=vert,target=target)
pile_inter = tt.transit(water_path=water_path,vert=vert,target=target)
#### Restart from sav of the previous run if interrupted 
#pile=pickle.load(gzip.open('pile2-sav-R.pkl','rb'))

#%% Central section
run_type = supertype+'-Box-'+vert+'-'+target
    
#pile_sav_name = os.path.join(out_dir,'pile-save-'+run_type+'.pkl')
pile_sav_stream = os.path.join(out_dir,'pile-save-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')
pile_sav_inter_stream = os.path.join(out_dir,'pile-save-stream-'+run_type+'-'+date+'-h'+str(hinter)+'.pkl')

# Manage the file that receives the print output
if quiet:
    # Output file
    print_file = os.path.join(out_dir,'out','pile-'+run_type+'-'+date+'.out')
    fsock = open(print_file,'w')
    sys.stdout=fsock
    
print("stat-forw> process "+run_type+' '+date)

run_dir = os.path.join(traj_dir,'FORW-'+supertype+'-'+str(year)+'-'+date)
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
print("stat-forw-uniq > date "+ date,np.sum(sources['veryhigh']),np.sum(sources['silviahigh']))

# kill the sources of 30 August 2017 at 11h (227h)
if date == 'Aug-21':
    sources['live'][sources['ir_start'] == 3600*227] = False

# Loop on steps
for step in range(step_start,hmax + step_inc ,step_inc):
    print("stat-forw> step "+str(step))
    # Read the nactive parcels at current step
    data = io107.readpart107(step,run_dir,quiet=True)
    # Get the list of indexes for the active parcels after removing the offset
    idxsel = data['idx_back']-IDX_ORGN
    # Generate the list of ages of active parcels from their launch (in days)
    data['age'] = step/24 - sources['ir_start'][idxsel]/86400
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
    selage = np.all([data['age'] <= age_max, sources['live'][idxsel] == True],axis=0)
    # A second filter is applied to intermediate age
    selage_inter = np.all([data['age'] <= age_max_inter, sources['live'][idxsel] == True],axis=0)
    # Copy of data and reduction if the filter applies (only for the variables used in transit)
    if step <= hinter:
        data_inter = {}
        for var in ('x','y','p','t','age'):
            data_inter[var] = data[var][selage_inter]
        idxsel_inter = idxsel[selage_inter]
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
    if step <= hinter:
        data_inter['x0'] = sources['x'][idxsel_inter]
        data_inter['y0'] = sources['y'][idxsel_inter]
        data_inter['p0'] = sources['p'][idxsel_inter]
        data_inter['t0'] = sources['t'][idxsel_inter]
        data_inter['thet0'] = sources['thet'][idxsel_inter]
        data_inter['alt0'] = sources['alt'][idxsel_inter]    
        data_inter['veryhigh'] = sources['veryhigh'][idxsel_inter]
        data_inter['silviahigh'] = sources['silviahigh'][idxsel_inter]
    if water_path :
        sources['rv_s'][idxsel] = np.minimum(sources['rv_s'][idxsel],satratio(data['p'],data['t']))
        data['rv_t'] = sources['rv_s'][idxsel]
        if step <= hinter: data_inter['rv_t'] = sources['rv_s'][idxsel_inter]
    pile.update(data)
    data.clear()
    pile.transit['count'] +=1
    if step <= hinter:
        pile_inter.update(data_inter)
        data_inter.clear()
        pile_inter.transit['count'] +=1
    #pickle.dump(pile,gzip.open(pile_sav_name,'wb',pickle.HIGHEST_PROTOCOL))
    print('Memory used: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)+' (kb)')
    if step == hinter:
        pickle.dump(pile_inter,gzip.open(pile_sav_inter_stream,'wb',pickle.HIGHEST_PROTOCOL))
sources.clear()
sys.stdout.flush()
pickle.dump(pile,gzip.open(pile_sav_stream,'wb',pickle.HIGHEST_PROTOCOL))
    
# %%
# Final processing and save 
#pile_name = os.path.join(out_dir,'pile-'+run_type+'.pkl')
#pile.complete()        
#pickle.dump(pile,gzip.open(pile_name,'wb',pickle.HIGHEST_PROTOCOL))
#transit.transit_show(pile)
if quiet: fsock.close()