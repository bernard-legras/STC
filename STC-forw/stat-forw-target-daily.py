# -*- coding: utf-8 -*-
"""
Created on 20 May 2019

This code process the forward runs to provide transit statistics

This version is derived from stat-forw-uniq.py to provide statistics within
a given period in the target space, that is, instead of looping on all dates
of a stream, it loops on all streams for dates corresponding to a given
time interval.

As a matter of simplification there is one 10-day period associated with each
stream, as described in trueDates and the relevant streams for each period
are listed in preDates. The first period that can be investigated is Jul-11. 

@author: Bernard Legras
"""

import os
import numpy as np
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import gzip, pickle
import sys
import resource
import argparse
import socket
import transit as tt
import io107
from satratio import satratio
from zISA import z1,z2
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
parser.add_argument("-st","--stream",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')
parser.add_argument("-am","--age_max",type=int,help="max age to be processed")
parser.add_argument("-i","--inc",type=int,help='step inc for processing')

trueDates ={'Jun-01':[datetime(2017,6,1,0),datetime(2017,6,10)],
            'Jun-11':[datetime(2017,6,11,0),datetime(2017,6,20)],
            'Jun-21':[datetime(2017,6,21,0),datetime(2017,6,30)],
            'Jul-01':[datetime(2017,7,1,0),datetime(2017,7,10)],
            'Jul-11':[datetime(2017,7,11,0),datetime(2017,7,20)],
            'Jul-21':[datetime(2017,7,21,0),datetime(2017,7,31)],
            'Aug-01':[datetime(2017,8,1,0),datetime(2017,8,10)],
            'Aug-11':[datetime(2017,8,11,0),datetime(2017,8,20)],
            'Aug-21':[datetime(2017,8,21,0),datetime(2017,8,31)]}

preDates = {'Jun-01':None,'Jun-11':None,'Jun-21':None,'Jul-01':None,
            'Jul-11':['Jun-01','Jun-11','Jun-21','Jul-01','Jul-11'],
            'Jul-21':['Jun-11','Jun-21','Jul-01','Jul-11','Jul-21'],
            'Aug-01':['Jun-21','Jul-01','Jul-11','Jul-21','Aug-01'],
            'Aug-11':['Jul-01','Jul-11','Jul-21','Aug-01','Aug-11'],
            'Aug-21':['Jul-11','Jul-21','Aug-01','Aug-11','Aug-21']}

# default values for the date
year = 2017
#month = 7
#day = 1
stream = 'Jul-11'
supertype = 'EAD'
water_path = True
quiet = False
vert = 'theta'
target = 'FullAMA'

# ages max in days
age_max = 40
step_inc = 6
#step_start=240
#step_end=240
step_start  = step_inc

args = parser.parse_args()
if args.type is not None: supertype = args.type
if args.vert is not None: vert = args.vert
if args.quiet is not None:
        if args.quiet=='y': quiet=True
        else: quiet=False
if args.full is not None:
        if args.full=='y':
            target = 'global'
if args.stream is not None: stream = args.stream
if preDates[stream] == None:
    print('Not enough predates for this stream')
    exit()
if 'FULL' not in supertype: target = 'FullAMA'
if args.age_max is not None: age_max = args.age_max
if args.inc is not None: step_inc = args.inc
    
#hmax = 24*(age_max + 10)
#hmax = 1728

day1,day2 = trueDates[stream]

print('target',target)
print('stream and interval',stream,day1,day2)

# other parameters 
#traj_dir = "/data/legras/flexout/STC/FORW"+saf    
#out_dir = "/data/legras/STC/STC-FORW"+saf+"-OUT"
traj_dir = "/data/legras/flexout/STC/FORWBox-meanhigh"    
if 'ciclad' in socket.gethostname():
    out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"
elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
    out_dir = "/homedata/legras/STC/STC-FORWBox-meanhigh-OUT"
    
run_type = supertype+'-Box-'+vert+'-'+target

# initialize transit dictionary
pile = {}
dd = day1
while dd <= day2:
    pile[dd] = tt.transit(water_path=water_path,vert=vert,target=target)
    dd += timedelta(days=1)
#### Restart from sav of the previous run if interrupted 
#pile=pickle.load(gzip.open('pile2-sav-R.pkl','rb'))

#pile2_sav_stream = os.path.join(out_dir,'pile2-save-daily-target-'+run_type+'-'+stream+'-h'+str(hmax)+'.pkl')
pile_sav_stream = os.path.join(out_dir,'pile-save-daily-target-'+run_type+'-'+stream+'-age'+str(24*age_max)+'.pkl')

# Manage the file that receives the print output
if quiet:
    # Output file
    print_file = os.path.join(out_dir,'out','pile-daily-target-'+run_type+'-'+stream+'-age'+str(24*age_max)+'.out')
    fsock = open(print_file,'w')
    sys.stdout=fsock
    
print("stat-forw-target_daily> process "+run_type+' '+stream)


for date in preDates[stream]:    
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
    print("stat-forw-target-daily > date "+ date,np.sum(sources['veryhigh']),np.sum(sources['silviahigh']))

    # kill the sources of 30 August 2017 at 11h (227h)
    if date == 'Aug-21':
        sources['live'][sources['ir_start'] == 3600*227] = False
    #  filtering high lat / high alt (spurious sources above 360K at lat>40N)
    sources['live'][(sources['thet']>360) & (sources['y']>=40)] = False
    # additional filtering in the European region 
    sources['live'][(sources['thet']>360) & (sources['y']>=35) & (sources['x']<40)] = False

    #loop on the days
    dd = day1
    while dd <= day2:
        print('new day',dd)
        # Loop on steps
        dd1 = dd - trueDates[date][0]
        dd2 = dd1 + timedelta(days=1)
        #step1 = min(int(dd1.total_seconds()/3600),hmax)
        #step2 = min(int(dd2.total_seconds()/3600),hmax)
        step1 = int(dd1.total_seconds()/3600)
        step2 = int(dd2.total_seconds()/3600)
        print('step1, step2',step1,step2)
        for step in range(step1,step2,step_inc):
            if step ==0: 
                print("stat-forw_target-daily> date",date,"skip step 0")
                continue
            print("stat-forw_target-daily> date", date," step ",str(step))
            # Read the nactive parcels at current step
            data = io107.readpart107(step,run_dir,quiet=True)
            # Get the list of indexes for the active parcels after removing the offset
            idxsel = data['idx_back']-IDX_ORGN
            # Generate the list of ages of active parcels from their launch (in days)
            data['age'] = step/24 - sources['ir_start'][idxsel]/86400
            # Killing exited parcels if FULL in supertype and if target=FullAMA
            # The killed parcels remain killed afterwards, even if they re-enter the domain
            # ACHTUNG ACHTUNG
            # This does not work here since the part data are not read from the beginning of the
            # run. In order to avoid that, should read at the beginning of stream processing
            # a file that contains the maximum age of the parcel before it leaves the FullAMA 
            # domain. This information should be provided by another program.
            if ('FULL' in supertype) & (target == 'FullAMA'):
                selkill = np.any([data['x']<=-10,data['x']>=160,data['y']<=0,data['y']>=50],axis=0)
                idxkill = idxsel[selkill]
                sources['live'][idxkill] = False
                del selkill
                del idxkill
            # Selecting parcels which are both of age less than age_max (in days) and have never left the FullAMA box
            # and for which the source ia alive
            # selage is a boolean array of dimension nactive
            selage = np.all([data['age'] <= age_max, sources['live'][idxsel] == True],axis=0)
            if np.sum(selage) == 0:
                print('step beyond limiting age, no selected parcels, skipped')
                continue
            if np.sum(selage) < len(selage):
                print('selected parcels',np.sum(selage),'out of',len(selage))
                for var in ('x','y','p','t','age'):
                    data[var] = data[var][selage]
            else:
                print('all',len(selage),'parcels selected')
            idxsel = idxsel[selage]
            data['x0'] = sources['x'][idxsel]
            data['y0'] = sources['y'][idxsel]
            data['p0'] = sources['p'][idxsel]
            data['t0'] = sources['t'][idxsel]
            data['thet0'] = sources['thet'][idxsel]
            data['alt0'] = sources['alt'][idxsel]
            # generate a veryhigh boolean slice among live parcels 
            data['veryhigh'] = sources['veryhigh'][idxsel]
            data['silviahigh'] = sources['silviahigh'][idxsel]
            if water_path :
                sources['rv_s'][idxsel] = np.minimum(sources['rv_s'][idxsel],satratio(data['p'],data['t']))
                data['rv_t'] = sources['rv_s'][idxsel]   
            pile[dd].update(data)
            data.clear()
        dd += timedelta(days=1)
        sys.stdout.flush()
        print('Memory used: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)+' (kb)')
    sources.clear()    
    # save the pile to restart the run if needed
    #pickle.dump(pile,gzip.open(pile2_sav_stream,'wb',pickle.HIGHEST_PROTOCOL))

# Writing final result
print('writing final result')
pickle.dump(pile,gzip.open(pile_sav_stream,'wb',pickle.HIGHEST_PROTOCOL))
print('done')
# %%
# Final processing and save 
#pile_name = os.path.join(out_dir,'pile-'+run_type+'.pkl')
#pile.complete()        
#pickle.dump(pile,gzip.open(pile_name,'wb',pickle.HIGHEST_PROTOCOL))
#transit.transit_show(pile)
if quiet: fsock.close()
