#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun 25 March 2018

This script generates the part_000 file containing parcels initialized from the top of
convective clouds for a forward run. The parcels are generated from the SAFNWC data on a 0.1° grid
with a cutoff at 175 hPa.

The output is produced as a part_000 format 107 file which can be used in a StratoClim
backward run. The format 107 is used for convenience even if the flag field is of no use
here. flag is set to the uniform value 0x35 (=53) which means pressure coordinate (to be changed 
in a diabatic run), new parcel, SAFNWC grid and times relative to stamp-date. In addition mode=3     

Version N: uses the SFNWC-PTOP products calculated from the reprocessed archive in the SAFBox

@author: Bernard Legras
"""

from __future__ import division, print_function

import os
import socket
from sys import stdout
from datetime import datetime, timedelta
import numpy as np
import math
import pickle,gzip
#import pygrib
#from scipy.interpolate import RectBivariateSpline
import io107
from ECMWF_N import ECMWF
import argparse
pref = 101325.
# bbox, note it differs from transit.py, y is first 
target_range=np.array([[0.,50.],[-10.,160.]])
target_binx = 1700; target_biny = 500
deltay = (target_range[0,1]-target_range[0,0])/target_biny
deltax = (target_range[1,1]-target_range[1,0])/target_binx
ycent = np.arange(target_range[0,0] + 0.5*deltay,target_range[0,1],deltay)
xcent = np.arange(target_range[1,0] + 0.5*deltax,target_range[1,1],deltax)

x0=-10.
y0=0.
dxy=0.25
pixel_per_degree = 4
max_ix = 680
max_jy = 200

# definition zone for default parameters

R = 287.04
Cp = 1005.7
kappa = R/Cp
p0 = 100000.

# %%    
if __name__ == '__main__': 
    """ Produce the sequence of initial positions and initial times to be used in
    the backward run. Generates output as part_000 file."""
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-fd","--first_day",type=int,choices=1+np.arange(31),help="first day")
    parser.add_argument("-ld","--last_day",type=int,choices=1+np.arange(31),help="last day")      
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
    
    # Default values
    base_year = 2017
    base_month = 8
    base_day1 = 1
    base_day2 = 1
    cut_level = 175.1
    quiet = False
    
    #saf_levels = [175.,150.,125.,100.,75.]
    #p_interp = [17500.,15000.,12500.,1000.,7500.]
    
    args = parser.parse_args()
    if args.year is not None:
        base_year = args.year
    if args.month is not None:
        base_month = args.month
    if args.first_day is not None:
        base_day1 = args.first_day
    if args.last_day is not None:
        base_day2 = args.last_day
    if args.quiet is not None:
        if args.quiet=='y':
            quiet=True
        else:
            quiet=False
    
    # Needs to start at 0h next month to be on an EN date, avoids complicated adjustment
    # for negligible effect. 
    date_beg = datetime(year=base_year, month=base_month, day=base_day1, hour=0)
    date_end = datetime(year=base_year, month=base_month, day=base_day2, hour=23)
    # Test on gort
    #if 'gort' in socket.gethostname():
    #    date_beg = datetime(year=2005,month=1,day=1,hour=15)
    #    date_end = datetime(year=2005,month=1,day=1,hour=21)
    
    # Define main directories
    if 'ciclad' in socket.gethostname():
            SAF_dir = '/homedata/legras/STC/STC-SAFNWC-OUT'
            flexout = '/data/legras/flexout/STC/FORWN'
            #part_dir = os.path.join(flexout,'BACK-EAZ-'+date_beg.strftime('%b')\
            #                        +'-'+str(int(base_level))+'hPa')
            part_dir = os.path.join(flexout,'FORW-EAZ-'+date_beg.strftime('%Y-%b-%d'))
            try:
                os.mkdir(part_dir)
            except:
                pass
    else:
        'THIS CODE CANNOT RUN ON THIS COMPUTER'
        exit()
        
    # Generate the dictionary to be used to write part_000
    part0 = {}
    # Heading data
    part0['lhead'] = 3
    part0['outnfmt'] = 107
    part0['mode'] = 3   # modify that
    part0['stamp_date'] = date_beg.year*10**10 + date_beg.month*10**8 + \
        date_beg.day*10**6 + date_beg.hour*10**4 + date_beg.minute*100
    part0['itime'] = 0
    part0['step'] = 450
    part0['idx_orgn'] = 1
    part0['nact_lastO'] = 0
    part0['nact_lastNM'] = 0
    part0['nact_lastNH'] = 0
    part0['ir_start'] = np.empty(0,dtype=int)
    part0['x'] = np.empty(0,dtype=float)
    part0['y'] = np.empty(0,dtype=float)
    part0['t'] = np.empty(0,dtype=float)
    part0['p'] = np.empty(0,dtype=float)
    part0['flag'] = np.empty(0,dtype=int)
    
    # Generate the grid of points that will be used with the masked array
    # trick to generate the coordinates of points to be launched
    xg = np.tile(xcent,(target_biny,1))
    yg = np.tile(ycent,(target_binx,1)).T
    bloc_size = target_binx * target_biny
    
    # first date
    current_date = date_beg
    
    # No need to preread ECMWF files as we read also hourly SAFNWC data 
    # and no temporal interpolation is required
        
    # Initialize numpart
    numpart = 0
    
    # loop on the time
    # made easy because we use hourly data for both sources
    while current_date <= date_end:
        #print('date',current_date)
        FAFfile = current_date.strftime('SAFNWC-PTOP-%Y-%m-%d-%H:%M.pkl')
        fullname = os.path.join(SAF_dir,current_date.strftime('%Y/%m/%Y-%m-%d'),FAFfile)
        # First try to get the SAFNWC file
        try:
            with gzip.open(fullname,'rb') as file:
                [ptop,cloud_flag] = pickle.load(file)
        except IOError: 
            print('date',current_date,' file not found ',FAFfile)
            print(current_date,'is skipped')
            current_date += timedelta(hours=1)
            continue
        # The get the ECMWF data
        data = ECMWF('STC',current_date)
        data._get_T()
        data._mkp()
        data.close()
        # find time
        ir_start = int((current_date - date_beg).total_seconds())
        # process the cloud tops above 175 hPa
        ptop_temp = np.ma.masked_greater(ptop,cut_level)
        # apply same mask to other fields
        x_temp = np.ma.array(data=xg,mask=ptop_temp.mask)
        y_temp = np.ma.array(data=yg,mask=ptop_temp.mask)
        cloud_flag_temp = np.ma.array(cloud_flag,mask=ptop_temp.mask)
        num_temp = ptop_temp.count()
        # compressed 1d data
        # pressure converted to Pa)
        p_1d = ptop_temp.compressed()*100
        x_1d = x_temp.compressed()
        y_1d = y_temp.compressed()
        cloud_flag_1d = cloud_flag_temp.compressed()
        # calculate the index on the grid, clipping the round-off errors
        ix = np.clip(np.floor((x_1d-x0)*pixel_per_degree),0,max_ix).astype(np.int)
        jy = np.clip(np.floor((y_1d-y0)*pixel_per_degree),0,max_jy).astype(np.int)
        tt_1d = np.empty(num_temp)
        # vertical interpolation of temperature on a close column of the ECMWF grid
        # (lower left corner of the embedding mesh cell)
        for i in range(num_temp):
            tt_1d[i] = np.interp(math.log(p_1d[i]),np.log(data.var['P'][:,jy[i],ix[i]]),data.var['T'][:,jy[i],ix[i]])           
        numpart += num_temp
        part0['x'] = np.append(part0['x'],x_1d) 
        part0['y'] = np.append(part0['y'],y_1d)
        part0['t'] = np.append(part0['t'],tt_1d)
        part0['p'] = np.append(part0['p'],p_1d)
        part0['ir_start'] = np.append(part0['ir_start'],np.full(num_temp,ir_start,dtype=int))
        # add the flag for good and opaq clouds to the standard value 53=0x35 
        # SAFNWC + new parcel + time relative to stampdate 
        part0['flag'] = np.append(part0['flag'],cloud_flag_1d*0x8000000+53)
        # increment of the date
        print('date',current_date,' numpart',numpart)
        current_date = current_date + timedelta(hours=1)
        stdout.flush()
    
    # final size information 
    part0['idx_back'] = (1+np.arange(numpart)).astype(np.int)
    part0['numpart'] = numpart
    part0['nact'] = numpart
    
    # write the result as part_000 file 
    if not os.path.exists(part_dir):
        os.makedirs(part_dir)
    newpart0=os.path.join(part_dir,'part_000')
    io107.writeidx107(newpart0,part0)
