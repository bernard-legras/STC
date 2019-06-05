#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu 30 ay 2019

This script generates a part_000 file containing initial positions from subvisible cirrus clouds
as described in SEL2 files produced by Vincent Noel
    
@author: Bernard Legras
"""
import os
import socket
from datetime import datetime, timedelta
import numpy as np
import pickle, gzip
from itertools import chain
from scipy.interpolate import interp1d
from ECMWF_N import ECMWF, curtain
import io107
from netCDF4 import Dataset
import argparse

# %%    
if __name__ == '__main__': 
    """ Produce the sequence of initial positions and initial times to be used in
    the backward run. Generates output as part_000 file."""
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-d1","--day1",type=int,choices=1+np.arange(31),help="day1")
    parser.add_argument("-d2","--day2",type=int,choices=1+np.arange(31),help="day2")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
    
    # Default values
    base_year = 2017
    base_month = 7
    day1 = 1
    day2 = 1
    quiet = False
    verbose = False
    
    MISSING = -9999.0
    
    args = parser.parse_args()
    if args.year is not None: base_year = args.year
    if args.month is not None: base_month = args.month
    if args.day1 is not None: day1 = args.day1
    if args.day2 is not None: day2 = args.day2
    if args.quiet is not None:
        if args.quiet=='y':
            quiet=True
        else:
            quiet=False
            
    # Dates beginning and end
    date_beg = datetime(year=base_year, month=base_month, day=day1, hour=0)
    date_end = datetime(year=base_year, month=base_month, day=day2, hour=0)
    ref_date = datetime(1993,1,1,0)
    stamp_date = date_end
    # test that date_end > date_beg
    
    # Define main directories
    if 'ciclad' in socket.gethostname():
            SVC_Dir = '/bdd/CFMIP/SEL2'
            flexout = '/data/akottayil/flexout/STC/BACK'           
    elif socket.gethostname() == 'satie':
            SVC_Dir = '/home/legras/sandbox/SVC'
            flexout = '/home/legras/sandbox/SVC'
    part_dir = os.path.join(flexout,'BACK-SVC-EAD-'+date_beg.strftime('%b-%Y-day%d-')+date_end.strftime('%d-D01'))
    try:
        os.mkdir(part_dir)
    except:
        pass        
    
    # Generate the dictionary to be used to write part_000
    part0 = {}
    # Heading data
    part0['lhead'] = 3
    part0['outnfmt'] = 107
    part0['mode'] = 3   # modify that
    part0['stamp_date'] = date_end.year*10**10 + date_end.month*10**8 + \
        date_end.day*10**6 + date_end.hour*10**4 + date_end.minute*100
    part0['itime'] = 0
    part0['step'] = 450
    part0['idx_orgn'] = 1
    part0['nact_lastO'] = 0
    part0['nact_lastNM'] = 0
    part0['nact_lastNH'] = 0
    #part0['flag'] = np.empty(0,dtype=int)
    # we use lists as appending data is vey fast in this case
    for var in ['x','y','t','p','ptop','pbase','t','tE','tEbase','tEtop','ctop','cbase','idx_back','ir_start']:
        part0[var] = []

    numpart = 0
    cumulate_depth = 0
    
    # loop on the days
    date_p = date_beg
    while date_p <= date_end:
        try:
            os.chdir(os.path.join(SVC_Dir,date_p.strftime('%Y/%Y_%m_%d')))
        except:
            print('Directory is missing ',date_p)
            date_p += timedelta(days=1)
            continue
        ll = os.listdir()
        # loop on the orbits
        while len(ll)>0:
            file = ll.pop()
            print(file)
            ncid = Dataset(file)
            data = {}
            data['latitude'] = ncid.variables['latitude'][:]
            data['longitude'] = ncid.variables['longitude'][:]
            data['time'] = ncid.variables['Profile_Time'][:]
            # select part of the orbit uin the FullAMA domain
            selec = (data['latitude']>0) & (data['latitude']<50) & (data['longitude']<160) & (data['longitude']>0)
            # process only orbits that intersect the FullAMA domain
            if np.sum(selec)>0:
                data['latitude'] = data['latitude'][selec]
                data['longitude'] = data['longitude'][selec]
                data['time'] = data['time'][selec]
                
                for var in ['layer_top_altitude','layer_base_altitude','layer_temperature',\
                        'integrated_particulate_color_ratio','integrated_particulate_depolarization_ratio',\
                        'layer_optical_depth']:
                    data[var] = ncid.variables[var][selec,:]
                
                mean_date = ref_date + timedelta(seconds=data['time'].mean())
                # surrounding dates for ERA-Interim
                date1 = datetime(mean_date.year,mean_date.month,mean_date.day,mean_date.hour - mean_date.hour%3)
                date2 = date1 + timedelta(hours=3)
                dtt = (mean_date-date1).total_seconds()/10800

                # get ECMWF data for the orbit and generate a curtain over the retained segment
                dat1 = ECMWF('FULL-EI',date1)
                dat1._get_T()
                dat1._mkp()
                dat1._mkz()
                dat1.close()
                sect1 = dat1.interpol_orbit(data['longitude'],data['latitude'],varList=['P','T','Z'])
                del dat1
                dat2 = ECMWF('FULL-EI',date2)
                dat2._get_T()
                dat2._mkp()
                dat2._mkz()
                dat2.close()
                sect2 = dat2.interpol_orbit(data['longitude'],data['latitude'],varList=['P','T','Z'])
                del dat2
                sect = curtain()
                sect.x = sect1.x
                sect.y = sect1.y
                for var in ['P','T','Z']:
                    sect.var[var] = (1-dtt)*sect1.var[var] + dtt*sect2.var[var]
                del sect1
                del sect2
                # generate a list of interpolating functions (one for each profile)
                interpP = []
                interpT = []
                for i in range(len(data['longitude'])):
                    interpP.append(interp1d(sect.var['Z'][::-1,i],np.log(sect.var['P'][::-1,i]),copy=False,bounds_error=True,assume_sorted=True))
                    interpT.append(interp1d(sect.var['Z'][::-1,i],sect.var['T'][::-1,i],copy=False,bounds_error=True,assume_sorted=True))
                del sect
                
                # loop on the layers
                for layer in range(20):
                    # apply selection criteria for each layer
                    # diagnostics
                    if verbose:
                        selec_color = (data['integrated_particulate_color_ratio'][:,layer]>0.7) & \
                                 (data['integrated_particulate_color_ratio'][:,layer]<1.5)
                        selec_depolar = (data['integrated_particulate_depolarization_ratio'][:,layer]>0.1) & \
                                 (data['integrated_particulate_depolarization_ratio'][:,layer]<0.7)
                        selec_temperature = (data['layer_temperature'][:,layer]<-40) & (data['layer_temperature'][:,layer]>MISSING)
                        selec_optical_depth = (data['layer_optical_depth'][:,layer]<0.03) & (data['layer_temperature'][:,layer]>MISSING)
                        selec_color_MISSING = (data['integrated_particulate_color_ratio'][:,layer] == MISSING)
                        selec_depolar_MISSING = (data['integrated_particulate_depolarization_ratio'][:,layer] == MISSING)
                        selec_temperature_MISSING = (data['layer_temperature'][:,layer] == MISSING)
                        selec_optical_depth_MISSING = (data['layer_optical_depth'][:,layer] == MISSING)
                    # selection applying Martins et al criteria
                    selec2 = (data['integrated_particulate_color_ratio'][:,layer]>0.7) & \
                             (data['integrated_particulate_color_ratio'][:,layer]<1.5) & \
                             (data['integrated_particulate_depolarization_ratio'][:,layer]>0.1) & \
                             (data['integrated_particulate_depolarization_ratio'][:,layer]<0.7) & \
                             (data['layer_temperature'][:,layer]<-40) & \
                             (data['layer_optical_depth'][:,layer]<0.03) & \
                             (data['layer_base_altitude'][:,layer]>7.5)
                    # removing pixels with missing values 
                    selec2 = selec2 & (data['layer_temperature'][:,layer]>MISSING) & \
                             (data['layer_optical_depth'][:,layer]>MISSING)
                    # process the layer if the selection retains some features
                    if np.sum(selec2) >0:
                        # diagnostic prints
                        print('layer {:d} {:6.2f}%'.format(layer,100*np.sum(selec2)/len(selec2)))
                        if verbose:
                            print('{:6.2f}% {:6.2f}% {:6.2f}% {:6.2f}% M {:6.2f}% {:6.2f}% {:6.2f}% {:6.2f}%'.format(
                              100*np.sum(selec_color)/len(selec2), 100*np.sum(selec_depolar)/len(selec2),
                              100*np.sum(selec_temperature)/len(selec2), 100*np.sum(selec_optical_depth)/len(selec2),
                              100*np.sum(selec_color_MISSING)/len(selec2), 100*np.sum(selec_depolar_MISSING)/len(selec2),
                              100*np.sum(selec_temperature_MISSING)/len(selec2), 100*np.sum(selec_optical_depth_MISSING)/len(selec2)))
                        part0['x'].append(list(data['longitude'][selec2]))
                        part0['y'].append(list(data['latitude'][selec2]))
                        part0['t'].append(list(data['layer_temperature'][selec2,layer]))
                        part0['ctop'].append(list(data['layer_top_altitude'][selec2,layer]))
                        part0['cbase'].append(list(data['layer_base_altitude'][selec2,layer]))
                        cumulate_depth += int(np.sum(np.floor(1000*data['layer_top_altitude'][selec2,layer]-1000*data['layer_base_altitude'][selec2,layer])))
                        part0['ir_start'].append([(ref_date+timedelta(seconds=tt)-stamp_date).total_seconds() for tt in data['time'][selec2]])
                        numpart += np.sum(selec2)
                        # use the interpolation functions to get pressure and temperature from ECMWF data
                        lons = np.array(data['longitude'][selec2])
                        lats = np.array(data['latitude'][selec2])
                        # mean layer height (in m)
                        ctop = 1000*data['layer_top_altitude'][selec2,layer]
                        cbase = 1000*data['layer_base_altitude'][selec2,layer]
                        cmean = 0.5*(cbase+ctop)
                        # selection of interpolation fuctions
                        interpPsel = [interpP[k] for k in np.where(selec2)[0]]
                        interpTsel = [interpT[k] for k in np.where(selec2)[0]]
                        for i in range(len(lons)):
                            part0['p'].append(np.exp(interpPsel[i](cmean[i])))
                            part0['ptop'].append(np.exp(interpPsel[i](ctop[i])))
                            part0['pbase'].append(np.exp(interpPsel[i](cbase[i])))
                            part0['tE'].append(interpTsel[i](cmean[i]))
                            part0['tEtop'].append(interpTsel[i](ctop[i]))
                            part0['tEbase'].append(interpTsel[i](ctop[i]))
                    
        date_p += timedelta(days=1)
    
    # change the lists into np.arrays
    print ('final collection and write to part_000')
    print ('numpart ',numpart)
    print('cumulate_depth ',cumulate_depth)
    
    for var in ['x','y','t','ctop','cbase']:
        part0[var] = np.array(list(chain.from_iterable(part0[var]))).astype(np.float32)

    for var in ['p','ptop','pbase','tE','tEtop','tEbase']:
        part0[var] = np.array(part0[var])
    
    part0['ir_start'] = np.array(list(chain.from_iterable(part0['ir_start']))).astype(np.int32)
    part0['idx_back'] = np.arange(1,numpart+1)
    part0['flag'] = np.full(numpart,31,dtype=np.int32)
    print ('t vs tE',np.mean(part0['t']-part0['tE']+273.15),np.var(part0['t']-part0['tE']+273.15))
                  
    # final size information 
    part0['numpart'] = numpart
    part0['nact'] = numpart
    
    # Here we can decide to write the result or to go one step further by splitting each layer into a number of parcels 
    # In practice, we can launch one parcel per m over the depth of each layer
    # Generate the new dictionary to be used to write part_000
    partR = {}
    # Heading data
    partR['lhead'] = 3
    partR['outnfmt'] = 107
    partR['mode'] = 3   # modify that
    partR['stamp_date'] = date_end.year*10**10 + date_end.month*10**8 + \
        date_end.day*10**6 + date_end.hour*10**4 + date_end.minute*100
    partR['itime'] = 0
    partR['step'] = 450
    partR['idx_orgn'] = 1
    partR['nact_lastO'] = 0
    partR['nact_lastNM'] = 0
    partR['nact_lastNH'] = 0
    #partR['flag'] = np.empty(0,dtype=int)
    # we use lists as appending data is vey fast in this case
    partR['ir_start'] = np.empty(cumulate_depth,dtype=np.int32)
    partR['x'] = np.empty(cumulate_depth,dtype=np.float32)
    partR['y'] = np.empty(cumulate_depth,dtype=np.float32)
    partR['t'] = np.empty(cumulate_depth,dtype=np.float32)
    partR['p'] = np.empty(cumulate_depth,dtype=np.float32)
    
    # Loop over the number of parcels
    numpartR = 0
    for i in range(numpart):
        n = int(1000*(part0['ctop'][i]-part0['cbase'][i]))
        partR['x'][numpartR:numpartR+n] = part0['x'][i]
        partR['y'][numpartR:numpartR+n] = part0['y'][i]
        partR['ir_start'][numpartR:numpartR+n] = part0['ir_start'][i]
        partR['p'][numpartR:numpartR+n] = np.exp(np.linspace(np.log(part0['pbase'][i]),np.log(part0['ptop'][i]),n,dtype=np.float32))
        if verbose:
            if part0['tE'][i]<min(part0['tEtop'][i],part0['tEbase'][i]):
                print('ACHTUNG tropopause straddling ',part0['ctop'][i],part0['cbase'][i])
        partR['t'][numpartR:numpartR+n] = np.linspace(part0['tEbase'][i],part0['tEtop'][i],n,dtype=np.float32)
        numpartR += n
        
    partR['idx_back'] = np.arange(1,numpartR+1)
    partR['flag'] = np.full(numpartR,31,dtype=np.int32)
    
    # final size information 
    partR['numpart'] = numpartR
    partR['nact'] = numpartR
    print('numpartR ',numpartR)
    
    # write the result as part_000 file that can be used by TRACZILLA
    # Attention: ctop, cbase or tE are not written
    # make another output (by pickle dump) to svae such data
    if not os.path.exists(part_dir):
        os.makedirs(part_dir)
    newpart0 = os.path.join(part_dir,'part_000')
    io107.writeidx107(newpart0,partR)
    arkiv = os.path.join(part_dir,'arkiv.pkl')
    with gzip.open(arkiv,'wb') as f:
        pickle.dump([part0,partR],f)
