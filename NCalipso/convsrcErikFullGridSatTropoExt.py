#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the air sampled by CALIOP/CLOUDSAT
using the DARDAR sampling geometry to produce a curtain of parcels.

This version is based on the GridSat product and compares the temperature of the
parcel to the IR brightness temperature.
It calculates the tropopause location every 3 hours and does not consider cloud
intersections for parcels which are above the cold point

The selected data from DARDAR are organized as a list of parcels with initial
location (x,y,p) and temperature and time of observation (with respect to a
reference time). The time ordering is from the end of the month to the beginning
but it is not strictly monotonic because it is increasing for each orbit, so that
is has a sawtooth profile. However this ordering does not matter here.

This program determines the hits with convective clouds as described by GridSat
brightness temperature. The criterion is that the parcel must be in the troposphere,
that is below the WMO tropopause from ERA5 and at a temperature warmer than the
top of the cloud - a threshold (5K corresponding to 1 km (Minnis and Sherwood
criterion)).

As the traczilla program only stores active parcels in the outputs for compactness
the parcels are indexed with the ir_start field which makes possible to follow them
in time. A shiff by IDX_ORGN which can be 1 or 0 is necessary to use python addressing
since the index starts at 1. Sorry, this was initially thought to be exploited
with matlab.

In addition the program makes moisture calculations to check maximum and minimum
qsat along the path. It also tests whether the parcel becomes clear during the path
if it is initially cloudy with several thresholds (this is to test convective
versus in situ origin of cirrus). The DARDAR IWC is used here and also the SVC
data from V. Noel & H. Chepfer (subvisible cirrus).

Various other tasks are performed to handle deadborne parcels, parcels exiting
by the bottom or the top and old parcels.

Created on Wednesday 6 June 2019
Adapted on Tuesday 14 September 2021
Adapted again on 16 October 2025 from a version modified by Erik Johansson
   - flammkuchen replaced by h5py and two dedicated functions to write and
   read prod0, the archive contains meta data providing a description of its
   content
   - a detailed diagnostic of moisture evolution along the trajectory, introducing
   new flags and diagnostics
   - ages now all in days (correcting a previous error that made them
                           inconsistant for old parcels)
   - improvement and reorganisation of the diagnostics for clarification
   - better organisation or the restart

@author: Bernard Legras
"""

import socket
import numpy as np
from collections import defaultdict
from numba import jit
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import os
import sys
import argparse
import psutil
#import flammkuchen as fk  # obsolete and deprecated
import h5py
import gzip, pickle
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator  #
import time
import glob
#sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
from ECMWF_N import ECMWF
import geosat  #
from io107 import readpart107, readidx107

""" Global fixed parameters """
from createFlexpartFiles import max_life_time
p0 = 100000.
I_DEAD = 0x200000 #: Dead
I_HIT = 0x400000 #: Hit a cloud
I_CROSSED = 0x2000000 #: outside domain
I_DBORNE =  0x1000000 #: Lost after first step
I_OLD = 0x800000 #: Reached end of time without encounter a cloud
I_STOP = I_HIT + I_DEAD
# New flags for this version
I_CLOUD = 0x4000000
I_SVC = 0x8000000
''' Flags for clear encounter. The highest flag sets the other. '''
I_CLEAR = {'C06':0x10000000, 'C08':0x30000000,\
         'C10':0x70000000, 'C12':0xF0000000, 'C14':0x1F0000000}

# low p cut in the traczilla runs
lowpcut = 50
# highpcut in the traczilla runs
highpcut = 50000

# if True print a lot oj junk
verbose = False
debug = False

# idx_orgn was set to 1 in the prep file
IDX_ORGN = 1

# Error handling
class BlacklistError(Exception):
    pass

# Enter here the dates of the blacklisted GridSat files
blacklist = []

#%%
"""@@@@@@@@@@@@@@@@@@@@@@@@   FUNCTIONS   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

#%%
def satratio(p,T):
    """ Calculate the mass saturation ratio from pressure (in Pa) and temperature
    (in K).
    Taken from the book of K. Emanuel with corrections (origin?)
    estar is in hPa
    This code is less accurate at low temperatures near tropical tropopause than
    the formula used in esati_murphy that should be preferred."""
    estar = 1.0008*np.exp(23.33086-(6111.72784/T)+0.15215*np.log(T))
    satr = 0.622 * 1.0021 * estar/(0.01*p-estar)
    return satr
vsatratio = np.vectorize(satratio)

def esati_murphy(p, T):
    """ei in Pa saturation vapor pressure with respect to hexagonal (most stable)
    ice
    source: Murphy and Koop 2005, QJRMS """
    lnP = 9.550426-5723.265/T+3.53068*np.log(T)-0.00728332*T
    estar = np.exp(lnP)
    esati_murphy = estar / (p-estar)
    #: Convert to kg/kg
    esati_murphy = esati_murphy * 0.622
    return esati_murphy
vesati_murphy = np.vectorize(esati_murphy)

#%%
""" Functions related to the parcel data """

def get_slice_part(part_a,part_p,live_a,live_p,current_date,dstep,slice_width):
    """ Generator to generate 5' slices along flight track
        If nb_slices is 36, the position of the particles on a segemnt between
        a being 0 and p being 36 are at 0.5, 1.5, 2.5, ... 35.5
    """
    nb_slices = int(dstep/slice_width)
    ta = current_date + dstep
    tp = current_date
    tf = ta
    # test if no data available
    empty_live = (live_a.sum() == 0)
    for i in range(nb_slices):  # @UnusedVariable
        ti = tf- slice_width
        # note that 0.5*(ti+tf) cannot be calculated as we are adding two dates
        # tmid is the middle of a segment of size slice_width
        tmid = ti+0.5*(tf-ti)
        # interpolation at tmid from ante and post
        coefa = (tmid-tp)/dstep
        coefp = (ta-tmid)/dstep
        dat = {}
        dat['time'] = tmid
        if empty_live:
            dat['idx_back'] = dat['x'] = dat['y'] = dat['p'] = dat['t'] = []
            dat['itime']= None
        else:
            dat['idx_back'] = part_a['idx_back'][live_a]
            dat['x'] = coefa*part_a['x'][live_a] + coefp*part_p['x'][live_p]
            dat['y'] = coefa*part_a['y'][live_a] + coefp*part_p['y'][live_p]
            dat['p'] = coefa*part_a['p'][live_a] + coefp*part_p['p'][live_p]
            dat['t'] = coefa*part_a['t'][live_a] + coefp*part_p['t'][live_p]
            dat['itime'] = int(coefa*part_a['itime'] + coefp*part_p['itime'])

        tf -= slice_width
        yield dat

#%%
""" Function managing the exiting parcels
    In this version, the box range is not used.
    No lateral exit test is performed.
    Notice that even if the origin of the
    exit is not determinde, the parcel is
    flagged as crossed."""

@jit(nopython=True)
def exiter(itime, x,y,p,t,idx_back, flag,xc,yc,pc,tc,age,\
           ir_start, rr):
    nexits = 0
    for i in range(len(x)):
        i0 = idx_back[i]-IDX_ORGN
        if flag[i0] & I_DEAD == 0:
            nexits += 1
            xc[i0] = x[i]
            yc[i0] = y[i]
            tc[i0] = t[i]
            pc[i0] = p[i]
            age[i0] = (ir_start[i0] - itime)/86400
            #if   y[i] < rr[1,0] + 4.: excode = 6
            #elif x[i] < rr[0,0] + 4.: excode = 3
            #elif y[i] > rr[1,1] - 4.: excode = 4
            #elif x[i] > rr[0,1] - 4.: excode = 5
            if p[i] > highpcut - 150: excode = 1
            elif p[i] < lowpcut  + 15 : excode = 2
            else:                   excode = 7
            flag[i0] |= (excode << 13) + I_DEAD + I_CROSSED
    return nexits

#%%
""" Function providing the spherical distance
x1, y1: lon & lat first series of points
x2, y2: lon & lat second series of points
Using haversine formula and returning 2*sqrt(haversine)
We do not do the arc sine and keep the sinus of theta/2 
multiplied per 2 as an estimate of theta.
Warning: the initial version before the output to prod0
was performed, contained 0.5*sqrt(haversine)
"""

fff = np.pi/180
@jit(nopython=True)
def spher_dist(x1,y1,x2,y2):
    haversine = lambda theta: (1-np.cos(theta))/2
    y1r = fff * y1
    y2r = fff * y2
    return 2*np.sqrt(haversine(y2r-y1r)\
            + np.cos(y1r)*np.cos(y2r)*haversine(fff*(x2-x1)))

#%%
""" Function doing the comparison between parcels and clouds and setting the result field
In this version, we check only whether the temperature is warmer than the brightness temperature
from the infrared window. """

@jit(nopython=True)
def convbirth(itime, x,y,p,t,idx_back, flag,ptrop, xc,yc,pc,tc,age, BT, ir_start, x0,y0,stepx,stepy,binx,biny):
    #all_arguments = locals()
    #for (key,val) in all_arguments.items():
    #    print(key,type(val))
    nhits = 0
    for i in range(len(x)):
        # do not consider parcels which are above the wmo tropopause
        if ptrop[i] > p[i]: continue
        # find corresponding
        idx = min(int(np.floor((x[i]-x0)/stepx)),binx-1)
        idy = int(np.floor((y[i]-y0)/stepy))
        # does not process parcels outside of the GridSat domain
        # should not happen but as a precaution
        if (idy<0) | (idy>biny-1): continue
        # Test that the (shifted) brightness temperature is colder
        # than parcel temperature
        if BT[idy,idx] < t[i]:
            # find index of the parcel in the source list
            i0 = idx_back[i]-IDX_ORGN
            # if not already dead
            # Notice this filter could be applied ahead
            if flag[i0] & I_DEAD == 0:
                nhits += 1
                flag[i0] |= I_HIT + I_DEAD
                xc[i0] = x[i]
                yc[i0] = y[i]
                tc[i0] = t[i]
                pc[i0] = p[i]
                age[i0] = (ir_start[i0] - itime)/86400
    return nhits

#%%
""" Function related to satellite and ERA5 read """


def read_sat(t0,dtRange,pre=False,vshift=0):
    """ Generator reading the satellite data.
    The loop is infinite; sat data are called when required until the end of
    the parcel loop. """
    # initial time for internal loop
    current_time = t0
    while True:
        try:
            if (current_time in blacklist): raise BlacklistError()
            # defining a new object is necessary to avoid messing dat if an error occurs
            dat0 = geosat.GridSat(current_time)
            # Handle the return from geosat if the file is not found
            if dat0 is None: raise FileNotFoundError('WARNING')
            dat0._get_IR0()
            dat0.close()
            print('read GridSat for ',current_time)
            # The data are read as a masked array which is not convenient
            # for convbirth
            # extract the mask and the data
            mask = dat0.var['IR0'].mask
            dat0.var['IR0'] = dat0.var['IR0'].data
            if np.sum(mask) >0:
                print(np.sum(mask),' pixels masked')
                dat0.var['IR0'][mask] = 9999
            selneg = dat0.var['IR0'] < 0
            if np.sum(selneg) >0:
                print(np.sum(selneg),' negative pixels')
                dat0.var['IR0'][selneg] = 9999

            # remove dat and make it a view of dat0, try to avoid errors on first attempt

            try:
                del dat
            except:
                pass # intended for the first pass when dat undefined
            # Make shift below if needed, presently shift not used
            # Start shift at 10
            if vshift == 10:
                # make temperatures colder by 5K, approximately +1km
                dat0.var['IR0'] -= 5
            dat = dat0
            if pre:
                dat.tf = current_time + dtRange
                dat.ti = current_time
            else:
                dat.tf = current_time
                dat.ti = current_time - dtRange
        except BlacklistError:
            print('blacklisted date for GridSat',current_time)
            print('extend the lease')
            # extend the lease while keeping the old dat
            dat.ti -= dtRange
        except FileNotFoundError:
            print('GridSat file not found ',current_time)
            print('extend the lease')
            dat.ti -= dtRange
        # To collect other exceptions of if the planed management of exceptions
        # does not work as intended.
        except:
            print('Other exception hindering GridSat read ',current_time)
            print('extend the lease')
            dat.ti -= dtRange
        current_time -= dtRange

        yield dat

def read_ERA5(t0,dtRange,pre=False,vshift=0):
    """ Generator reading the ERA5 data.
    The loop is infinite; ERA5 data are called when required until the end of
    the parcel loop.
    The output contain a function that gives the pressure of the tropo^pause as
    a function of lon and lat in the FullAMA domain"""
    # initial time for internal loop
    current_time = t0
    while True:
        try:
            if (current_time in blacklist): raise BlacklistError()
            # defining a new object is necessary to avoid messing dat if an error occurs
            dat = ECMWF('FULL-EA',current_time)
            print('read ERA5 current-time', current_time)
            dat._get_T()
            dat._mkp()
            dat._mkz()
            # shift the origin to -180
            dats = dat.shift2west(-180)
            # extraction encompasses the GridSat domain
            datr0 = dats.extract(latRange=[-70,70],varss=['P','T','Z'])
            del dat, dats
            print('get ERA5 tropopause for ',current_time)
            usePcold = False
            if usePcold:
                print("Use coldpoint tropopause")
                datr0._CPT()
                # add a right column to 'pcold'
                ptropp = np.append(datr0.d2d['pcold'].T,[datr0.d2d['pcold'][:,0]],0).T
            else:
                print("Use WMO tropopause")
                datr0._WMO()
                # add an extra meridian to fill the [-180 180] domain
                ptropp = np.append(datr0.d2d['pwmo'].T,[datr0.d2d['pwmo'][:,0]],0).T
            datr0.fP = RegularGridInterpolator((np.arange(-70,71), np.arange(-180,181)), ptropp, bounds_error=True)
            try:
                del datr
            except:
                pass # intended for the first pass when datr undefined
            datr = datr0 # datr as a view of datr0
            # We reproduce here the same sequence as for the satellite file
            # although perhaps not appropriate (the ERA5 data are instantaneous)
            if pre:
                datr.tf = current_time + dtRange
                datr.ti = current_time
            else:
                datr.tf = current_time
                datr.ti = current_time - dtRange
        except BlacklistError:
            print('blacklisted date for GridSat',current_time)
            # extend the lease while keeping the old dat
            datr.ti -= dtRange
        except FileNotFoundError:
            print('ERA5 file not found ',current_time)
            # extend the lease while keeping the old dat
            datr.ti -= dtRange
        current_time -= dtRange

        yield datr

def check(dat,t):
        # check that the zone is not expired
    try:
        test = (t > dat.ti) and (t <= dat.tf)
    # Exception for the first usage when the time keys are not defined
    except KeyError:
        print('check KeyError')
        test = False
    except AttributeError:
        print('check AttributeError, normal at first step only')
        test = False
    return test

#%%
def findLastFile(dn):
    """ Finds the last part file available in the output directory of
    traczilla. This is meant to be compared to hmax. """
    fs = glob.glob(os.path.join(dn,"part_*"))
    fs.sort()
    i = -1
    for f in fs:
        i = i + 1
        num = int(os.path.basename(f).replace('.gz', '').split('_')[-1])
        if i == 0:
            stp0 = num
        elif i == 1:
            stp = num - stp0
        if i == 0:
            maxnum = 0
        else:
            maxnum = np.max((maxnum, num))
    return maxnum, stp
#%%
def h5prod0_write(file):
    """ Output of prod0 to an HDF5 file mimicking the structure of the
    dictionary with a description in the attributes of the file and groups."""
    h5f = h5py.File(file, 'w')
    h5f.attrs['creator'] = 'convsrcErikFullGridSatTropoExt.py'
    # This is complicated because HDF5 has no date type (or it is broken)
    arr = np.array([np.datetime64(datetime.now())])
    h5f.attrs['date'] = arr.astype(h5py.opaque_dtype(arr.dtype))
    h5f.attrs['run'] = outnames
    diag = h5f.create_group('diag')
    diag.attrs['object'] = 'Diagnostic group'
    diag.attrs.update(prod0['attr'])
    nhitsl = diag.create_dataset('nhitsl',data=prod0['diag']['nhitsl'])
    nhitsl.attrs['object'] = 'List of hits during the analysis of the run'
    nexitsl = diag.create_dataset('nexitsl',data=prod0['diag']['nexitsl'])
    nexitsl.attrs['object'] = 'List of exits during the analysis of the run'
    noldl = diag.create_dataset('noldl',data=prod0['diag']['noldl'])
    noldl.attrs['object'] = 'List of old during the analysis of the run'
    nnewl = diag.create_dataset('nnewl',data=prod0['diag']['nnewl'])
    nnewl.attrs['object'] = 'List of new parcel during the analysis of the run'
    dmean = diag.create_dataset('dist_mean',data=prod0['diag']['dist_mean'])
    dmax  = diag.create_dataset('dist_max', data=prod0['diag']['dist_max'])
    dstd  = diag.create_dataset('dist_std', data=prod0['diag']['dist_std'])
    dmean.attrs['object'] = 'Mean distance between ante/post pairs'
    dmax.attrs['object']  = 'Max distance between ante/post pairs'
    dstd.attrs['object']  = 'Std distance between ante/post pairs'
    dmean.attrs['units'] = '2*sqrt(haversine)'
    flags = diag.create_dataset('flags', data=prod0['flag_source'],\
                        compression="gzip", compression_opts=6)
    flags.attrs['object'] = 'Flags for run history and analysis. Check doc'
    for group in ('src','lowest_sat','highest_sat'):
        grp = h5f.create_group(group)
        for var in prod0[group]:
            grp.create_dataset(var, data=prod0[group][var],\
                               compression="gzip", compression_opts=6)
    h5f['src'].attrs['object'] = 'Group describing the hits (location and rvs value)'
    h5f['src'].attrs['content'] = 'lon (x), lat (y), pressure (p), temperature (t) '\
                                  +'age since release (age), saturation ratio (rvs)'
    h5f['src'].attrs['units'] = 'lat and lon: degrees, pressure: Pa, temperature: K '\
                                +'age: day, rvs: kg/kg'
    h5f['src'].attrs['vshift'] = vshift
    h5f['lowest_sat'].attrs['object'] = 'Group describing the lowest saturation met '\
                                    +'by the trajectory (location, value and age)'
    h5f['lowest_sat'].attrs['content'] = 'Same content as src group'
    h5f['highest_sat'].attrs['object'] = 'Group describing the highest saturation met '\
                                    +'by the trajectory (location, value and age)'
    h5f['highest_sat'].attrs['content'] = 'Same content as src group'
    fc = h5f.create_group('first_clear')
    for off in qoff.keys():
        grp = fc.create_group(off)
        for var in prod0['first_clear'][off]:
            grp.create_dataset(var, data=prod0['first_clear'][off][var],
                               compression="gzip", compression_opts=6)
    fc.attrs['Object'] = 'Group describing the first encounter of clear sky for '\
                         +'particles emitted from the cloudy pixels acoording to '\
                         +'DARDAR or the SVC database'
    fc.attrs['Content'] = 'There are 5 subgroups according to the threshold, each '\
                          +'containing the same data as the src group to characterize '\
                          +'the first encounter with the given threshold.'
    fc.attrs['Method'] = 'The test is qsat > threshold * q0 where q0 is the initial total '\
                            +'moisture, that is iwc + qsat0 for DARDAR cloudy pixels '\
                            +'and qsat0 for SVC pure pixels (not cloudy in DARDAR).'
    fc.attrs['Thresholds'] = '0.6, 0.8, 1.0, 1.2, 1.4'
    h5f.close()

def h5prod0_read(file):
    """ Revert  h5prod0_write by reconstructing prod0 from the HDF5 file """

    h5f = h5py.File(file, 'r')
    prod0 = defaultdict(dict)
    prod0['attr'] = {}
    for var in h5f['diag'].attrs.keys():
        prod0['attr'][var] = h5f['diag'].attrs[var]
    prod0['diag'] = {}
    for var in h5f['diag'].keys():
        if var != 'flags':
            prod0['diag'][var] = list(h5f['diag'][var][:])
    prod0['flag_source'] = h5f['diag']['flags'][:]
    for group in ('src','lowest_sat','highest_sat'):
        prod0[group] = {}
        for var in h5f[group].keys():
            prod0[group][var] = h5f[group][var][:]
    prod0['first_clear'] = {}
    for off in qoff.keys():
        prod0['first_clear'][off] = {}
        for var in h5f['first_clear'][off].keys():
            prod0['first_clear'][off][var] = h5f['first_clear'][off][var][:]
    h5f.close()
    return prod0

#%%
"""@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   MAIN   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

if __name__ == '__main__':
    # global IDX_ORGN
    parser = argparse.ArgumentParser()
    # No longer optional
    #parser.add_argument("-d","--use_dardar", action='store_true', default = True,
    #                    help = "Use DARDAR data")
    parser.add_argument("-y","--year", type = int, default = 2008,\
                        help = "year. Default = 2008")
    parser.add_argument("-m", "--month", type=int, choices=np.arange(1, 13),\
                        default=7, help = "Month. Default = 7")
    parser.add_argument("-n", "--night", type=str, choices=["d","n", "a"],\
        default='n',help = "pixlar with day (d), night (n) or all (a). Default = n")
    #parser.add_argument("-c","--clean0",type=bool,help="clean part_000")
    parser.add_argument("-t","--step", default=3, type=int,\
                        help="step in hour between two part files. Default = 3")
    #parser.add_argument("-k","--diffus",type=str,choices=['0','01','1','001'],\
    #    help='diffusivity parameter')
    parser.add_argument("-v","--vshift",type=int,choices=[0,10], default=10,\
                        help='vertical shift. Default = 10')
    #parser.add_argument("-hm","--hmax",type=int, default = 5640,\
    parser.add_argument("-hm","--hmax",type=int, default = 99,\
                        help='maximum considered integration time')
    # hmax =  1800 #: 75 days (75 * 24 = 1800)
    #5640 : 235 days (235 * 24 = 5640)
    # Boolean type arguments
    parser.add_argument("-q","--quiet",action=argparse.BooleanOptionalAction,\
        default=False,help="quiet with --quiet and --no-quiet for screen output")
    parser.add_argument("-b","--bak", action=argparse.BooleanOptionalAction,\
        default=True,\
        help="Create backup files with --bak and avoid with --no-bak. Default --bak")
    parser.add_argument("-c","--cont",action=argparse.BooleanOptionalAction,\
        default=False,help = "Continue from backup with --cont. Default --no-cont")
    parser.add_argument("-bs", "--bakstep", type=int, default=39,\
        help='backup step (hour) to be multiply with step. Default = 39')

    args = parser.parse_args()

    # to be updated
    # Define main directories
    if 'spirit' in socket.gethostname():
        #main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
            #SVC_Dir = '/bdd/CFMIP/SEL2'
        datPath = os.environ['HOME'].replace('/home/', '/data/')
        ekjDir = '/proju/flexpart/flexpart_in/EKJ/ejohansson'
        traj_dir = f"{ekjDir}/flexout/STC/Calipso"
        temp_dir = f"{traj_dir}/TempFiles"
        out_dir = "/data/legras/STC/Calipso-OUT"
        #out_dir = '/data/legras/STC'
    elif 'Mentat' in socket.gethostname():
        datPath = "C:\\cygwin64\\home\\berna\\data\\erikj\\flexout\\STC"
        traj_dir= os.path.join(datPath,'Calipso')
        temp_dir = os.path.join(traj_dir,'TempFiles')
        out_dir = os.path.join(datPath,'Calipso-OUT-NN')
    else:
        print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')
        exit()

    """ Parameters """
    # Parsing
    year = args.year
    month = args.month
    hmax = args.hmax
    quiet = args.quiet
    vshift = args.vshift
    step = args.step
    backup = args.bak
    restart = args.cont

    # Fixed parameters
    # max_life_time is imported from createFlexpartFiles
    # not necessarily a good idea
    age_bound = max_life_time - 0.25 #199.75
    advect = 'EAD'
    suffix =''
    clean0 = True
    diffus = '0'
    super =''

    # Derived parameters
    backup_step = args.bakstep * step
    dstep = timedelta(hours=step)
    dtRange = timedelta(hours=step)
    # time width of the parcel slice
    slice_width = timedelta(minutes=5)
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)


    """ date_end is not really the end of this run. It defines
    the month for which the analysed trajectories have been
    calculated and is the first day at 0h of this month. It is
    only used to define sdate and to name the input and output
    files. The name of this variable was recycled from a previous
    version where there was a beginning and an ending date. It
    could be changed to something more meaningful."""
    date_end = datetime(year=year, month=month, day=1, hour=0)
    """ Initial time is first day at 0h of the next month
    assuming date_end is the first day at 0h of the month.
    It is checked below that this corresponds to the stamp_date of the
    traczilla run and the run aborts if not."""
    sdate = date_end + relativedelta(months=1)
    # patch to fix the problem of the data hole on the evening of 2 August 2017
    # document this problem !! kept here as a post-it since not relevant
    #if sdate == datetime(2017,8,2): sdate = datetime(2017,8,3,6)

    if diffus == '0':
        diffus = ''
    else:
        diffus = '-D' + diffus

    if vshift > 0:
        super = '-super'+str(vshift)
    # Update the out_dir with the cloud type and the super paramater
    #out_dir = os.path.join(out_dir,'SVC-OUT-GridSat'+super)
    if vshift == 0:
        out_dir = os.path.join(out_dir,'Vshift_0')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.isdir(out_dir + '/out'):
        os.makedirs(out_dir + '/out')
    else:
        print('out_dir directory already created')

    #: filename of outfiles
    outnames = f"CALIOP-{advect}-{date_end.strftime('%Y%m')}{diffus}-{args.night}"
    #if args.use_dardar:
    outnames = f"{outnames}-DD"
    # Manage the file that receives the print output and start diverting print
    if quiet:
        # Output file
        print_file = os.path.join(out_dir,'out',outnames+'.out')
        if restart:
            fsock = open(print_file,'a')
        else:
            fsock = open(print_file,'w')

        sys.stdout=fsock

    # print all the arguments
    print('\nRun convsrcErikFullGridSatTropoExt performed on '+\
          socket.gethostname())
    if restart: print('RESTART RUN')
    else: print('PRISTINE RUN')

    print(datetime.now())
    print('\n Arguments \n')
    print('year    ',args.year)
    print('month   ',args.month)
    print('night   ',args.night)
    print('vshift  ',args.vshift)
    print('bak     ',args.bak)
    print('cont    ',args.cont)
    print('bakstep ',args.bakstep,'\n')
    print(f"age_bound = {age_bound}")
    print(f"Out_dir = {out_dir}\n")
    print('outnames ',outnames)

    # Directories of the backward trajectories and name of the output file
    ftraj = os.path.join(traj_dir,outnames)
    out_file2 = os.path.join(out_dir,outnames+'.h5')
    # backup file
    bak_file_prod0 = os.path.join(out_dir,outnames + '-K-backup-prod0.h5')
    bak_file_params = os.path.join(out_dir,outnames +'-K-backup-params.pkl')

    """ Initialization of the calculation """
    #: Control hmax and step to fit the content of the directory
    # Add a switch to avoid exit if hmax is meant to be smaller than the date
    # of last part file or change this as a warning
    emax, estep = findLastFile(ftraj)
    if (emax != hmax) or (estep != step):
        print('ERROR: There is a mismatch')
        print('hmax = %d' %hmax)
        print('emax = %d' %emax)
        print('step = %d' %step)
        print('extep = %d' %estep)
        print('run abort')
        sys.exit()
    # Initialize the dictionary of the parcel dictionaries
    partStep={}
    ''' Initialize a the GridSat grid parameters which are used in exiter
        to localize the parcel in the GridSat grid. This grid runs from -180
        in longitude and -70 in latitude.'''
    gg = geosat.GeoGrid('GridSat')

    # Read the index file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'part_000'),quiet=False)
    # Reset the longitude in a -180 to 180 range
    part0['x'] = np.where(part0['x']>180, part0['x'] - 360, part0['x'])

    # check that the stamp_date in the part0 file is the same as sdate and exit
    # if not as this is a serious mismatch that hinders the validity of this run
    if datetime.strptime(str(part0['stamp_date']),'%Y%m%d%H%M%S') != sdate:
        print('ERROR: the sdate and stamp_date are not identical')
        print('sdate      ',sdate)
        print('stamp_date ',str(part0['stamp_date']))
        print('run abort')
        sys.exit()
    current_date = sdate
    # check flag is clean, the result should be (0, 0, 0)
    print('flag is checked ',\
          ((part0['flag']&I_HIT)!=0).sum(), ((part0['flag']&I_DEAD)!=0).sum(),\
          ((part0['flag']&I_CROSSED)!=0).sum() )
    # check idx_orgn
    if part0['idx_orgn'] != 1:
        print('IDX_ORGN NOT 1 AS ASSUMED, CORRECTED WITH READ VALUE')
        print('VALUE ',part0['idx_orgn'])
        IDX_ORGN = part0['idx_orgn']

    # Open the init file produced by Erik
    InitFile = os.path.join(temp_dir,f"{outnames}-init.h5")
    init = Dataset(InitFile)
    # Read the iwc from DARDAR data and the initial saturation ratio
    part0['iwc'] = init.variables['iwc0'][:].astype(np.float32)
    part0['rvs'] = init.variables['rvs0'][:].astype(np.float32)
    init.close()

    # Get here the SVC mask for the launched particles
    SVCFile = os.path.join(temp_dir,\
              'Hist',f"svc_mask-global-{date_end.strftime('%Y%m')}-n-DD.h5")
    SVC = Dataset(SVCFile)
    part0['svc_mask'] = SVC.variables['svc_cld'][:].astype(np.bool)
    SVC.close()

    # Build a dictionary to host the results and initialize its numerous arrays
    # The numerous nan that will be left at the end are compressed in the output
    # file. This is easier than handling explicitely the spareness.
    prod0 = defaultdict(dict)
    prod0['attr']['numpart'] = part0['numpart']
    prod0['diag'] = {}
    # Information about the source, the saturation events and the clearings
    qoff = {'C06':0.6, 'C08':0.8, 'C10':1., 'C12':1.2, 'C14': 1.4}
    for off in qoff.keys():
        prod0['first_clear'][off] ={}
    for var in ('x','y','p','t','age','rvs'):
        prod0['src'][var] = np.full(part0['numpart'],np.nan,dtype='float32')
        prod0['lowest_sat'][var] = np.full(part0['numpart'],np.nan,dtype='float32')
        prod0['highest_sat'][var] = np.full(part0['numpart'],np.nan,dtype='float32')
        for off in qoff.keys():
            prod0['first_clear'][off][var] = np.full(part0['numpart'],np.nan,dtype='float32')
    # Set rvs at initial values for lowest sat and highest sat
    prod0['lowest_sat']['rvs'] = np.copy(part0['rvs'])
    prod0['highest_sat']['rvs'] = np.copy(part0['rvs'])

    # flag: first copy from part and add new flags
    prod0['flag_source'] = np.copy(part0['flag'])
    # Add flags according to the initial situation of clouds
    sel_cloud = part0['iwc']>0
    sel_svc = part0['svc_mask']>0
    prod0['flag_source'][sel_cloud] |= I_CLOUD
    prod0['flag_source'][sel_svc] |= I_SVC
    sel_puresvc = (sel_svc) & (~sel_cloud)
    # select clear sky at the origin of the parcels
    sel_clear = (~sel_svc) & (~sel_cloud)
    if np.sum(sel_cloud)+np.sum(sel_puresvc)+np.sum(sel_clear)!=prod0['attr']['numpart']:
        print('mismatch in the cloud mask')
        print('FIX IT')
    prod0['attr']['nice'] = int(np.sum(sel_cloud))
    prod0['attr']['nsvc'] = int(np.sum(sel_svc))
    prod0['attr']['npuresvc'] = int(np.sum(sel_puresvc))
    prod0['attr']['nclear'] = int(np.sum(sel_clear))

    # Initial statistics
    print('\nInitial statistics\n')
    print('year',year,'month',month)
    print('advect',advect)
    print('suffix',suffix)
    print('number of particles  ',part0['numpart'])
    print('number of ice part.  ',prod0['attr']['nice'],\
          '{:.2f}%'.format(100*prod0['attr']['nice']/part0['numpart']))
    print('number of SVC part.  ',prod0['attr']['nsvc'],\
          '{:.2f}%'.format(100*prod0['attr']['nsvc']/part0['numpart']))
    print('number of pure SVC   ',prod0['attr']['npuresvc'],\
          '{:.2f}%'.format(100*prod0['attr']['npuresvc']/part0['numpart']))
    print('number of clear part.',prod0['attr']['nclear'],\
          '{:.2f}%'.format(100*prod0['attr']['nclear']/part0['numpart']))

    # Set q0 according to initial conditions
    # For cloud conditions with IWC,
    # IWV read from Init file, rvs too (checked that it matches a recalculation)
    part0['q0'] = np.full(part0['numpart'],0.,dtype='float32')
    part0['q0'][sel_cloud] = part0['iwc'][sel_cloud] + part0['rvs'][sel_cloud]
    part0['q0'][sel_puresvc] = part0['rvs'][sel_puresvc]
    part0['q0'][sel_clear] = np.nan
    # Arrays no longer needed
    del part0['svc_mask'],part0['iwc']

    # number of hists and exits
    prod0['attr']['nhits'] = 0
    prod0['diag']['nhitsl'] = []
    prod0['attr']['nexits'] = 0
    prod0['diag']['nexitsl'] = []
    prod0['attr']['nnew'] = 0
    prod0['diag']['nnewl'] = []
    prod0['attr']['nold'] = 0
    prod0['diag']['noldl'] = []
    # Distance ante-post diagnostics
    prod0['diag']['dist_mean'] = []
    prod0['diag']['dist_max'] = []
    prod0['diag']['dist_std'] = []
    offset = 0 # must be zero for a non restart run

    # initialize datsat to None to force first read
    datsat = None

    print('\nInitialization completed')

    # Restart
    if restart:
        #: Make sure that there are backupfiles
        if (not os.path.isfile(bak_file_params)) | (not os.path.isfile(bak_file_prod0)):
            print('Params ao backup file missing.')
            print('Are you sure about restarting this run?')
            print('run abort')
            sys.exit()
        print('restart run')
        try:
            with gzip.open(bak_file_params,'r') as handle:
                [[offset, current_date], sel_clear, new] = pickle.load(handle)
            prod0 = h5prod0_read(bak_file_prod0)
            # start_from_begin = False
        except:
            print('cannot load backup')
            print('run abort')
            sys.exit()
        # Restore the last partante file
        partStep[offset] = readpart107(offset,ftraj,quiet=True)
        # Unclear why this is needed
        if isinstance(partStep[offset]['x'], list):
            partStep[offset]['x'] = np.asarray(partStep[offset]['x'])
        # Restore x within [-180,180]
        partStep[offset]['x'][partStep[offset]['x']>180] -= 360
        # Initialize sat and ERA5 yield one step ahead as a precaution
        get_sat = read_sat(current_date + dstep,dtRange,pre=True,vshift=vshift)
        get_ERA5 = read_ERA5(current_date + dstep,dtRange,pre=True)

    # Pristine run
    else:
        # reassign the part_000 file as a copy
        # part0 still needed in the sequel should not disappear
        partStep[0] = part0.copy()
        # partStep[0] = readpart107(0,ftraj,quiet=False)
        # cleaning is necessary for traczilla runs starting in the fake restart
        # mode, which is now the most common,
        # otherwise all parcels are thought to exit at the first step
        if clean0:
            partStep[0]['idx_back']=[]
        # parcels with longitude east of zero degree are set to negative values
        partStep[0]['x'][partStep[0]['x']>180] -= 360
        # used to get non borne parcels
        # no parcel is initially new
        new = np.full(part0['numpart'],False,dtype='bool')
        # Build the satellite field generator
        get_sat = read_sat(sdate,dtRange,pre=True,vshift=vshift)
        # Build the ECMWF field generator
        # (both GridSat and ERA5 are available at the same dtRange = 6h interval
        # in this application)
        get_ERA5 = read_ERA5(sdate,dtRange,pre=True)

    """ Main loop on the output time steps """
    for hour in range(step + offset, hmax + 1, step):
        pid = os.getpid()
        py = psutil.Process(pid)
        memoryUse = py.memory_info()[0]/2**30
        print('memory use: {:4.2f} gb'.format(memoryUse))
        # Get rid of dictionary no longer used
        if hour >= 2*step:
            try:
                del partStep[hour-2*step]
            except:
                print('nothing to delete, probably due to restart')
        # Read the new data, to be assignes as partpost
        partStep[hour] = readpart107(hour,ftraj,quiet=True)
        # Link the names
        partante = partStep[hour-step]
        partpost = partStep[hour]
        # parcels with longitude east of zero degree are set to negative values
        # done with try to prevent errors when the file is empty
        try:
            partpost['x'][partpost['x']>180] -= 360
            #print("partpost x ",partpost['x'].min(),partpost['x'].max())
        except: pass

        if partpost['nact']>0:
            print('hour ',hour,'  numact ', partpost['nact'], '  max p ',partpost['p'].max())
        else:
            print('hour ',hour,'  numact ', partpost['nact'])
        # New date valid for partpost
        current_date -= dstep

        # Processing of water mixing ratios
        # Select non stopped parcels in partante
        if len(partante['idx_back']) >0:
            selec = (prod0['flag_source'][partante['idx_back']-IDX_ORGN] & I_DEAD) == 0
            # index of the active non stopped parcels in the full list
            idx0 = partante['idx_back'][selec]
            # index of the active non stopped parcels in the reduced list
            idxt = np.where(selec)[0]
            for i in range(len(idx0)):
                idloc = idxt[i]
                idglo = idx0[i]-IDX_ORGN
                rvs = esati_murphy(partante['p'][idloc],partante['t'][idloc])
                # tests for the min and the max rvs
                if rvs < prod0['lowest_sat']['rvs'][idglo]:
                    prod0['lowest_sat']['rvs'][idglo] = rvs
                    prod0['lowest_sat']['x'][idglo] = partante['x'][idloc]
                    prod0['lowest_sat']['y'][idglo] = partante['y'][idloc]
                    prod0['lowest_sat']['p'][idglo] = partante['p'][idloc]
                    prod0['lowest_sat']['t'][idglo] = partante['t'][idloc]
                    prod0['lowest_sat']['age'][idglo] = (part0['ir_start'][idglo] \
                                                        - partante['itime'])/86400
                if rvs > prod0['highest_sat']['rvs'][idglo]:
                    prod0['highest_sat']['rvs'][idglo] = rvs
                    prod0['highest_sat']['x'][idglo] = partante['x'][idloc]
                    prod0['highest_sat']['y'][idglo] = partante['y'][idloc]
                    prod0['highest_sat']['p'][idglo] = partante['p'][idloc]
                    prod0['highest_sat']['t'][idglo] = partante['t'][idloc]
                    prod0['highest_sat']['age'][idglo] =\
                        (part0['ir_start'][idglo] - partante['itime'])/86400
                """ test for the first encounter with clear air
                 tested in increasing order of the threshold
                 the lowest is the more easily satisfied
                 the flag for the a given threshold sets also the flag
                 for lower thresholds if not already set """
                flag = prod0['flag_source'][idglo]
                for off in qoff:
                    # Select cloudy parcels not DEAD and which have not
                    # already passed the threshold
                    if ((flag & I_CLEAR[off]) != I_CLEAR[off]) &\
                       ((flag & I_DEAD) == 0) & ~sel_clear[idglo]:
                        if rvs > qoff[off] * part0['q0'][idglo]:
                            prod0['flag_source'][idglo] |= I_CLEAR[off]
                            prod0['first_clear'][off]['x'][idglo] =\
                                partante['x'][idloc]
                            prod0['first_clear'][off]['y'][idglo] =\
                                partante['y'][idloc]
                            prod0['first_clear'][off]['p'][idglo] =\
                                partante['p'][idloc]
                            prod0['first_clear'][off]['t'][idglo] =\
                                partante['p'][idloc]
                            prod0['first_clear'][off]['rvs'][idglo] = rvs
                            prod0['first_clear'][off]['age'][idglo] =\
                                (part0['ir_start'][idglo] - partante['itime'])/86400

        """ Select the parcels that are common to the two steps
        ketp_a is a logical field with same length as partante
        kept_p is a logical field with same length as partpost
        new diagnostic simplified in this application
        Once the beginning of investigated month is reached, there
        should not be any new parcel.
        """
        kept_a = np.isin(partante['idx_back'],partpost['idx_back'],assume_unique=True)
        kept_p = np.isin(partpost['idx_back'],partante['idx_back'],assume_unique=True)
        # in1d deprecated
        #new_p = ~np.in1d(partpost['idx_back'],partpost['idx_back'],assume_unique=True)
        new_p = len(partpost['x'])-kept_p.sum()
        print('kept a, p ',len(kept_a),len(kept_p),kept_a.sum(),kept_p.sum(),
              '  new ',len(partpost['x'])-kept_p.sum())
        """ The first two items are the number of parcels at ante and post time
         The third and fourth are the number of common parcesl. The two numbers
         should be identical.
         The fifth is the number of new parcels in post. Should fall to 0 outside
         the time interval of parcel craetion. """

        # new counters below

        """ IDENTIFY AND TAKE CARE OF DEADBORNE AS NON BORNE PARCELS
        """
        # Flag as new all the created parcels which are seen at least once
        if (hour <= 31*24) & (partpost['nact']>0):
            # Flag new particles created since last output and alive
            new[partpost['idx_back'][~kept_p]-IDX_ORGN] = True
            # increment the counter of new particles using new_p above
            prod0['attr']['nnew'] += new_p
            prod0['diag']['nnewl'].append(new_p)
            # after the end of the parcel creation interval, there is no need
            # to update these counters

        """ Check the missing parcels after the end of the creation period
            and flag them as deadborne with 0 age and exit at their birthplace
            This is done once. It could be done at each step but is would be
            much more complicated."""
        if hour == 31*24 + step:
            prod0['attr']['ndborne'] = np.sum(~new)
            prod0['flag_source'][~new] |= I_DBORNE + I_DEAD
            prod0['src']['x'][~new] = part0['x'][~new]
            prod0['src']['y'][~new] = part0['y'][~new]
            prod0['src']['p'][~new] = part0['p'][~new]
            prod0['src']['t'][~new] = part0['t'][~new]
            prod0['src']['age'][~new] = 0
            # rvs here attributed but recalculated before final
            # output from p and t
            prod0['src']['rvs'][~new] = part0['rvs'][~new]
            print('number of dead borne',prod0['attr']['ndborne'],\
                  part0['numpart']-prod0['attr']['nnew'])
            ndborne_cloud = np.sum((~new) & sel_cloud)
            ndborne_puresvc = np.sum((~new) & sel_puresvc)
            ndborne_clear = np.sum((~new) & sel_clear)
            print('among which cloudy',ndborne_cloud,' svc',ndborne_puresvc,
                  ' clear',ndborne_clear)
            prod0['attr']['ndborne_cloud'] = ndborne_cloud
            prod0['attr']['ndborne_svc'] = np.sum((~new) & sel_svc)
            prod0['attr']['ndborne_puresvc'] = ndborne_puresvc
            prod0['attr']['ndborne_clear'] = ndborne_clear
            # Delete arrays no longer used
            del new
            del part0['rvs']
            new = None #: Needs to save something in Backupfiles

        """ PROCESSING OF CROSSED PARCELS
        This section processes the parcel found at ante time which have
        disappeared at post time and tries to find a reason.
        Parcels already DEAD for another reason are not labelled as
        CROSSED."""
        if len(~kept_a)>0:
            exits = exiter(int((partante['itime']+partpost['itime'])/2), \
                partante['x'][~kept_a],partante['y'][~kept_a],partante['p'][~kept_a],\
                partante['t'][~kept_a],partante['idx_back'][~kept_a],\
                prod0['flag_source'],prod0['src']['x'],prod0['src']['y'],\
                prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],\
                part0['ir_start'], gg.box_range)
            prod0['attr']['nexits'] += exits
            prod0['diag']['nexitsl'].append(exits)

            print('exit ',prod0['attr']['nexits'], exits, np.sum(~kept_a),\
                  len(kept_a) - len(kept_p))
            """ The first number is the total number number of crossed parcels.
            The second is the number of parcel crossed in this step. It should
            be inferior to the second number representing the number of missing
            ante parcels in post. The third number is the difference of number
            of parcel in ante and kept. In the initial period it should be <=
            the number of new parcels. At large time, it should be equal to the
            second.
            """
        else:
            prod0['diag']['nexitsl'].append(0)
            print('exit ',prod0['attr']['nexits'], 0, 0,\
                  len(kept_a) - len(kept_p))

        # Notice that rvs is not calculated inside this function and cannot be
        # if accelerated with numba. It is done before final output.

        """ PROCESSING OF PARCELS WHICH ARE COMMON TO THE TWO OUTPUTS """
        """ Select the kept parcels which have not been hit yet
            !!! Never use "and" between two lists, the result is wrong
            IMPORTANT WARNING: a parcel is live when it is present in the
            two outputs. Hocwever it can be declared DEAD later when it hits
            a cloud and is then no longer processed.
        """

        if kept_p.sum()==0:
            live_a = live_p = kept_p
        else:
            live_a = np.logical_and(kept_a,(prod0['flag_source'][partante['idx_back']-IDX_ORGN] & I_DEAD) == 0)
            live_p = np.logical_and(kept_p,(prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_DEAD) == 0)
        print('live a, b ',live_a.sum(),live_p.sum())

        """ Correction step that modifies partante['x'] to avoid big jumps
        at the periodicity frontier on the x-axis. """
        # Done this way as double slicing forbidden in assignation
        if live_p.sum()>0:
            xante = partante['x'][live_a]
            sel1 = (partpost['x'][live_p] - xante) < -180
            sel2 = (partpost['x'][live_p] - xante) > 180
            xante[sel1] += 360
            xante[sel2] -= 360
            partante['x'][live_a] = xante
            print("Date line crossing correction", np.sum(sel1)+np.sum(sel2))
            """ Produce here a few diagnostics on the distance between partante
                and partpost """
            ddd = spher_dist(partante['x'][live_a],partante['y'][live_a],\
                             partpost['x'][live_p],partpost['y'][live_p])
            ddd_mean = np.mean(ddd) # in (pseudo) radian (we do not do the arc sin in spher_dist)
            prod0['diag']['dist_mean'].append(ddd_mean/fff)
            prod0['diag']['dist_max'].append(ddd.max()/fff)
            prod0['diag']['dist_std'].append(np.std(ddd,mean=ddd_mean)/fff)
            print('ante-post (°) max {:4.2f} mean {:4.2f} std {:4.2f}'.\
                  format(prod0['diag']['dist_max'][-1],ddd_mean/fff,\
                         prod0['diag']['dist_std'][-1]))
        else:
            prod0['diag']['dist_mean'].append(0)
            prod0['diag']['dist_max'].append(0)
            prod0['diag']['dist_std'].append(0)

        gsp = get_slice_part(partante,partpost,live_a,live_p,current_date,dstep,slice_width)
        if verbose: print('built parcel generator for ',current_date)

        """  MAIN LOOP ON THE PARCEL TIME SLICES  """
        # counter over the loop
        nb_hits = 0
        for i in range(nb_slices):
            # get the slice for the particles
            datpart = next(gsp)
            # Check x within (-180,180), necessary when GridSat is used
            #if datpart['x'] != []:
            if len(datpart['x'])>0:
                #: datpart['x'] = np.where(datpart['x']>180, datpart['x'] - 360, datpart['x'])
                datpart['x'] = (datpart['x']+180)%360 - 180
            if verbose: print('part slice ',i, datpart['time'])
            # Check whether the present satellite image is valid
            # The while should ensure that the run synchronizes
            # when it starts.

            swap = False
            while check(datsat, datpart['time']) is False:
                # if not get next satellite image
                datsat = next(get_sat)
                datera = next(get_ERA5)
                swap = True

            """ Select the not DEAD parcels located within the domain """
            # The domain is here based on the GridSat grid which combines GEO and LEO data
            # at high latitude with filling methods. Although the grid goes to 70 degrees
            # we limit the domain here to 65°
            # TODO: remove this hard coding
            if len(datpart['x'])>0 :
                # reset longitudes in the -180,180 domain (for slices coming
                # from pairs crossing the dateline)
                datpart['x'][datpart['x']>180] -= 360
                datpart['x'][datpart['x']<-180] += 360
                indomain = np.all((datpart['y']>-65,datpart['y']<65),axis=0)
                indomain = indomain & ((prod0['flag_source'][datpart['idx_back']-IDX_ORGN] & I_DEAD)==0)
            else:
                indomain = np.zeros(1)

            """ PROCESS THE COMPARISON OF PARCEL PRESSURES TO CLOUDS """
            if indomain.sum()>0:
                # tropopause pressure at the location of the parcels
                #print('ptropo x ',datpart['x'].min(),datpart['x'].max())
                #print('ptropo y ',datpart['y'].min(),datpart['y'].max())
                ''' The number of hits is collected in a list to check for
                anomalies that might point to corrupted gridsat images '''
                ptrop = datera.fP(np.transpose([datpart['y'][indomain],\
                                                datpart['x'][indomain]]))
                new_hits = convbirth(datpart['itime'],
                    datpart['x'][indomain],datpart['y'][indomain],datpart['p'][indomain],\
                    datpart['t'][indomain],datpart['idx_back'][indomain],\
                    prod0['flag_source'],ptrop,prod0['src']['x'],prod0['src']['y'],\
                    prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],\
                    datsat.var['IR0'], part0['ir_start'],\
                    datsat.geogrid.box_range[0,0],datsat.geogrid.box_range[1,0],\
                    datsat.geogrid.stepx,datsat.geogrid.stepy,\
                    datsat.geogrid.box_binx,datsat.geogrid.box_biny)
            else: new_hits = 0
            if swap:
                print('post swap hits', new_hits)
                swap_hits = new_hits
                swap = False
            nb_hits += new_hits

        if nb_hits>0:
            prod0['diag']['nhitsl'].append(nb_hits)
            prod0['attr']['nhits'] += nb_hits
            print('total hits {} prop swap (%) {:.2f}'\
                  .format(nb_hits,100*swap_hits/nb_hits))
        else: print('total hits',nb_hits)
        """ End of of loop on slices """

        """ INSERT HERE CODE FOR PARCELS ENDING BY AGE -> TO BE CHECKED
         Important: if the age limit is set in the run, the age_bound needs to
         be decremented by at least one output step. Otherwise, the parcel disappear
         from the output before being processed here and is wrongly processed
         as crossed."""
        # Check the age limit (easier to do it here)
        # Revert partante correction        # Revert partante correction

        if len(partante['idx_back']) > 0:
            if np.sum(live_a)>0:
                xante[sel1] -= 360
                xante[sel2] += 360
                partante['x'][live_a] = xante
                del xante
            print("Manage age limit",flush=True)
            # select here the non DEAD parcels which have reached the age limit
            age_sec =\
                part0['ir_start'][partante['idx_back']-IDX_ORGN]-partante['itime']
            IIold_o = age_sec > (age_bound-(step/24)) * 86400
            IIold_o = IIold_o &\
                ((prod0['flag_source'][partante['idx_back']-IDX_ORGN] & I_DEAD)==0)
            idx_IIold = partante['idx_back'][IIold_o]
            j_IIold_o = np.where(IIold_o)
            # set the flag
            prod0['flag_source'][idx_IIold-IDX_ORGN] =\
                prod0['flag_source'][idx_IIold-IDX_ORGN] | I_DEAD+I_OLD
            prod0['src']['x'][idx_IIold-IDX_ORGN] = partante['x'][j_IIold_o]
            prod0['src']['y'][idx_IIold-IDX_ORGN] = partante['y'][j_IIold_o]
            prod0['src']['p'][idx_IIold-IDX_ORGN] = partante['p'][j_IIold_o]
            prod0['src']['t'][idx_IIold-IDX_ORGN] = partante['t'][j_IIold_o]
            prod0['src']['age'][idx_IIold-IDX_ORGN] = \
                (part0['ir_start'][idx_IIold-IDX_ORGN]- partante['itime'])/86400
            # rvs calculated from t and p before final output
            print("number of IIold ",len(idx_IIold))
            prod0['attr']['nold'] += len(idx_IIold)
            prod0['diag']['noldl'].append(len(idx_IIold))
            # if len(prod0['src']['x']) > 0:
            #     if ((prod0['src']['x'] > 180).any() or (prod0['src']['x']<-180).any()):
            #         print('Check x')
            #         pdb.set_trace()
        # find parcels still alive       if kept_p.sum()==0:
        try:
            prod0['attr']['nlive'] =\
                ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_DEAD) == 0).sum()
            prod0['attr']['n_nohit'] =\
                ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT) == 0).sum()
        except:
            prod0['attr']['nlive'] = 0
            prod0['attr']['n_nohit'] = 0
        print('end hour ',hour,'  numact', partpost['nact'],\
              ' nexits',prod0['attr']['nexits'],' nhits',prod0['attr']['nhits'],\
              ' nlive',prod0['attr']['nlive'],' nohit',prod0['attr']['n_nohit'],\
              ' nold',prod0['attr']['nold'])
        # check that nlive + nhits + nexits = numpart, should be true after
        # the first day
        # does not make sense before ndborne has been calculated
        if hour >= 31*24 + step:
            # updated to correct spurious warnings
            # ndprobe which is not restablished after restart replaced by prod0['attr']['ndborne']
            # nhold added in the sum
            if part0['numpart'] != prod0['attr']['nexits']\
                                 + prod0['attr']['nhits']\
                                 + prod0['attr']['nold']\
                                 + prod0['attr']['nlive']\
                                 + prod0['attr']['ndborne']:
                print('@@@ ACHTUNG numpart not equal to sum ',part0['numpart'],\
                      prod0['attr']['nexits']+prod0['attr']['nhits']\
                      +prod0['attr']['nlive']+prod0['attr']['ndborne'])

        sys.stdout.flush()

        if backup & (hour % backup_step == 0):
            print("Save backup")
            print(bak_file_prod0)
            tic = time.time()
            h5prod0_write(bak_file_prod0)
            with gzip.open(bak_file_params,'w') as handle:
                pickle.dump([[hour, current_date], sel_clear, new],\
                            handle, protocol=pickle.HIGHEST_PROTOCOL)
            #fk.save(bak_file_prod0, (prod0,part0['q0'],sel_clear,new))
            #fk.save(bak_file_params,{'params': [hour, nhits, nexits, nold, nnew, current_date]})
            toc = time.time()
            print('It takes %d s to save backupfiles' %(int(toc - tic)))
    """ End of the procedure and storage of the result """

    # Calculation of rvs at the final location of each parcel
    prod0['src']['rvs'] = esati_murphy(prod0['src']['p'], prod0['src']['t'])
    #output file
    # prod0['src']['x'] = np.where(prod0['src']['x'] > 180, prod0['src']['x'] - 360, prod0['src']['x'])
    print('Max and Min of x')
    print(np.nanmax(prod0['src']['x']))
    print(np.nanmin(prod0['src']['x']))
    #fk.save(out_file2,prod0,compression='zlib')
    h5prod0_write(out_file2)
    try:
        os.remove(bak_file_prod0)
        os.remove(bak_file_params)
    except:pass
    print('CONGRATULATION, ALL THE JOB IS DONE')
    #pickle.dump(prod0,gzip.open(out_file,'wb'))
    # close the print file
    if quiet:
        sys.stdout.flush()
        fsock.close()

"""@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""
