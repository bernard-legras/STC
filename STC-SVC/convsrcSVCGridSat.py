# -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the air sampled during the StratoClim
campaign.
This version is based on the GridSat product and compares the temperature of the parcel to the IR brightness temperature.

Notice that the part time is supposed to start at day+1 0h where day is the day of the flight.

Created on Wednesday 6 June 2019

@author: Bernard Legras
"""

import socket
import numpy as np
from collections import defaultdict
from numba import jit
from datetime import datetime, timedelta
import os
import sys
import argparse
import psutil
import deepdish as dd
#import SAFNWCnc
import geosat
#import constants as cst
from io107 import readpart107, readidx107

p0 = 100000.
I_DEAD = 0x200000
I_HIT = 0x400000
I_CROSSED = 0x2000000
I_DBORNE =  0x1000000
# ACHTUNG I_DBORNE has been set to 0x10000000 (one 0 more) in a number of earlier analysis 
# prior to 18 March 2018

# misc parameters
# step in the cloudtop procedure
cloudtop_step = timedelta(hours=12)
# low p cut in the M55 traczilla runs
lowpcut = 3000
# highpcut in the M55 traczilla runs
highpcut = 50000

# if True print a lot oj junk
verbose = False
debug = False

# idx_orgn was not set to 1 but to 0 in M55 and GLO runs
IDX_ORGN = 0

# Error handling
class BlacklistError(Exception):
    pass

blacklist = []

#%%
"""@@@@@@@@@@@@@@@@@@@@@@@@   MAIN   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

def main():
    global IDX_ORGN
    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-d1","--day1",type=int,choices=1+np.arange(31),help="day1")
    parser.add_argument("-d2","--day2",type=int,choices=1+np.arange(31),help="day2")
    parser.add_argument("-a","--advect",type=str,choices=["OPZ","EAD","EAZ","EID","EIZ"],help="source of advecting winds")
    parser.add_argument("-s","--suffix",type=str,help="suffix for special cases")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
    #parser.add_argument("-c","--clean0",type=bool,help="clean part_000")
    parser.add_argument("-t","--step",type=int,help="step in hour between two part files")
    parser.add_argument("-ct","--cloud_type",type=str,choices=["meanhigh","veryhigh","silviahigh"],help="cloud type filter")
    parser.add_argument("-k","--diffus",type=str,choices=['01','1','001'],help='diffusivity parameter')
    parser.add_argument("-v","--vshift",type=int,choices=[0,10],help='vertical shift')
    parser.add_argument("-hm","--hmax",type=int,help='maximum considered integration time')
    
    # to be updated
    # Define main directories
    if 'ciclad' in socket.gethostname():
            #main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
            #SVC_Dir = '/bdd/CFMIP/SEL2'
            traj_dir = '/data/akottayil/flexout/STC/BACK'
            out_dir = '/data/akottayil/STC'
    else:
         print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')
         exit()   

    """ Parameters """
    step = 3
    hmax = 732
    dstep = timedelta (hours=step)
    # time width of the parcel slice
    slice_width = timedelta(minutes=5)
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)
    # default values of parameters
    # date of the flight
    year=2017
    month=7
    day1=1
    day2=10
    advect = 'EAD'
    suffix =''
    quiet = False
    clean0 = True
    cloud_type = 'silviahigh'
    diffus = '01'
    vshift = 0
    super =''
    dtRange = timedelta(hours=6)
    args = parser.parse_args()
    if args.year is not None: year=args.year
    if args.month is not None: month=args.month
    if args.advect is not None: advect=args.advect
    if args.day1 is not None: day1 = args.day1
    if args.day2 is not None: day2 = args.day2
    if args.hmax is not None: hmax = args.hmax
    if args.suffix is not None: suffix='-'+args.suffix
    if args.quiet is not None:
        if args.quiet=='y': quiet=True
        else: quiet=False
    #if args.clean0 is not None: clean0 = args.clean0
    if args.cloud_type is not None: cloud_type = args.cloud_type
    if args.step is not None: step = args.step
    if args.diffus is not None: diffus = args.diffus
    diffus = '-D' + diffus
    if args.vshift is not None: 
        vshift = args.vshift
        if vshift > 0:
            super = 'super'+str(vshift)

    # Update the out_dir with the cloud type and the super paramater
    out_dir = os.path.join(out_dir,'STC-BACK-OUT-GridSat-'+super+cloud_type)
    try:
        os.mkdir(out_dir)
        os.mkdir(out_dir+'/out')
    except:
        print('out_dir directory already created')    
    
    # Dates beginning and end
    date_beg = datetime(year=year, month=month, day=day1, hour=0)
    date_end = datetime(year=year, month=month, day=day2, hour=0)

    # Manage the file that receives the print output
    if quiet:
        # Output file      
        print_file = os.path.join(out_dir,'out','BACK-SVC-EAD-'+date_beg.strftime('%b-%Y-day%d-')+date_end.strftime('%d-D01')+'.out')
        fsock = open(print_file,'w')
        sys.stdout=fsock

    # initial time to read the sat files
    # should be after the end of the flight
    sdate = date_end + timedelta(days=1)
    print('year',year,'month',month,'day',day2)
    print('advect',advect)
    print('suffix',suffix)
    
    # patch to fix the problem of the data hole on the evening of 2 August 2017
    if sdate == datetime(2017,8,2): sdate = datetime(2017,8,3,6)

    # Directories of the backward trajectories and name of the output file
    ftraj = os.path.join(traj_dir,'BACK-SVC-EAD-'+date_beg.strftime('%b-%Y-day%d-')+date_end.strftime('%d-D01'))    
    out_file2 = os.path.join(out_dir,'BACK-SVC-EAD-'+date_beg.strftime('%b-%Y-day%d-')+date_end.strftime('%d-D01')+'.hdf5')

    """ Initialization of the calculation """
    # Initialize the dictionary of the parcel dictionaries
    partStep={}
    # initialize a grid that will be used to before actually doing any read
    gg = geosat.GeoGrid('GridSat')

    # Build the satellite field generator
    get_sat = read_sat(sdate,dtRange,pre=True,vshift=vshift)
    # Read the index file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'part_000'),quiet=True)
    print('numpart',part0['numpart'])
    # stamp_date not set in these runs
    # current_date actually shifted by one day / sdate
    # We assume here that part time is defined from this day at 0h
    # sdate defined above should be equal or posterior to current_date
    current_date = sdate
    # check flag is clean
    print('check flag is clean ',((part0['flag']&I_HIT)!=0).sum(),((part0['flag']&I_DEAD)!=0).sum(),\
                                 ((part0['flag']&I_CROSSED)!=0).sum())
    # check idx_orgn
    if part0['idx_orgn'] != 0:
        print('MINCHIA, IDX_ORGN NOT 0 ASSTC-M55-OUT-SAF-super1silviahigh ASSUMED, CORRECTED WITH READ VALUE')
        print('VALUE ',part0['idx_orgn'])
        IDX_ORGN = part0['idx_orgn']

    # Build a dictionary to host the results
    prod0 = defaultdict(dict)
    prod0['src']['x'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['y'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['p'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['t'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['age'] = np.empty(part0['numpart'],dtype='int')
    prod0['flag_source'] = part0['flag']
    # truncate eventually to 32 bits at the output stage

    # read the part_000 file
    partStep[0] = readpart107(0,ftraj,quiet=True)
    # cleaning is necessary for runs starting in the fake restart mode
    # otherwise all parcels are thought to exit at the first step
    if clean0:
        partStep[0]['idx_back']=[]

    # number of hists and exits
    nhits = 0
    nexits = 0
    ndborne = 0
    nnew = 0
    
    # initialize datsat to None to force first read
    datsat = None

    # used to get non borne parcels
    new = np.empty(part0['numpart'],dtype='bool')
    new.fill(False)

    print('Initialization completed')

    """ Main loop on the output time steps """
    for hour in range(step,hmax+1,step):
        pid = os.getpid()
        py = psutil.Process(pid)
        memoryUse = py.memory_info()[0]/2**30
        print('memory use: {:4.2f} gb'.format(memoryUse))
        # Get rid of dictionary no longer used
        if hour >= 2*step: del partStep[hour-2*step]
        # Read the new data
        partStep[hour] = readpart107(hour,ftraj,quiet=True)
        # Link the names
        partante = partStep[hour-step]
        partpost = partStep[hour]
        if partpost['nact']>0:
            print('hour ',hour,'  numact ', partpost['nact'], '  max p ',partpost['p'].max())
        else:
            print('hour ',hour,'  numact ', partpost['nact'])
        # New date valid for partpost
        current_date -= dstep
        """ Select the parcels that are common to the two steps
        ketp_a is a logical field with same length as partante
        kept_p is a logical field with same length as partpost
        After the launch of the earliest parcel along the flight track, there
        should not be any member in new
        The parcels
        """
        kept_a = np.in1d(partante['idx_back'],partpost['idx_back'],assume_unique=True)
        kept_p = np.in1d(partpost['idx_back'],partante['idx_back'],assume_unique=True)
        #new_p = ~np.in1d(partpost['idx_back'],partpost['idx_back'],assume_unique=True)
        print('kept a, p ',len(kept_a),len(kept_p),kept_a.sum(),kept_p.sum(),'  new ',len(partpost['x'])-kept_p.sum())

        """ IDENTIFY AND TAKE CARE OF DEADBORNE AS NON BORNE PARCELS """
        if (hour < ((day2-day1+5)*24)) & (partpost['nact']>0):
            new[partpost['idx_back'][~kept_p]-IDX_ORGN] = True
            nnew += len(partpost['x'])-kept_p.sum()
        if hour == ((day2-day1+5)*24):
            ndborne = np.sum(~new)
            prod0['flag_source'][~new] |= I_DBORNE + I_DEAD
            prod0['src']['x'][~new] = part0['x'][~new]
            prod0['src']['y'][~new] = part0['y'][~new]
            prod0['src']['p'][~new] = part0['p'][~new]
            prod0['src']['t'][~new] = part0['t'][~new]
            prod0['src']['age'][~new] = 0
            print('number of dead borne',ndborne,part0['numpart']-nnew)
            del new

        """ INSERT HERE CODE FOR NEW PARCELS """
        # nothing to be done for new parcels, just wait and see

        """ PROCESSING OF CROSSED PARCELS """
        if len(kept_a)>0:
            exits = exiter(int((partante['itime']+partpost['itime'])/2), \
                partante['x'][~kept_a],partante['y'][~kept_a],partante['p'][~kept_a],\
                partante['t'][~kept_a],partante['idx_back'][~kept_a],\
                prod0['flag_source'],prod0['src']['x'],prod0['src']['y'],\
                prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],\
                part0['ir_start'], gg.box_range)
            nexits += exits
            print('exit ',nexits, exits, np.sum(~kept_a), len(kept_a) - len(kept_p))

        """ PROCESSING OF PARCELS WHICH ARE COMMON TO THE TWO OUTPUTS  """
        # Select the kept parcels which have not been hit yet
        # !!! Never use and between two lists, the result is wrong

        if kept_p.sum()==0:
            live_a = live_p = kept_p
        else:
            live_a = np.logical_and(kept_a,(prod0['flag_source'][partante['idx_back']-IDX_ORGN] & I_DEAD) == 0)
            live_p = np.logical_and(kept_p,(prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_DEAD) == 0)
        print('live a, b ',live_a.sum(),live_p.sum())

        # Build generator for parcel locations of the 5' slices
        gsp = get_slice_part(partante,partpost,live_a,live_p,current_date,dstep,slice_width)
        if verbose: print('built parcel generator for ',current_date)

        """  MAIN LOOP ON THE PARCEL TIME SLICES  """

        for i in range(nb_slices):
            # get the slice for the particles
            datpart = next(gsp)
            if verbose: print('part slice ',i, datpart['time'])
            # Check whether the present satellite image is valid
            # The while should ensure that the run synchronizes
            # when it starts.
            while check(datsat, datpart['time']) is False:
                # if not get next satellite image 
                datsat = next(get_sat)
            
            """ PROCESS THE COMPARISON OF PARCEL PRESSURES TO CLOUDS """
            if len(datpart['x'])>0:
                nhits += convbirth(datpart['itime'],
                    datpart['x'],datpart['y'],datpart['p'],datpart['t'],datpart['idx_back'],\
                    prod0['flag_source'],prod0['src']['x'],prod0['src']['y'],\
                    prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],\
                    datsat.var['IR0'], part0['ir_start'],\
                    datsat.geogrid.box_range[0,0],datsat.geogrid.box_range[1,0],datsat.geogrid.stepx,\
                    datsat.geogrid.stepy,datsat.geogrid.box_binx,datsat.geogrid.box_biny)

            sys.stdout.flush()

        """ End of of loop on slices """
        # find parcels still alive       if kept_p.sum()==0:
        try:
            nlive = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_DEAD) == 0).sum()
            n_nohit = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT) == 0).sum()
        except:
            nlive = 0
            n_nohit =0
        print('end hour ',hour,'  numact', partpost['nact'], ' nexits',nexits,' nhits',nhits, ' nlive',nlive,' nohit',n_nohit)
        # check that nlive + nhits + nexits = numpart, should be true after the first day
        if part0['numpart'] != nexits + nhits + nlive + ndborne:
            print('@@@ ACHTUNG numpart not equal to sum ',part0['numpart'],nexits+nhits+nlive+ndborne)

    """ End of the procedure and storage of the result """
    #output file
    dd.io.save(out_file2,prod0,compression='zlib')
    #pickle.dump(prod0,gzip.open(out_file,'wb'))
    # close the print file
    if quiet: fsock.close()

"""@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

#%%
""" Functions related to the parcel data """

def get_slice_part(part_a,part_p,live_a,live_p,current_date,dstep,slice_width):
    """ Generator to generate 5' slices along flight track """
    nb_slices = int(dstep/slice_width)
    ta = current_date + dstep
    tp = current_date
    tf = ta
    empty_live = (live_a.sum() == 0)
    for i in range(nb_slices):
        ti = tf- slice_width
        # note that 0.5*(ti+tf) cannot be calculated as we are adding two dates
        tmid = ti+0.5*(tf-ti)
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
""" Function managing the exiting parcels """

@jit(nopython=True)
def exiter(itime, x,y,p,t,idx_back, flag,xc,yc,pc,tc,age, ir_start, rr):
    nexits = 0
    for i in range(len(x)):
        i0 = idx_back[i]-IDX_ORGN
        if flag[i0] & I_DEAD == 0:
            nexits += 1
            xc[i0] = x[i]
            yc[i0] = y[i]
            tc[i0] = t[i]
            pc[i0] = p[i]
            age[i0] = ir_start[i0] - itime
            if   y[i] < rr[1,0] + 4.: excode = 6
            elif x[i] < rr[0,0] + 4.: excode = 3
            elif y[i] > rr[1,1] - 4.: excode = 4
            elif x[i] > rr[0,1] - 4.: excode = 5
            elif p[i] > highpcut - 150: excode = 1
            elif p[i] < lowpcut  + 15 : excode = 2
            else:                   excode = 7
            flag[i0] |= (excode << 13) + I_DEAD + I_CROSSED
    return nexits

#%%
""" Function doing the comparison between parcels and clouds and setting the result field 
In this version, we check only whether the temperature is warmer than the brightness temperature
from the infrared window. """

@jit(nopython=True)
def convbirth(itime, x,y,p,t,idx_back, flag,xc,yc,pc,tc,age, BT, ir_start, x0,y0,stepx,stepy,binx,biny):
    nhits = 0
    for i in range(len(x)):
        idx = min(int(np.floor((x[i]-x0)/stepx)),binx-1)
        idy = min(int(np.floor((y[i]-y0)/stepy)),biny-1)
        if BT[idy,idx] < tc[i]:
            i0 = idx_back[i]-IDX_ORGN
            if flag[i0] & I_DEAD == 0:
                nhits += 1
                flag[i0] |= I_HIT + I_DEAD
                xc[i0] = x[i]
                yc[i0] = y[i]
                tc[i0] = t[i]
                pc[i0] = p[i]
                age[i0] = ir_start[i0] - itime
    return nhits

#%%
""" Function related to satellite read """

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
            dat0._get_IR0()
            dat0.close()
            # remove dat and make it a view of dat0, try to avoid errors on first attempt
            print('read GridSat for ',current_time)
            try: del dat
            except: pass
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
            # extend the lease while keeping the old dat
            dat.ti -= dtRange
        except FileNotFoundError:
            print('GridSat file not found ',current_time)
            # extend the lease while keeping the old dat
            dat.ti -= dtRange
        current_time -= dtRange
        
        yield dat

def check(dat,t):
        # check that the zone is not expired
    try:
        test = (t > dat.ti) and (t <= dat.tf)
    # Exception for the first usage when the time keys are not defined
    except KeyError:
        print('check KeyError')
        test = False
    except AttributeError:
        print('check AttributeError')
        test = False
    return test
#

if __name__ == '__main__':
    main()