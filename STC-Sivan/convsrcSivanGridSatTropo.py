# -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the air sampled during the StratoClim
campaign.
This version is based on the GridSat product and compares the temperature of the parcel to the IR brightness temperature.
This updated version calculates the tropopause location every 3 hours and does not consider cloud intersections for
parcels which are above the cold point

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
import flammkuchen as fl
from scipy.interpolate import RegularGridInterpolator
#import SAFNWCnc
import geosat
from ECMWF_N import ECMWF
#import constants as cst
from io107 import readpart107, readidx107

p0 = 100000.
I_DEAD = 0x200000
I_HIT = 0x400000
I_OLD = 0x800000
I_CROSSED = 0x2000000
I_DBORNE =  0x1000000
I_STOP = I_HIT + I_DEAD
# ACHTUNG I_DBORNE has been set to 0x10000000 (one 0 more) in a number of earlier analysis
# prior to 18 March 2018

# misc parameters
# low p cut in the M55 traczilla runs
lowpcut = 3000
# highpcut in the M55 traczilla runs
highpcut = 50000

# if True print a lot of junk
verbose = False
debug = False

# idx_orgn was set 1
# this is checked
IDX_ORGN = 1

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
    parser.add_argument("-m1","--month1",type=int,choices=1+np.arange(12),help="month1")
    parser.add_argument("-m2","--month2",type=int,choices=1+np.arange(12),help="month2")
    parser.add_argument("-d1","--day1",type=int,choices=1+np.arange(31),help="day1")
    parser.add_argument("-d2","--day2",type=int,choices=1+np.arange(31),help="day2")
    parser.add_argument("-a","--advect",type=str,choices=["OPZ","EAD","EAZ","EID","EIZ"],help="source of advecting winds")
    parser.add_argument("-s","--suffix",type=str,help="suffix for special cases")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
    parser.add_argument("-l","--level",type=int,help="PT level")
    #parser.add_argument("-c","--clean0",type=bool,help="clean part_000")
    parser.add_argument("-t","--step",type=int,help="step in hour between two part files")
    #parser.add_argument("-ct","--cloud_type",type=str,choices=["meanhigh","veryhigh","silviahigh"],help="cloud type filter")
    #parser.add_argument("-k","--diffus",type=str,choices=['01','1','001'],help='diffusivity parameter')
    parser.add_argument("-v","--vshift",type=int,choices=[0,10],help='vertical shift')
    parser.add_argument("-hm","--hmax",type=int,help='maximum considered integration time (hour)')
    parser.add_argument("-b","--bak",type=str,choices=["y","n"],help="backup (y) or not (n)")
    parser.add_argument("-c","--cont",type=str,choices=["y","n"],help="continuation (y) or not (n)")
    parser.add_argument("-bs","--bakstep",type=int,help='backup step (hour)')

    # to be updated
    # Define main directories
    if 'ciclad' in socket.gethostname():
        #main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
            #SVC_Dir = '/bdd/CFMIP/SEL2'
        traj_dir = '/data/legras/flexout/STC/Sivan'
        out_dir = '/data/legras/STC'
        #out_dir = '/data/legras/STC'
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        traj_dir = '/data/legras/flexout/STC/Sivan'
        out_dir = '/homedata/legras/STC'
    else:
        print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')
        exit()

    """ Parameters """
    # Output step of the trajectories in hour
    step = 6
    hmax = 3360
    #hmax = 18
    dstep = timedelta (hours=step)
    # time width of the parcel slice
    # might be better passed as a parameter
    age_bound = 44.5
    # time width of the parcel slice
    slice_width = timedelta(minutes=5)
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)
    # Number of parcels launched per time slot (1 per degree on a 170x50 grid)
    granule_size = 8500
    # number of granules in step hours
    granule_step = 6
    # size of granules launched during step hours
    granule_quanta = granule_size * granule_step
    # default values of parameters
    # date of the flight
    year=2017
    month1 = 7
    month2 = 9
    day1=1
    day2=30
    advect = 'EAD'
    suffix =''
    quiet = False
    clean0 = True
    backup = False
    backup_step = 40 * step
    restart = False
    #cloud_type = 'silviahigh'
    #diffus = '01'
    vshift = 0
    level = 380
    super =''
    # Step of GridSat
    dtRange = timedelta(hours=3)
    args = parser.parse_args()
    if args.year is not None: year=args.year
    if args.month1 is not None: month1=args.month1
    if args.month2 is not None: month2=args.month2
    if args.advect is not None: advect=args.advect
    if args.day1 is not None: day1 = args.day1
    if args.day2 is not None: day2 = args.day2
    if args.level is not None: level=args.level
    if args.hmax is not None: hmax = args.hmax
    if args.suffix is not None: suffix='-'+args.suffix
    if args.quiet is not None:
        if args.quiet=='y': quiet=True
        else: quiet=False
    #if args.clean0 is not None: clean0 = args.clean0
    #if args.cloud_type is not None: cloud_type = args.cloud_type
    if args.step is not None: step = args.step
    #if args.diffus is not None: diffus = args.diffus
    #diffus = '-D' + diffus
    if args.vshift is not None:
        vshift = args.vshift
        if vshift > 0:
            super = '-super'+str(vshift)
    if args.bak is not None:
        if args.bak=='y': backup=True
        else: backup=False
    if args.bakstep is not None: backup_step = args.bakstep
    if args.cont is not None:
        if args.cont=='y': restart=True
        else: restart=False

    # Update the out_dir with the cloud type and the super paramater
    out_dir = os.path.join(out_dir,'STC-Sivan-OUT-GridSat-WMO'+super)
    try:
        os.mkdir(out_dir)
        os.mkdir(out_dir+'/out')
    except:
        print('out_dir directory already created')

    # Dates beginning and end
    date_beg = datetime(year=year, month=month1, day=day1, hour=0)
    date_end = datetime(year=year, month=month2, day=day2, hour=0)

    # Manage the file that receives the print output
    if quiet:
        # Output file
        #print_file = os.path.join(out_dir,'out','BACK-SVC-EAD-'+date_beg.strftime('%b-%Y-day%d-')+date_end.strftime('%d-D01')+'.out')
        print_file = os.path.join(out_dir,'out','Sivan-T2-'+advect+'-JAS-'+date_beg.strftime('%Y-')+str(level)+'K.out')
        fsock = open(print_file,'w')
        sys.stdout=fsock

    # initial time to read the sat files
    # should be after the latest launch date
    sdate = date_end + timedelta(days=1)
    print('year',year,'month1',month1,'day1',day1,'month2',month2,'day2',day2)
    print('advect',advect)
    print('suffix',suffix)

    # patch to fix the problem of the data hole on the evening of 2 August 2017
    # should not happen
    #if sdate == datetime(2017,8,2): sdate = datetime(2017,8,3,6)

    # Directory of the backward trajectories (input) and name of the output file
    ftraj = os.path.join(traj_dir,'Sivan-'+advect+'-JAS-'+date_beg.strftime('%Y-')+str(level)+'K')
    out_file2 = os.path.join(out_dir,'Sivan-T2-'+advect+'-JAS-'+date_beg.strftime('%Y-')+str(level)+'K.hdf5')
    # backup file
    if backup:
        bak_file_prod0 = os.path.join(out_dir,'Sivan-T2-'+advect+'-JAS-'+date_beg.strftime('%Y-')
                                      +str(level)+'K-backup-prod0.hdf5')
        bak_file_params = os.path.join(out_dir,'Sivan-T2-'+advect+'-JAS-'+date_beg.strftime('%Y-')
                                       +str(level)+'K-backup-params.hdf5')
    """ Initialization of the calculation """
    # Initialize the dictionary of the parcel dictionaries
    partStep={}
    # initialize a grid that will be used to before actually doing any read
    gg = geosat.GeoGrid('GridSat')

    # Build the satellite field generator
    get_sat = read_sat(sdate,dtRange,pre=True,vshift=vshift)
    # Build the ECMWF field generator (both available at 3 hour interval)
    get_ERA5 = read_ERA5(sdate,dtRange,pre=True)
    # Read the index file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'part_000'),quiet=False)
    print('numpart',part0['numpart'])
    numpart = part0['numpart']
    numpart_s = granule_size # because of parcels launched at time 0
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
    idx1 = IDX_ORGN

    # Build a dictionary to host the results
    prod0 = defaultdict(dict)
    prod0['src']['x'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['y'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['p'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['t'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['age'] = np.empty(part0['numpart'],dtype='int')
    prod0['flag_source'] = part0['flag']
    # Set rvs at arbitrary large value
    prod0['rvs'] = np.full(part0['numpart'],0.01,dtype='float')
    # truncate eventually to 32 bits at the output stage

    # read the part_000 file
    if not restart:
        partStep[0] = readpart107(0,ftraj,quiet=False)
        # cleaning is necessary for runs starting in the fake restart mode
        # otherwise all parcels are thought to exit at the first step
        if clean0:
            partStep[0]['idx_back']=[]
        # parcels with longitude east of zero degree are set to negative values
        partStep[0]['x'][partStep[0]['x']>180] -= 360

    # number of hists and exits
    nhits = 0
    nexits = 0
    ndborne = 0
    nnew = granule_size # because of parcels launched at time 0
    nold = 0
    offset = 0

    # initialize datsat to None to force first read
    datsat = None

    # used to get non borne parcels (presently unused)
    #new = np.empty(part0['numpart'],dtype='bool')
    #new.fill(False)

    print('Initialization completed')

    if restart:
        print('restart run')
        try:
            params = fl.load(bak_file_params)
            prod0 = fl.load(bak_file_prod0)
        except:
            print('cannot load backup')
            return -1
        #[offset,nhits,nexits,nold,ndborne,nnew,idx1,numpart_s,current_date] = params['params']
        [offset,nhits,nexits,nold,ndborne,nnew,current_date] = params['params']
        idx1 = 620501
        numpart_s = 620500
        partStep[offset] = readpart107(offset,ftraj,quiet=True)
        partStep[offset]['x'][partStep[offset]['x']>180] -= 360
        # Initialize sat and ERA5 yield one step ahead as a precaution
        get_sat = read_sat(current_date + dstep,dtRange,pre=True,vshift=vshift)
        get_ERA5 = read_ERA5(current_date + dstep,dtRange,pre=True)

    """ Main loop on the output time steps """
    for hour in range(step+offset,hmax+1,step):
        pid = os.getpid()
        py = psutil.Process(pid)
        memoryUse = py.memory_info()[0]/2**30
        print('memory use: {:4.2f} gb'.format(memoryUse))
        # Get rid of dictionary no longer used
        if hour >= 2*step:
            try: del partStep[hour-2*step]
            except : print('nothing to delete')

        # Read the new data
        partStep[hour] = readpart107(hour,ftraj,quiet=True)
        # Link the names as views
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

        # Processing of water mixing ratio
        # Select non stopped parcels in partante
        if len(partante['idx_back']) >0:
            selec = (prod0['flag_source'][partante['idx_back']-IDX_ORGN] & I_STOP) == 0
            idx = partante['idx_back'][selec]
            prod0['rvs'][idx-IDX_ORGN] = np.minimum(prod0['rvs'][idx-IDX_ORGN],\
                satratio(partante['p'][selec],partante['t'][selec]))

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
        nnew += len(partpost['x'])-kept_p.sum()

        """ PROCESSING OF DEADBORNE PARCELS
        Manage the parcels launched during the last 6-hour which have already
        exited and do not appear in posold or posact (borne dead parcels).
        These parcels are stored in the last part of posact, at most
        the last granule_quanta parcels. """
        if numpart_s < numpart :
            print("manage deadborne",flush=True)
            # First index of the current quanta """
            numpart_s += granule_quanta
            print("idx1",idx1," numpart_s",numpart_s)
            # Extract the last granule_size indexes from posacti, this is where the new parcels should be
            if hour==step:
                idx_act = partpost['idx_back']
            else:
                idx_act = partpost['idx_back'][-granule_quanta:]
            # Generate the list of indexes that should be found in this range
            # This should works for both values of IDX_orgn
            idx_theor = np.arange(idx1,numpart_s+IDX_ORGN)
            # Find the missing indexes in idx_act (make a single line after validation)
            kept_borne = np.in1d(idx_theor,idx_act,assume_unique=True)
            idx_deadborne = idx_theor[~kept_borne]
            # Process these parcels by assigning exit at initial location
            prod0['flag_source'][idx_deadborne-IDX_ORGN] = prod0['flag_source'][idx_deadborne-IDX_ORGN] | I_DEAD+I_DBORNE
            prod0['src']['x'][idx_deadborne-IDX_ORGN] = part0['x'][idx_deadborne-IDX_ORGN]
            prod0['src']['y'][idx_deadborne-IDX_ORGN] = part0['y'][idx_deadborne-IDX_ORGN]
            prod0['src']['p'][idx_deadborne-IDX_ORGN] = part0['p'][idx_deadborne-IDX_ORGN]
            prod0['src']['t'][idx_deadborne-IDX_ORGN] = part0['t'][idx_deadborne-IDX_ORGN]
            prod0['src']['age'][idx_deadborne-IDX_ORGN] = 0.
            print("number of deadborne ",len(idx_deadborne))
            ndborne += len(idx_deadborne)
            idx1 = numpart_s + IDX_ORGN

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

        """ Correction step that moves partante['x'] to avoid big jumps at the periodicity frontier on the x-axis
        To be activated in the FULL version"""
        # if kept_p.sum()>0:
        #     diffx = partpost['x'][live_p] - partante['x'][live_a]
        #     bb = np.zeros(len(diffx))
        #     bb[diffx>180] = 360
        #     bb[diffx<-180] = -360
        #     partante['x'][live_a] += bb
        #     del bb
        #     del diffx
        # del kept_p
        # del kept_a
        # DOES not work in the following WAY
        #partante['x'][live_a][diffx>180] += 360
        #partante['x'][live_a][diffx<-180] -= 360

        # Build generator for parcel locations of the 5' slices
        gsp = get_slice_part(partante,partpost,live_a,live_p,current_date,dstep,slice_width)
        if verbose: print('built parcel generator for ',current_date)

        """  MAIN LOOP ON THE PARCEL TIME SLICES  """

        for i in range(nb_slices):
            # get the slice for the particles
            datpart = next(gsp)
            # Check x within (-180,180), necessary when GridSat is used
            #if datpart['x'] != []:
            if len(datpart['x'])>0:
                datpart['x'] = (datpart['x']+180)%360 - 180
            if verbose: print('part slice ',i, datpart['time'])
            # Check whether the present satellite image is valid
            # The while should ensure that the run synchronizes
            # when it starts.
            while check(datsat, datpart['time']) is False:
                # if not get next satellite image
                datsat = next(get_sat)
                datera = next(get_ERA5)

            """ PROCESS THE COMPARISON OF PARCEL PRESSURES TO CLOUDS """
            if len(datpart['x'])>0:
                # tropopause pressure at the location of the parcels
                #print('ptropo x ',datpart['x'].min(),datpart['x'].max())
                ptrop = datera.fP(np.transpose([datpart['y'],datpart['x']]))
                nhits += convbirth(datpart['itime'],
                    datpart['x'],datpart['y'],datpart['p'],datpart['t'],datpart['idx_back'],\
                    prod0['flag_source'],ptrop,prod0['src']['x'],prod0['src']['y'],\
                    prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],\
                    datsat.var['IR0'], part0['ir_start'],\
                    datsat.geogrid.box_range[0,0],datsat.geogrid.box_range[1,0],datsat.geogrid.stepx,\
                    datsat.geogrid.stepy,datsat.geogrid.box_binx,datsat.geogrid.box_biny)

            sys.stdout.flush()

        """ End of of loop on slices """

        """ INSERT HERE CODE FOR PARCELS ENDING BY AGE
        Important: if the age limit is set in the run, the age_bound needs to
        be decremented by one output step. Otherwise, the parcel disappeared
        and is wrongly processed as crossed."""
        # Check the age limit (easier to do it here)
        if len(partante['idx_back']) >0:
            print("Manage age limit",flush=True)
            age_sec = part0['ir_start'][partante['idx_back']-IDX_ORGN]-partante['itime']
            IIold_o = age_sec > (age_bound-(step/24)) * 86400
            IIold_o = IIold_o & ((prod0['flag_source'][partante['idx_back']-IDX_ORGN] & I_STOP)==0)
            idx_IIold = partante['idx_back'][IIold_o]
            j_IIold_o = np.where(IIold_o)
            prod0['flag_source'][idx_IIold-IDX_ORGN] = prod0['flag_source'][idx_IIold-IDX_ORGN] | I_DEAD+I_OLD
            prod0['src']['x'][idx_IIold-IDX_ORGN] = partante['x'][j_IIold_o]
            prod0['src']['y'][idx_IIold-IDX_ORGN] = partante['y'][j_IIold_o]
            prod0['src']['p'][idx_IIold-IDX_ORGN] = partante['p'][j_IIold_o]
            prod0['src']['t'][idx_IIold-IDX_ORGN] = partante['t'][j_IIold_o]
            prod0['src']['age'][idx_IIold-IDX_ORGN] = ((part0['ir_start'][idx_IIold-IDX_ORGN]- partante['itime'])/86400)
            print("number of IIold ",len(idx_IIold))
            nold += len(idx_IIold)

        # find parcels still alive       if kept_p.sum()==0:
        try:
            nlive = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_DEAD) == 0).sum()
            n_nohit = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT) == 0).sum()
        except:
            nlive = 0
            n_nohit =0
        print('end hour ',hour,'  numact', partpost['nact'], ' nexits',nexits,' nhits',nhits, ' nlive',nlive,' nohit',n_nohit,' nold',nold)
        # check that nlive + nhits + nexits = numpart, should be true after the first day
        if part0['numpart'] != nexits + nhits + nlive + ndborne:
            print('@@@ ACHTUNG numpart not equal to sum ',part0['numpart'],nexits+nhits+nlive+ndborne)

        if backup & (hour%backup_step == 0):
            fl.save(bak_file_prod0,prod0)
            fl.save(bak_file_params,{'params':[hour,nhits,nexits,nold,ndborne,nnew,idx1,numpart_s,current_date]})

    """ End of the procedure and storage of the result """
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use before clean: {:4.2f} gb'.format(memoryUse))
    del partante
    del partpost
    del live_a
    del live_p
    del datpart
    prod0['rvs'] = prod0['rvs'].astype(np.float32)
    for var in ['age','p','t','x','y']:
        prod0['src'][var] = prod0['src'][var].astype(np.float32)
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use after clean: {:4.2f} gb'.format(memoryUse))
    #output file
    fl.save(out_file2,prod0,compression='zlib')
    #pickle.dump(prod0,gzip.open(out_file,'wb'))
    # close the print file
    if quiet: fsock.close()

"""@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

#%%
def satratio(p,T):
    """ Calculate the volume saturation ratio from pressure (in Pa) and temperature
    (in K). Output in vmr.
    usual factor 0.622 must be applied to get mass mixing ratio
    estar is in hPa """
    estar = 1.0008*np.exp(23.33086-(6111.72784/T)+0.15215*np.log(T))
    satr = 1.0021 * estar/(0.01*p-estar)
    return satr
vsatratio = np.vectorize(satratio)

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
def convbirth(itime, x,y,p,t,idx_back, flag,ptrop,xc,yc,pc,tc,age, BT, ir_start, x0,y0,stepx,stepy,binx,biny):
    nhits = 0
    for i in range(len(x)):
        # do not consider parcels as long as they are above the cold point
        if ptrop[i] > p[i]: continue
        idx = min(int(np.floor((x[i]-x0)/stepx)),binx-1)
        idy = min(int(np.floor((y[i]-y0)/stepy)),biny-1)
        if BT[idy,idx] < t[i]:
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
            dat0.var['IR0'][dat0.var['IR0']<0] = 9999
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
            dat._get_T()
            dat._mkp()
            dats = dat.shift2west(-179)
            # extraction in a domain that encompasses FullAMA
            datr0 = dats.extract(latRange=[-5,55],lonRange=[-15,165],varss=['P','T'])
            del dat, dats
            datr0._WMO()
            datr0.fP = RegularGridInterpolator((np.arange(-5,56),np.arange(-15,166)),datr0.d2d['pwmo'],bounds_error=True)
            print('get ERA5 tropopause for ',current_time)
            try: del datr
            except: pass # intended for the first pass when datr undefined
            datr = datr0 # datr as a view of datr0
            # We reproduce here the same sequence as for the satellite file
            # although perhapsnot appropriate (the ERA5 data are instantaneous)
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
        print('check AttributeError')
        test = False
    return test
#

if __name__ == '__main__':
    main()
