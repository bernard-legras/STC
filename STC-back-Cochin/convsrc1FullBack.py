# -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the air sampled during the StratoClim
campaign.
This version is based on the cloudtop files calculated during summer 2017. It does not
exploit the full resolution of the original satellite data and does not account properly
of the weight of each satellite pixel.
This is partially compensated by using a fairly coarse mesh 0.1°x0.1°. When several
cloudtops fall within a mesh, we retain the highest one. This can be modified.

Created on Sun Oct  8 14:03:20 2017
Modified to process BACK runs 2 March 2018 
Modified to process the Cochin runs 9 December 2018

@author: Bernard Legras
"""

import socket
import numpy as np
from collections import defaultdict
from numba import jit
from datetime import datetime, timedelta
import os
from os.path import join
import deepdish as dd
import pickle
#import pickle, gzip
import sys
import argparse
import psutil

from io107 import readpart107, readidx107
p0 = 100000.
I_DEAD = 0x200000
I_HIT = 0x400000
I_OLD = 0x800000
I_CROSSED = 0x2000000
I_DBORNE =  0x1000000
I_STOP = I_HIT + I_DEAD

# if True print a lot oj junk
verbose = False
debug = False

# idx_orgn was not set to 1 but to 0 in M55 and GLO runs
IDX_ORGN = 0

#%%
"""@@@@@@@@@@@@@@@@@@@@@@@@   MAIN   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

def main():
    global IDX_ORGN
    global highpcut, lowpcut
    global cloudtop_step
    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-a","--advect",type=str,choices=["EID","EIZ"],help="source of advecting winds")
    #parser.add_argument("-l","--level",type=int,help="PT level")
    parser.add_argument("-t","--step",type=int,help="step in hours between two part files")
    parser.add_argument("-s","--suffix",type=str,help="suffix for special cases")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
    parser.add_argument("-d","--duration",type=int,help="duration of integration in hours")
    parser.add_argument("-gs","--granule_size",type=int,help="size of the granule")
    parser.add_argument("-gn","--granule_step",type=int,help="number of granules in a step")
    parser.add_argument("-ab","--age_bound",type=int,help="age_bound")
    parser.add_argument("-binx","--binx",type=int,help="number of grid points in longitude direction")
    parser.add_argument("-biny","--biny",type=int,help="number of grid points in latitude direction")
    parser.add_argument("-hmax","--hmax",type=int,help='maximum number of hours in traczilla simulation')
    parser.add_argument("-username","--username",type=str,help='username')
    parser.add_argument("-userout","--userout",type=str,help='userout')    
    
    """ Parsed parameters"""
    # Parsed parameters
    # Interval between two part files (in hours)
    step = 6
    # Largest time to be processed (in hours)
    hmax = 1440
    # Age limit in days
    age_bound = 30
    # start date of the backward run, corresponding to itime=0 
    year=2017
    # 8 +1 means we start on September 1 at 0h and cover the whole month of August 
    month=6+1 
    advect = 'EIZ'
    suffix ='_150_150hPa_500'
    quiet = False
    #level = 150
    username = "sbucci"
    userout = username.copy() 
    # Bound on the age of the parcel (in days)
    age_bound = 30
    # Number of parcels launched per time slot (grid size)
    binx = 320
    biny = 224
    # Number of granules in a step betwwen two part files
    granule_step = 6
    
    """ Non parsed parameters"""
    # Time width of the parcel slice
    slice_width = timedelta(minutes=5)
    # dtRange
    dtRange={'MSG1':timedelta(minutes=15),'Hima':timedelta(minutes=20)}
    # day=1 should not be changed
    day=1
    # low p cut applied in exiter
    lowpcut = 3000
    # high p cut applied in exiter
    highpcut = 50000
    # step in the cloudtop procedure
    cloudtop_step = timedelta(hours=12)
    
    args = parser.parse_args()
    if args.year is not None: year=args.year
    if args.month is not None: month=args.month+1
    if args.advect is not None: advect=args.advect
    if args.suffix is not None: suffix=args.suffix
    #if args.level is not None: level=args.level
    if args.quiet is not None:
        if args.quiet=='y': quiet=True
        else: quiet=False
    if args.step is not None: step = args.step
    if args.hmax is not None: hmax = args.hmax
    if args.binx is not None: binx = args.binx
    if args.biny is not None: biny = args.biny
    granule_size = binx*biny
    if args.granule_size is not None: granule_size = args.granule_size
    if args.duration is not None: hmax = args.duration
    if args.age_bound is not None: age_bound = args.age_bound
    if args.granule_step is not None: granule_step = args.granule_step
    if args.username is not None: username = args.username
    if args.userout is not None: userout = args.userout
        
    # Define main directories
    main_sat_dir = '/bdd/STRATOCLIM/flexpart_in'
    if 'ciclad' in socket.gethostname():
        traj_dir = join('/data/',username,'flexout','COCHIN','BACK')
        out_dir = join('/data',userout,'STC')
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):      
        traj_dir = join('/homedata/',username,'flexout','COCHIN','BACK')
        out_dir = join('/homedata',userout,'STC')
    else:
         print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')
         exit()
    # Update the out_dir with the platform
    out_dir = join(out_dir,'STC-BACK-OUT-Cochin')

    sdate = datetime(year,month,day)
    # fdate defined to make output under the name of the month where parcels are released 
    fdate= sdate - timedelta(days=1)
    
    # Number of slices between two outputs
    dstep = timedelta (hours=step)
    nb_slices = int(dstep/slice_width)
    
    # size of granules launched during a step
    granule_quanta = granule_size * granule_step
    
    # Manage the file that receives the print output
    if quiet:
        # Output file
        print_file = join(out_dir,'out','BACK-'+advect+fdate.strftime('-%b-%Y')+suffix+'.out')
        fsock = open(print_file,'w') 
        sys.stdout=fsock    

    print('year',year,'month',month,'day',day)
    print('advect',advect)
    print('suffix',suffix)
    
    # Directory of the backward trajectories
    ftraj = join(traj_dir,'BACK-'+advect+fdate.strftime('-%b-%Y')+suffix)

    # Output file
    out_file = join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y')+suffix+'.hdf5')
    out_file1 = join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y')+suffix+'.hdf5b')
    out_file2 = join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y')+suffix+'.pkl')

    # Directories for the satellite cloud top files
    satdir ={'MSG1':join(main_sat_dir,'StratoClim+1kmD_msg1-c'),\
             'Hima':join(main_sat_dir,'StratoClim+1kmD_himawari-d')}

    """ Initialization of the calculation """
    # Initialize the slice map to be used as a buffer for the cloudtops
    satmap = pixmap()
    satfill = {}
    datsat = {}
    # Initialize the dictionary of the parcel dictionaries
    partStep={}

    # Build the satellite field generator
    get_sat = {'MSG1': read_sat(sdate,dtRange['MSG1'],satdir['MSG1'],pre=True),\
               'Hima': read_sat(sdate,dtRange['Hima'],satdir['Hima'],pre=True)}

    # Open the part_000 file that contains the initial positions
    part0 = readidx107(join(ftraj,'part_000'),quiet=True)
    print('numpart',part0['numpart'])
    numpart = part0['numpart']
    numpart_s = granule_size
  
    # stamp_date not set in these runs
    # current_date actually shifted by one day / sdate
    current_date = sdate
    # check flag is clean
    print('check flag is clean ',((part0['flag']&I_HIT)!=0).sum(),((part0['flag']&I_DEAD)!=0).sum(),\
                                 ((part0['flag']&I_CROSSED)!=0).sum())
    # check idx_orgn
    if part0['idx_orgn'] != 0:
        print('MINCHIA, IDX_ORGN NOT 0 AS ASSUMED, CORRECTED WITH READ VALUE')
        print('VALUE ',part0['idx_orgn'])
        IDX_ORGN = part0['idx_orgn']
    idx1 = IDX_ORGN

    # Build a dictionary to host the results
    prod0 = defaultdict(dict)
    prod0['src']['x'] = np.full(part0['numpart'],fill_value=np.nan,dtype='float')
    prod0['src']['y'] = np.full(part0['numpart'],fill_value=np.nan,dtype='float')
    prod0['src']['p'] = np.full(part0['numpart'],fill_value=np.nan,dtype='float')
    prod0['src']['t'] = np.full(part0['numpart'],fill_value=np.nan,dtype='float')
    prod0['src']['age'] = np.full(part0['numpart'],fill_value=np.nan,dtype='int')
    prod0['flag_source'] = part0['flag']
    prod0['rvs'] = np.full(part0['numpart'],0.01,dtype='float')
    
    # truncate eventually to 32 bits at the output stage

    # read the part_000 file for the first granule
    partStep[0] = {}
    partStep[0]['x']=part0['x'][:granule_size]
    partStep[0]['y']=part0['y'][:granule_size]
    partStep[0]['t']=part0['t'][:granule_size]
    partStep[0]['p']=part0['p'][:granule_size]
    partStep[0]['t']=part0['t'][:granule_size]
    partStep[0]['idx_back']=part0['idx_back'][:granule_size]
    partStep[0]['ir_start']=part0['ir_start'][:granule_size]
    partStep[0]['itime'] = 0

    # number of hists and exits
    nhits = 0
    nexits = 0
    ndborne = 0
    nnew = granule_size
    nold = 0

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
        
        # Processing of water mixing ratio
        # Select non stopped parcels in partante
        selec = (prod0['flag_source'][partante['idx_back']-IDX_ORGN] & I_STOP) == 0
        idx = partante['idx_back'][selec]
        prod0['rvs'][idx-IDX_ORGN] = np.minimum(prod0['rvs'][idx-IDX_ORGN],\
                satratio(partante['p'][selec],partante['t'][selec]))
             
        """ Select the parcels that are common to the two steps
        ketp_a is a logical field with same length as partante
        kept_p is a logical field with same length as partpost
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
            # ACHTUNG ACHTUNG : this works because IDX_ORGN=1, FIX THAT
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
                part0['ir_start'], satmap.range)
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
        print('live a, p ',live_a.sum(),live_p.sum())
        del kept_a
        del kept_p
        
        """ Correction step that moves partante['x'] to avoid big jumps at the periodicity frontier on the x-axis """
        diffx = partpost['x'][live_p] - partante['x'][live_a]
        bb = np.zeros(len(diffx))
        bb[diffx>180] = 360
        bb[diffx<-180] = -360
        partante['x'][live_a] += bb
        del bb
        del diffx
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
            if verbose: print('part slice ',i, datpart['time'])
            # Make sure the present satellite slice is OK
            # The while should ensure that the run synchronizes
            # when it starts.
            while satmap.check('MSG1',datpart['time']) is False:
                # if not get next satellite slice
                try:
                    void = next(satfill['MSG1'])
                # read new satellite file if the slice generator is over
                # make a new slice generator and get first slice
                except:
                    datsat['MSG1'] = next(get_sat['MSG1'])
                    satfill['MSG1'] = satmap.fill('MSG1',datsat)
                    void = next(satfill['MSG1'])
                finally:
                    if verbose: print('check MSG1 ',satmap.check('MSG1',datpart['time']),'##',datpart['time'],
                          '##',satmap.zone['MSG1']['ti'],'##',satmap.zone['MSG1']['tf'])
            while satmap.check('Hima',datpart['time']) is False:
                try:
                    void = next(satfill['Hima'])
                except:
                    datsat['Hima'] = next(get_sat['Hima'])
                    satfill['Hima'] = satmap.fill('Hima',datsat)
                    void = next(satfill['Hima'])
                finally:
                    if verbose: print('check Hima ',satmap.check('Hima',datpart['time']),'##',datpart['time'],
                          '##',satmap.zone['Hima']['ti'],'##',satmap.zone['Hima']['tf'])
            
            """ Select the parcels located within the domain """
             # TODO TODO the values used here should be derived from parameters defined above
            indomain = np.all((datpart['x']>-10,datpart['x']<160,datpart['y']>0,datpart['y']<50),axis=0)
            
            """ PROCESS THE COMPARISON OF PARCEL PRESSURES TO CLOUDS """
            if indomain.sum() >0:
                nhits += convbirth(datpart['itime'],
                    datpart['x'][indomain],datpart['y'][indomain],datpart['p'][indomain],\
                    datpart['t'][indomain],datpart['idx_back'][indomain],\
                    prod0['flag_source'],prod0['src']['x'],prod0['src']['y'],\
                    prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],\
                    satmap.ptop, part0['ir_start'],\
                    satmap.range[0,0],satmap.range[1,0],satmap.stepx,satmap.stepy,satmap.binx,satmap.biny)

            sys.stdout.flush()
            
            """ End of of loop on slices """

        # Check the age limit (easier to do it here)
        print("Manage age limit",flush=True)
        age_sec = part0['ir_start'][partante['idx_back']-IDX_ORGN]-partante['itime']
        IIold_o = age_sec > (age_bound-0.25) * 86400
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
            #n_nohit = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT) == 0).sum()
        except:
            nlive = 0
            #n_nohit =0
        print('end hour ',hour,'  numact', partpost['nact'],' nnew',nnew, ' nexits',nexits,' nhits',nhits, ' nlive',nlive,\
              ' nold',nold,' ndborne',ndborne)
        # check that nold + nlive + nhits + nexits = nnew
        if nnew != nexits + nhits + nlive + nold:
            print('@@@ ACHTUNG nnew not equal to sum ',nnew,nexits+nhits+nlive+nold)

    """ End of the procedure and storage of the result """
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use before clean: {:4.2f} gb'.format(memoryUse))
    del partante
    del partpost
    del live_a
    del live_p
    del datsat
    del datpart
    del satfill
    prod0['rvs'] = prod0['rvs'].astype(np.float32)
    for var in ['age','p','t','x','y']:
        prod0['src'][var] = prod0['src'][var].astype(np.float32)
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use after clean: {:4.2f} gb'.format(memoryUse))
    #output file
    try:
        print('output with deepdish/zlib')
        dd.io.save(out_file,prod0,compression='zlib')
    except:
        print('error with dd zlib')
        try:
           print('try with deepdish/blosc')
           dd.io.save(out_file1,prod0,compression='blosc')
        except:
           print('error with blosc')
           try:
              print('last attempt with pickle no compression')
              pickle.dump(prod0,open(out_file2,'wb'))
           except:
              print('pickle failed, NO OUTPUT')
    # close the print file
    if quiet: fsock.close()

"""@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

#%%
def satratio(p,T):
    """ Calculate the mass saturation ratio from pressure (in Pa) and temperature 
    (in K). Output in ppmv 
    usual factor 0.622 multiplied per 29/18 """
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
           dat['idx_back'] = dat['x'] = dat['y'] = dat['p'] = dat['t'] = np.array([])
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
""" Function managing the exiting parcels.
    We only consider exiting through top or bottom"""

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
            if p[i] > highpcut - 150: excode = 1
            elif p[i] < lowpcut  + 15 : excode = 2
            else:                   excode = 7
            flag[i0] |= (excode << 13) + I_DEAD + I_CROSSED
    return nexits

#%%
""" Function doing the comparison between parcels and clouds and setting the result field 
    Parcels outside the domain are not accounted."""

@jit(nopython=True)
def convbirth(itime, x,y,p,t,idx_back, flag,xc,yc,pc,tc,age, ptop, ir_start, x0,y0,stepx,stepy,binx,biny):
    nhits = 0
    for i in range(len(x)):
        idx = int(np.floor((x[i]-x0)/stepx))
        idy = int(np.floor((y[i]-y0)/stepy))
        #if any([idx<0, idy<0, idx>binx-1, idy>biny-1]):
        if (idx<0) or (idy<0) or (idx>binx-1) or (idy>biny-1):
            continue
        if ptop[idy,idx] < p[i]:
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

def read_sat(t0,dtRange,satdir,pre=False):
    """ Generator reading the satellite data.
    The loop is infinite; sat data are called when required until the end of
    the parcel loop. """
    # get dt from satmap
    dt = dtRange
    # initial time
    current_time = t0
    while True:
        fname = os.path.join(satdir,current_time.strftime('%Y%m%d%H_TB230'))
        dat = readidx107(fname,quiet=True)
        """ Generate the sequence of time ranges.
        This procedure works with empty time slots """
        # get the list of discontinuities, note that it is turned to a list
        id = list(np.where(dat['ir_start'][1:]-dat['ir_start'][:-1])[0])
        #test print(dat['ir_start'][id+1],dat['ir_start'][0])
        # append and prepend last and pre-first positions
        id.append(len(dat['ir_start'])-1)
        id[:0] = [-1]
        #test  print(id)
        dat['dtRange'] = dt
        dat['nt'] = int(cloudtop_step/dtRange)
        #test print('nt ',dat['nt'])
        tf = current_time + cloudtop_step/2
        # Generate list of time intervals
        dat['time'] = [[tf-dt,tf]]
        if pre:          
           off = timedelta()
           while tf >= current_time - cloudtop_step/2:
               dat['time'].append([tf,tf+dt])
               if debug: print('times',[tf,tf+dt])
               tf -= dt
        else:   
           off = -dt
           while tf > current_time - cloudtop_step/2 :
               tf -= dt
               dat['time'].append([tf-dt,tf])
               
        dat['time'].reverse()
        #test print(len(dat['time']))
        dat['numRange'] = np.zeros(dat['nt'],dtype='int')
        dat['indexRange'] = np.empty(shape=(dat['nt'],2),dtype='int')
        dat['indexRange'].fill(-999)
        # index in the list of time segments
        nc = dat['nt']-1
        # process the list of crossing from the last one, backward in time
        while len(id) > 1:
            idc = id.pop()
            # find the corresponding segment, skipping empty ones
            while off+current_time+timedelta(seconds=int(dat['ir_start'][idc])) \
                        != dat['time'][nc][0]:
                nc -= 1
            dat['indexRange'][nc,:] = [id[-1]+1,idc+1]
            dat['numRange'][nc] = idc - id[-1]
            nc -= 1
        # check that all parcels are sorted
        print('check sorting ',np.sum(dat['numRange']),dat['numpart'])
        # iterate time
        current_time -= cloudtop_step
        yield dat

#%%
""" Describe the pixel map that contains the 5' slice of cloudtop data used in
the comparison of parcel location """

class pixmap(object):

    def __init__(self):
        self.range = np.array([[-10.,160.],[0.,50.]])
        self.binx=1700; self.biny=500
        self.xedges = np.arange(self.range[0,0],self.range[0,1]+0.001,\
                               (self.range[0,1]-self.range[0,0])/self.binx)
        self.yedges = np.arange(self.range[1,0],self.range[1,1]+0.001,\
                               (self.range[1,1]-self.range[1,0])/self.biny)
        self.xcent = 0.5*(self.xedges[1:] + self.xedges[:-1])
        self.ycent = 0.5*(self.yedges[1:] + self.yedges[:-1])
        self.stepx = 0.1
        self.stepy = 0.1
        self.zone = defaultdict(dict)
        self.zone['MSG1']['range'] = np.array([[-10.,90.],[0.,50.]])
        self.zone['Hima']['range'] = np.array([[90.,160.],[0.,50.]])
        self.zone['MSG1']['binx'] = 1000
        self.zone['Hima']['binx'] = 700
        self.zone['MSG1']['biny'] = 500
        self.zone['Hima']['biny'] = 500
        self.zone['MSG1']['xi'] = 0
        self.zone['Hima']['xi'] = 1000
        self.zone['MSG1']['yi'] = 0
        self.zone['Hima']['yi'] = 0
        self.zone['MSG1']['dtRange'] = timedelta(minutes=30)
        self.zone['Hima']['dtRange'] = timedelta(minutes=20)
        # define the slice
        self.ptop = np.empty(shape=(self.biny,self.binx),dtype=np.float)
        self.ptop.fill(p0)
        self.num  = np.zeros(shape=(self.biny,self.binx),dtype=np.int)

    def set_mask(self):
         # define the regional mask of the pixmap
         pass

    def erase(self,zone):
         # erase the data in the zone
         x1 = self.zone[zone]['xi']
         x2 = x1 + self.zone[zone]['binx']
         y1 = self.zone[zone]['yi']
         y2 = y1 + self.zone[zone]['biny']
         self.ptop[y1:y2,x1:x2].fill(p0)
         self.num[y1:y2,x1:x2].fill(0)

    def check(self,zone,t):
         # check that the zone is not expired
         try:
             test = (t > self.zone[zone]['ti']) and (t <= self.zone[zone]['tf'])
         # Exception for the first usage when the time keys are not defined
         except KeyError:
             test = False
         return test

    def fill(self,zone,dat):
        """ Generator filling the slice with new data from the satellite dictionary
        of cloudtop pixels, using the precalculated time ranges.
        The data are read from the end to the beginning to fit the backward analysis
        The indexRange is such that the first value point to the end of the last
        range, the second to the end of the last-1 range and so on
        """
        for i in reversed(range(dat[zone]['nt'])):
            self.erase(zone)
            self.zone[zone]['tf'] = dat[zone]['time'][i][1]
            self.zone[zone]['ti'] = dat[zone]['time'][i][0]
            if dat[zone]['numRange'][i] >0:
                selec = range(dat[zone]['indexRange'][i,0],dat[zone]['indexRange'][i,1])

                #idx = np.floor((dat[zone]['x'][selec] - self.range[0,0])/self.stepx).astype('int')
                #idy = np.floor((dat[zone]['y'][selec] - self.range[1,0])/self.stepy).astype('int')
                fillfast(self.ptop,self.num,dat[zone]['x'][selec], \
                         dat[zone]['y'][selec],dat[zone]['p'][selec], \
                         self.range[0,0],self.range[1,0],self.stepx,self.stepy,self.binx,self.biny)
                if debug:
                    sel = self.ptop < p0
                    nbact = 100*sel.sum()/(self.binx*self.biny)
                    if verbose: print('fill ',zone,' #selec ',len(selec),' % {:4.2f}'.format(nbact),\
                          ' meanP {:6.0f} minP {:5.0f}'.format(self.ptop[sel].mean(),self.ptop.min()),\
                          ' nmax ',self.num.max(),\
                          ' minx {:7.2f} maxx {:7.2f}'.format(dat[zone]['x'][selec].min(),dat[zone]['x'][selec].max()))

            yield i

""" Functions related to slicing the satellite images """

@jit(nopython=True)
def fillfast(pSlice,numSlice,x,y,p,x0,y0,stepx,stepy,binx,biny):
    for i in range(len(x)):
        idx = min(int(np.floor((x[i]-x0)/stepx)),binx-1)
        idy = min(int(np.floor((y[i]-y0)/stepy)),biny-1)
        numSlice[idy,idx] += 1
        pSlice[idy,idx] = min(p[i],pSlice[idy,idx])


if __name__ == '__main__':
    main()
