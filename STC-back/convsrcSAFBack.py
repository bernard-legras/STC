# -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the backward trajectories from isentropiv or
isobaric levels.
This version is based on the SAFNWC product which is available at the same time and same 
resolution as the satellite data.
We use the reprocessed version which is produced on a truncated image.

Notice that the part time is supposed to start at day+1 0h where day is the day of the flight.

Created on Sun Oct  8 14:03:20 2017

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
from SAFNWCnc import SAFNWC_CTTH
import geosat

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
# step in the cloudtop procedure
cloudtop_step = timedelta(hours=12)
# low p cut in the STC traczilla runs
lowpcut = 3000
# highpcut in the STC traczilla runs
highpcut = 50000
# if True print a lot oj junk
verbose = False
debug = False

# idx_orgn was not set to 1 but to 0 in M55 and GLO runs
IDX_ORGN = 0

#%%
"""@@@@@@@@@@@@@@@@@@@@@@@@   MAIN   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

def main():
    global IDX_ORGN
    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-d","--day",type=int,choices=1+np.arange(31),help="day")
    parser.add_argument("-a","--advect",type=str,choices=["OPZ","EAD","EAZ","EID","EIZ"],help="source of advecting winds")
    parser.add_argument("-s","--suffix",type=str,help="suffix for special cases")
    parser.add_argument("-l","--level",type=int,help="PT level")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
    parser.add_argument("-o","--opaq",type=bool,help="keep only opaque clouds")
    parser.add_argument("-g","--good",type=bool,help="keep only good clouds") 
    parser.add_argument("-t","--step",type=int,help="step in hour between two part files")
    
    # to be updated
    if socket.gethostname() == 'graphium':
        pass
    elif 'ciclad' in socket.gethostname():
        #root_dir = '/home/legras/STC/STC-M55'
        main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
        traj_dir = '/data/legras/flexout/STC/BACK'
        out_dir = '/data/legras/STC'
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        #root_dir = '/home/legras/STC/STC-M55'
        main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
        traj_dir = '/data/legras/flexout/STC/BACK'
        out_dir = '/homedata/legras/STC'
    elif socket.gethostname() == 'grapelli':
        pass
    elif socket.gethostname() == 'gort':
        pass
    else:
         print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')
         exit()

    """ Parameters """
    # to do (perhaps) : some parameters might be parsed from command line
    # step and max output time
    step = 6
    hmax = 1824
    dstep = timedelta (hours=step)
    # age limit in days
    age_bound = 44.
    # time width of the parcel slice
    slice_width = timedelta(minutes=5)
    # dtRange
    dtRange={'MSG1':timedelta(minutes=15),'Hima':timedelta(minutes=20)}
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)
    # default values of parameters
    # start date of the backward run, corresponding to itime=0 
    year=2017
    # 8 +1 means we start on September 1 at 0h and cover the whole month of August 
    month=8+1
    # day=1 should not be changed
    day=1
    advect = 'EAZ'
    suffix =''
    quiet = False
    opaq = False
    good = False
    level = 380
    args = parser.parse_args()
    if args.year is not None:
        year=args.year
    if args.month is not None:
        month=args.month+1
    if args.advect is not None:
        advect=args.advect
    if args.suffix is not None:
        suffix='-'+args.suffix
    if args.level is not None:
        level=args.level
    if args.quiet is not None:
        if args.quiet=='y':
            quiet=True
        else:
            quiet=False
    if args.good is not None:
        good = args.good
    if args.opaq is not None:
        opaq = args.opaq
    if opaq:
        good = True
    if args.step is not None:
        step = args.step

    # Update the out_dir with the platform
    if opaq:
        out_dir = os.path.join(out_dir,'STC-BACK-OUT-SAF-OPAQ')
    elif good:
        out_dir = os.path.join(out_dir,'STC-BACK-OUT-SAF-GOOD')
    else:
        out_dir = os.path.join(out_dir,'STC-BACK-OUT-SAF-ALL')
    
    sdate = datetime(year,month,day)
    # fdate defined to make output under the name of the month where parcels are released 
    fdate= sdate - timedelta(days=1)
    
    # Number of parcels launched per time slot (1 per degree on a 170x50 grid)
    granule_size = 8500
    # number of granules in 6 hours
    granule_step = 4*6
    # size of granules launched during 6 hours
    granule_quanta = granule_size * granule_step
    # Modification of granule_size and granule_quanta
    if 'FULL' in advect:
        granule_size = 28800
        granule_step = 6
        granule_quanta = granule_size * granule_step
        
    # Manage the file that receives the print output
    if quiet:
        # Output file
        print_file = os.path.join(out_dir,'out','BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.out')
        fsock = open(print_file,'w') 
        sys.stdout=fsock
    
    print('year',year,'month',month,'day',day)
    print('advect',advect)
    print('suffix',suffix)

    # Directory of the backward trajectories
    ftraj = os.path.join(traj_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix)

    # Output file
    out_file = os.path.join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.hdf5z')
    #out_file1 = os.path.join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.hdf5b')
    #out_file2 = os.path.join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.pkl')

    # Directories for the satellite cloud top files
    satdir ={'MSG1':os.path.join(main_sat_dir,'msg1','S_NWC'),\
             'Hima':os.path.join(main_sat_dir,'himawari','S_NWC')}

    """ Initialization of the calculation """
    # Initialize the grid
    gg = geosat.GeoGrid('FullAMA_SAFBox')
    # Initialize the dictionary of the parcel dictionaries
    partStep={}
    satmap = pixmap(gg)

    # Build the satellite field generator
    get_sat = {'MSG1': read_sat(sdate,'MSG1',dtRange['MSG1'],satdir['MSG1']),\
               'Hima': read_sat(sdate,'Hima',dtRange['Hima'],satdir['Hima'])}

    # Open the part_000 file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'part_000'),quiet=True)
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
        # Link the names as views
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
        After the launch of the earliest parcel along the flight track, there
        should not be any member in new.
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
        print('live a, b ',live_a.sum(),live_p.sum())
        del kept_a
        del kept_p

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
            
            while satmap.check('MSG1',datpart['time']) is False:
                # if not get next satellite image 
                datsat1 = next(get_sat['MSG1'])
                # Check that the image is available
                if datsat1 is not None:
                    pm1 = geosat.SatGrid(datsat1,gg)
                    pm1._sat_togrid('CTTH_PRESS')
                    #print('pm1 diag',len(datsat1.var['CTTH_PRESS'][:].compressed()),
                    #                 len(pm1.var['CTTH_PRESS'][:].compressed()))
                    pm1._sat_togrid('ctth_quality')
                    pm1._sat_togrid('ctth_conditions')
                    pm1._sat_togrid('ctth_status_flag')
                    pm1.attr = datsat1.attr.copy()
                    satmap.fill('MSG1',pm1,good,opaq)
                    del pm1
                    del datsat1
                else:
                    # if the image is missing, extend the lease
                    try:
                        satmap.extend('MSG1')
                    except:
                        # This handle the unlikely case where the first image is missing
                        continue
            while satmap.check('Hima',datpart['time']) is False:
                # if not get next satellite image 
                datsath = next(get_sat['Hima'])
                # Check that the image is available
                if datsath is not None:
                    pmh = geosat.SatGrid(datsath,gg)
                    pmh._sat_togrid('CTTH_PRESS')
                    #print('pmh diag',len(datsath.var['CTTH_PRESS'][:].compressed()),
                    #                 len(pmh.var['CTTH_PRESS'][:].compressed()))
                    pmh._sat_togrid('ctth_quality')
                    pmh._sat_togrid('ctth_conditions')
                    pmh._sat_togrid('ctth_status_flag')
                    pmh.attr = datsath.attr
                    satmap.fill('Hima',pmh,good,opaq)
                    del datsath
                    del pmh
                else:
                    # if the image is missing, extend the lease
                    try:
                        satmap.extend('Hima')
                    except:
                        # This handle the unlikely case where the first image is missing
                        continue
          
            """ PROCESS THE COMPARISON OF PARCEL PRESSURES TO CLOUDS """
            if len(datpart['x'])>0:
                nhits += convbirth(datpart['itime'],
                    datpart['x'],datpart['y'],datpart['p'],datpart['t'],datpart['idx_back'],\
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
            n_nohit = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT) == 0).sum()
        except:
            nlive = 0
            n_nohit =0
        print('end hour ',hour,'  numact', partpost['nact'], ' nexits',nexits,' nhits',nhits, ' nlive',nlive,' nohit',n_nohit)
        # check that nlive + nhits + nexits = numpart, should be true after the first day
        if part0['numpart'] != nexits + nhits + nlive + ndborne:
            print('@@@ ACHTUNG numpart not equal to sum ',part0['numpart'],nexits+nhits+nlive+ndborne)           
            
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
    #pickle.dump(prod0,gzip.open(out_file,'wb'))
    #try:
    #    dd.io.save(out_file1,prod0,compression='blosc')
    #except:
    #    print('error with dd blosc')
    try:
        dd.io.save(out_file,prod0,compression='zlib')
    except:
        print('error with dd zlib')
    # close the print file

    """ End of the procedure and storage of the result """
    #output file
    #dd.io.save(out_file2,prod0,compression='zlib')
    #pickle.dump(prod0,gzip.open(out_file,'wb'))
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
""" Function doing the comparison between parcels and clouds and setting the result field """

@jit(nopython=True)
def convbirth(itime, x,y,p,t,idx_back, flag,xc,yc,pc,tc,age, ptop, ir_start, x0,y0,stepx,stepy,binx,biny):
    nhits = 0
    for i in range(len(x)):
        idx = min(int(np.floor((x[i]-x0)/stepx)),binx-1)
        idy = min(int(np.floor((y[i]-y0)/stepy)),biny-1)
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

def read_sat(t0,sat,dtRange,satdir):
    """ Generator reading the satellite data.
    The loop is infinite; sat data are called when required until the end of
    the parcel loop. """
    # get dt from satmap
    dt = dtRange
    # initial time
    current_time = t0
    namesat={'MSG1':'msg1','Hima':'himawari'}
    while True:
        fname = os.path.join(satdir,current_time.strftime('%Y/%Y_%m_%d'))
        if sat=='MSG1':
            fname = os.path.join(fname,current_time.strftime('S_NWC_CTTH_MSG1_FULLAMA-VISIR_%Y%m%dT%H%M00Z.nc'))
        elif sat=='Hima':
            fname = os.path.join(fname,current_time.strftime('S_NWC_CTTH_HIMAWARI08_FULLAMA-NR_%Y%m%dT%H%M00Z.nc'))
        else:
            print('sat should be MSG1 or Hima')
            return
        try:
            dat = SAFNWC_CTTH(current_time,namesat[sat],BBname='SAFBox')
            dat._CTTH_PRESS()
            # This pressure i left in hPa to allow masked with the fill_value in sat_togrid
            # The conversion to Pa is made in fill
            dat.attr['dtRange'] = dt
            dat.attr['lease_time'] = current_time - dtRange
            dat.attr['date'] = current_time
            dat._get_var('ctth_status_flag')
            dat._get_var('ctth_conditions')
            dat._get_var('ctth_quality')
            dat.close()
        except FileNotFoundError:
            print('SAF file not found ',current_time,namesat[sat])
            dat = None
        current_time -= dtRange
        yield dat

#%%
""" Describe the pixel map that contains the 5' slice of cloudtop data used in
the comparison of parcel location """

class pixmap(geosat.GridField):

    def __init__(self,gg):
        
        geosat.GridField.__init__(self,gg)
        
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
        self.zone['MSG1']['dtRange'] = timedelta(minutes=15)
        self.zone['Hima']['dtRange'] = timedelta(minutes=20)
        # define the slice
        self.ptop = np.empty(shape=self.geogrid.shapeyx,dtype=np.float)
        self.ptop.fill(p0)
        self.num  = np.zeros(shape=self.geogrid.shapeyx,dtype=np.int)
        self.range = self.geogrid.box_range
        self.binx = self.geogrid.box_binx
        self.biny = self.geogrid.box_biny
        self.stepx = (self.range[0,1]-self.range[0,0])/self.binx
        self.stepy = (self.range[1,1]-self.range[1,0])/self.biny
        print('steps',self.stepx,self.stepy)
    #def set_mask(self):
         # define the regional mask of the pixmap
    #     pass

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
            #print('check', zone,self.zone[zone]['ti'],self.zone[zone]['tf'])
        # Exception for the first usage when the time keys are not defined
        except KeyError:
            print('check KeyError')
            test = False
        return test
    
    def extend(self,zone):
        self.zone[zone]['ti'] -= self.zone[zone]['dtRange']

    def fill(self,zone,dat,good,opaq):
        """ Function filling the zone with new data from the satellite dictionary.
        """
        # Erase the zone
        self.erase(zone)
        # Mask outside the new data outside the zone
        if zone == 'MSG1':
            dat.var['CTTH_PRESS'][:,self.zone['Hima']['xi']:] = np.ma.masked
        elif zone == 'Hima':
            dat.var['CTTH_PRESS'][:,:self.zone['MSG1']['binx']] = np.ma.masked
        nbValidBeforeSel = len(dat.var['CTTH_PRESS'].compressed())
        # Filter according to quality keeping only good retrievals with all input fields
        if good:
            sel = ((dat.var['ctth_quality']&0x38)==0x8) & ((dat.var['ctth_conditions']&0xFF00)==0x5500)
            dat.var['CTTH_PRESS'][~sel] = np.ma.masked
        # Filter the non opaque clouds if required
        if opaq:
            sel = (dat.var['ctth_status_flag']&4)==4
            dat.var['CTTH_PRESS'][~sel] = np.ma.masked
        # test : count the number of valid pixels
        #nbValidAfterSel = len(dat.var['CTTH_PRESS'].compressed())
        #print('valid pixels before & after selection',zone,nbValidBeforeSel,nbValidAfterSel)
        # Inject the non masked new data in the pixmap
        # Conversion to Pa is done here
        self.ptop[~dat.var['CTTH_PRESS'].mask] = 100*dat.var['CTTH_PRESS'][~dat.var['CTTH_PRESS'].mask]
        # set the new expiration date 
        self.zone[zone]['tf'] = dat.attr['date']
        self.zone[zone]['ti'] = dat.attr['lease_time']
        
        if debug:
            sel = self.ptop < p0
            nbact = 100*sel.sum()/(self.geogrid.box_binx*self.geogrid.box_biny)
            if verbose: print('fill ',zone,' #selec ',len(sel),' % {:4.2f}'.format(nbact),\
                ' meanP {:6.0f} minP {:5.0f}'.format(self.ptop[sel].mean(),self.ptop.min()))
        return

if __name__ == '__main__':
    main()
