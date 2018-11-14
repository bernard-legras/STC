    # -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the air sampled during the StratoClim
campaign.
This version does not use satellite data but the detrainment rates provided by the ERA5.
It is derived from convsrc1 with heavy modifications.
TO DO: another version (convsrc2cc) that uses the cloud cover instead of the detrainement rates.
As this postprocessing uses ERA5 data, it is mostly consistent with runs made using the ERA5 winds
and heating rates.
However we use it here with ERA-interim trajectories in the FULL domain.
This requires a number of modification to restrain the analysis to part of the trajectories which are in 
the monsoon domain.

Created on Sat 3 February 2018

Modified 18 March 2018 to be adapted to BACK runs
Modified 22 March 2018 to be adapted to BACK FULL runs

@author: Bernard Legras
"""

import socket
import numpy as np
import math
from collections import defaultdict
from numba import jit, int64
from datetime import datetime, timedelta
import os
import pickle
import deepdish as dd
import sys
import argparse
import psutil
from sys import exit
from scipy.interpolate import RegularGridInterpolator
from ECMWF_N import ECMWF
from mki2d import tohyb

from io107 import readpart107, readidx107
p0 = 100000.
I_DEAD = 0x200000
I_HIT = 0x400000
I_OLD = 0x800000
I_CROSSED = 0x2000000
I_DBORNE =  0x1000000
I_STOP = I_HIT + I_DEAD

# misc parameters
# step in the ERA5 data
ERA5_step = timedelta(hours=1)
# low p cut in the M55 traczilla runs
lowpcut = 3000
# highpcut in the M55 traczilla runs
highpcut = 50000

# if True print a lot oj junk
verbose = False
debug = True

# idx_orgn was not set to 1 but to 0 in M55 and GLO runs
IDX_ORGN = 0

#%%
"""@@@@@@@@@@@@@@@@@@@@@@@@   MAIN   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

def main():
    global IDX_ORGN
    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-a","--advect",type=str,choices=["EID-FULL","EIZ-FULL"],help="source of advecting winds")
    parser.add_argument("-l","--level",type=int,help="PT level")
    parser.add_argument("-s","--suffix",type=str,help="suffix for special cases")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")

    # to be updated
    if socket.gethostname() == 'graphium':
        pass
    elif 'ciclad' in socket.gethostname():
        traj_dir = '/data/legras/flexout/STC/BACK'
        out_dir = '/data/legras/STC'
    elif socket.gethostname() == 'grapelli':
        pass
    elif socket.gethostname() == 'coltrane':
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
    slice_width = timedelta(hours=1)
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)
    # defines here the offset for the detrainment (100h)
    detr_offset = 1/(100*3600.)
    # defines the domain
    domain = np.array([[-10.,160.],[0.,50.]])
    
    # default values of parameters
    # start date of the backward run, corresponding to itime=0 
    year=2017
    # 8 +1 means we start on September 1 at 0h and cover the whole month of August
    month=8+1
    # Should not be changed
    day=1
    advect = 'EAD'
    suffix =''
    quiet = False
    level = 380
    args = parser.parse_args()
    if args.year is not None:
        year=args.year
    if args.month is not None:
        month=args.month+1
    if args.advect is not None:
        advect=args.advect
    if args.level is not None:
        level=args.level
    if args.suffix is not None:
        suffix='-'+args.suffix
    if args.quiet is not None:
        if args.quiet=='y':
            quiet=True
        else:
            quiet=False

    # Update the out_dir with the platform
    out_dir = os.path.join(out_dir,'STC-BACK-DETR-OUT')
    sdate = datetime(year,month,day)
    # fdate defined to make output under the name of the month where parcels are released 
    fdate = sdate - timedelta(days=1)
   
    """ Define granule_size and granule_quanta
    granule_size:  Number of parcels launched per time slot (1 per degree on a 170x50 grid)
    granule_step: number of granules in 6 hours
    granula_quanta:  size of granules launched during 6 hours """
    if 'FULL' in advect:
        granule_size = 28800
        granule_step = 6
        granule_quanta = granule_size * granule_step
    else:
        granule_size = 8500
        granule_step = 6*4
        granule_quanta = granule_size * granule_step

    # Manage the file that receives the print output
    if quiet:
        # Output file
        print_file = os.path.join(out_dir,'out','BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.out')
        fsock = open(print_file,'w')
        sys.stdout=fsock

    # initial time to read the sat files
    # should be after the end of the flight
    # and a 12h or 0h boundary
    print('year',year,'month',month,'day',day)
    print('advect',advect)
    print('suffix',suffix)

    # Directory of the backward trajectories
    ftraj = os.path.join(traj_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix)    

    # Output file
    out_file = os.path.join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.hdf5b')    
    out_file1 = os.path.join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.hdf5z')    
    out_file2 = os.path.join(out_dir,'BACK-'+advect+fdate.strftime('-%b-%Y-')+str(level)+'K'+suffix+'.pkl')    

    """ Initialization of the calculation """
    # Initialize the dictionary of the parcel dictionaries
    partStep={}   

    # Open the part_000 file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'part_000'),quiet=False)
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
    # Locations of the crossing and detrainement
    nsrc = 6
    prod0['src']['x'] = np.full(shape=(nsrc,part0['numpart']),fill_value=np.nan,dtype='float')
    prod0['src']['y'] = np.full(shape=(nsrc,part0['numpart']),fill_value=np.nan,dtype='float')
    prod0['src']['p'] = np.full(shape=(nsrc,part0['numpart']),fill_value=np.nan,dtype='float')
    prod0['src']['t'] = np.full(shape=(nsrc,part0['numpart']),fill_value=np.nan,dtype='float')
    prod0['src']['age'] = np.full(shape=(nsrc,part0['numpart']),fill_value=np.nan,dtype='float')
    # Flag is copied from index
    prod0['flag_source'] = part0['flag']
    # Make a source array to accumulate the chi 
    # Dimension is that of the ERA5 field (201,681)
    prod0['source'] = np.zeros(shape=(201,681),dtype='float')
    # truncate eventually to 32 bits at the output stage

    # Inintialize the erosion 
    prod0['chi'] = np.full(part0['numpart'],1.,dtype='float')
    prod0['passed'] = np.full(part0['numpart'],10,dtype='int')
   
    # Build the interpolator to the hybrid level
    fhyb, void = tohyb()
    #vfhyb = np.vectorize(fhyb)

    # Read the part_000 file for the first granule
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
    nhits = np.array([0,0,0,0,0,0])
    nexits = 0
    ndborne = 0
    nnew = granule_size
    nold = 0
    nradada = 0

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
        the last granule_quanta parcels. 
        PB: this does not process the first parcels launched at time 0 since initially
        numpart_s = granule_size"""
        if numpart_s < numpart :
            print("manage deadborne",flush=True)
            # First index of the current quanta """
            numpart_s += granule_quanta
            print("numpart_s ",numpart_s) 
            # Extract the last granule_size indexes from posact
            if hour==step:
                idx_act = partpost['idx_back']
            else:    
                idx_act = partpost['idx_back'][-granule_quanta:]
            # Generate the list of indexes that should be found in this range
            idx_theor = np.arange(idx1,numpart_s+IDX_ORGN)
            # Find the missing indexes in idx_act (make a single line after validation)
            kept_borne = np.in1d(idx_theor,idx_act,assume_unique=True)
            idx_deadborne = idx_theor[~kept_borne]
            # Process these parcels by assigning exit at initial location
            prod0['flag_source'][idx_deadborne-IDX_ORGN] = prod0['flag_source'][idx_deadborne-IDX_ORGN] | I_DEAD+I_DBORNE
            prod0['src']['x'][0,idx_deadborne-IDX_ORGN] = part0['x'][idx_deadborne-IDX_ORGN]
            prod0['src']['y'][0,idx_deadborne-IDX_ORGN] = part0['y'][idx_deadborne-IDX_ORGN]
            prod0['src']['p'][0,idx_deadborne-IDX_ORGN] = part0['p'][idx_deadborne-IDX_ORGN]
            prod0['src']['t'][0,idx_deadborne-IDX_ORGN] = part0['t'][idx_deadborne-IDX_ORGN]
            prod0['src']['age'][0,idx_deadborne-IDX_ORGN] = 0.
            print("number of deadborne ",len(idx_deadborne))
            ndborne += len(idx_deadborne)
            idx1 = numpart_s + IDX_ORGN

        """ PROCESSING OF CROSSED PARCELS """
        # last known location before crossing stored in the index 0 of src fields
        if len(kept_a)>0:
            exits = exiter(int((partante['itime']+partpost['itime'])/2), \
                partante['x'][~kept_a],partante['y'][~kept_a],partante['p'][~kept_a],\
                partante['t'][~kept_a],partante['idx_back'][~kept_a],\
                prod0['flag_source'],prod0['src']['x'],prod0['src']['y'],\
                prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],\
                part0['ir_start'], domain)
            nexits += exits
            #nhits[0] += exits
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

        # Build generator for live parcel locations of the 1h slices
        gsp = get_slice_part(partante,partpost,live_a,live_p,current_date,dstep,slice_width)
        if verbose: print('built parcel generator for ',current_date)

        """  MAIN LOOP ON THE 1H PARCEL TIME SLICES  """

        for i in range(nb_slices):
            # get the 1h slice for the particles
            datpart = next(gsp)         
            # skip if no particles
            if datpart['ti'] == None:
                continue
            print('current_date ',datpart['ti'])
            #@@ test
#            print('ti in main ',datpart['ti'])
#            print('pi in main ',np.min(datpart['pi']),np.max(datpart['pi']))
#            print('idx_back   ',np.min(datpart['idx_back']),np.max(datpart['idx_back']))
#            #@@ end test
            # as the ECMWF files are also available every hour
            datrean = read_ECMWF(datpart['ti'])
            """ 
             Calculate the -log surface pressure at parcel location at time ti  
             create a 2D linear interpolar from the surface pressure field 
             This interpolator is meant to generate the surface pressure which is in turn used to determine
             sigma and the hybrid level. 
             The interpolator is only defined in the area where ECMWF ERA5 data are available."""                            
            lsp = RegularGridInterpolator((datrean.attr['lats'],datrean.attr['lons']),\
                                          -np.log(datrean.var['SP']),method='linear')
            
            """ Select the pairs (i,f) entirely located within the domain """
            indomain = np.all((datpart['xi']>-10,datpart['xi']<160,datpart['yi']>0,datpart['yi']<50,
                           datpart['xf']>-10,datpart['xf']<160,datpart['yf']>0,datpart['yf']<50),axis=0)
            
            # perform the interpolation for the location of live parcels at time ti
            # ACHTUNG: this interpolation is only meaningful for parcels laying within
            # the domain of the ERA5 data
            lspi = lsp(np.transpose([datpart['yi'][indomain],datpart['xi'][indomain]]))
            #@@ test
#            print('surface pressure ',np.exp(-np.min(datpart['lspi'])),np.exp(-np.max(datpart['lspi'])))
#            print('particle pressure ',np.min(datpart['pi']),np.max(datpart['pi'])) 
            #@@ end test
            # get the closest hybrid level at time ti
            # define first -log sigma = -log(p) - -log(ps)
            lsig = - np.log(datpart['pi'][indomain]) - lspi
            #@@ test
#            print('sigma ',np.exp(-np.max(lsig)),np.exp(-np.min(lsig)))
            #@@ end test
            # get the hybrid level, the rank of the first retained level is substracted to have hyb starting from 0 
            hyb = np.floor(fhyb(np.transpose([lsig,lspi]))+0.5).astype(np.int64)-datrean.attr['levs'][0]
            #@@ test the extreme values of sigma end ps
            if np.min(lsig) < - np.log(0.95):
                print('large sigma detected ',np.exp(-np.min(lsig)))
            if np.max(lspi) > -np.log(45000):
                print('small ps detected ',np.exp(-np.max(lspi)))
                
            """ PROCESS THE PARCELS WHICH ARE TOO CLOSE TO GROUND
             These parcels are flagged as crossed and dead, their last location is stored in the
             index 0 of src fields.
             This test handles also the cases outside the interpolation domain as NaN produced by fhyb
             generates very large value of hyb. 
             The trajectories which are stopped here have exited the domain where winds are available to flexpart
             and therefore are wrong from this point. For this reason we label them from their last valid position.
             This last sentence is only valid in simulations using ERA5 as ERA-Interim contains data down to the ground.
            """
            if np.max(hyb)> 100 :
                selec = hyb>100
                nr = radada(datpart['itime'],
                     datpart['xf'][indomain][selec],datpart['yf'][indomain][selec],datpart['pf'][indomain][selec],
                     datpart['tempf'][indomain][selec],datpart['idx_back'][indomain][selec],
                     prod0['flag_source'],prod0['src']['x'],prod0['src']['y'],
                     prod0['src']['p'],prod0['src']['t'],prod0['src']['age'] ,prod0['source'],
                     part0['ir_start'])
                nradada += nr
                nhits[0] += nr
       
            """ PROCESS THE (ADJOINT) DETRAINMENT """
            n1 = detrainer(datpart['itime'], 
                datpart['xi'][indomain],datpart['yi'][indomain],datpart['pi'][indomain],datpart['tempi'][indomain],hyb,
                datpart['xf'][indomain],datpart['yf'][indomain], datrean.var['UDR'], datpart['idx_back'][indomain],\
                prod0['flag_source'],part0['ir_start'], prod0['chi'],prod0['passed'],\
                prod0['src']['x'],prod0['src']['y'],prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],prod0['source'],\
                datrean.attr['Lo1'],datrean.attr['La1'],datrean.attr['dlo'],datrean.attr['dla'],detr_offset)
            nhits += n1
            #@@ test
            # print('return from detrainer',nhits)
            #@@ end test 
            del datpart
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
        prod0['src']['x'][0,idx_IIold-IDX_ORGN] = partante['x'][j_IIold_o]
        prod0['src']['y'][0,idx_IIold-IDX_ORGN] = partante['y'][j_IIold_o]
        prod0['src']['p'][0,idx_IIold-IDX_ORGN] = partante['p'][j_IIold_o]
        prod0['src']['t'][0,idx_IIold-IDX_ORGN] = partante['t'][j_IIold_o]
        prod0['src']['age'][0,idx_IIold-IDX_ORGN] = ((part0['ir_start'][idx_IIold-IDX_ORGN]- partante['itime'])/86400)
        print("number of IIold ",len(idx_IIold)) 
        nold += len(idx_IIold)
 
        
        # find parcels still alive       if kept_p.sum()==0:
        try:
            # number of parcels still alive
            nlive = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_DEAD) == 0).sum()
            # number of parcels still alive and not hit
            nprist = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & (I_DEAD+I_HIT)) == 0).sum()
            # number of parcels which have hit and crossed
            nouthit = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT+I_CROSSED) == I_HIT+I_CROSSED).sum()
            # number of parcels which heve crossed without hit
            noutprist = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT+I_CROSSED) == I_CROSSED).sum()
            # number of parcels which have hit without crossing
            nhitpure = ((prod0['flag_source'][partpost['idx_back']-IDX_ORGN] & I_HIT+I_CROSSED) == I_HIT).sum()                  
        except:
            nlive = 0
            nprist =0
            nouthit = 0
            noutprist = 0
            nhitpure = 0
            nprist = part0['numpart']
            
        print('end hour ',hour,'  numact', partpost['nact'], ' nnew',nnew,' nexits',nexits,' nold',nold,' ndborne',ndborne)
        print('nhits',nhits)
        print('nlive', nlive,' nprist',nprist,' nouthit',nouthit,' noutprist',noutprist,' nhitpure',nhitpure)
        # check that nprist + nhits + nexits + nold = nnew
        #if partpost['nact'] != nprist + nouthit + noutprist + nhitpure + ndborne:
        #    print('@@@ ACHTUNG numact not equal to sum ',partpost['nact'],nprist + nouthit + noutprist + nhitpure + ndborne)
      

    """ End of the procedure and storage of the result """
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use before clean: {:4.2f} gb'.format(memoryUse))
    del partante
    del partpost
    del live_a
    del live_p
    # reduction of the size of prod0 by converting float64 into float32
    prod0['chi'] = prod0['chi'].astype(np.float32)
    prod0['passed'] = prod0['passed'].astype(np.int32)
    prod0['source'] = prod0['source'].astype(np.float32)
    for var in ['age','p','t','x','y']:
        prod0['src'][var] = prod0['src'][var].astype(np.float32)
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use after clean: {:4.2f} gb'.format(memoryUse))

    #output file
    try:
        dd.io.save(out_file1,prod0,compression='blosc')
    except:
        print('error with dd blosc')
    try:
        dd.io.save(out_file,prod0,compression='zlib')
    except:
        print('error with dd zlib')
    try:
        pickle.dump(prod0,open(out_file2,'wb'))
    except:
        print('error with pickle')

    # close the print file
    if quiet: fsock.close()

"""@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

#%%
""" Functions related to the parcel data """

def get_slice_part(part_a,part_p,live_a,live_p,current_date,dstep,slice_width):
    """ Generator to generate 1h slices along flight track. 
    Each slice contains the coordinates of the beginning and the end of the interval.
    tp <= ti < tf <= ta """
    nb_slices = int(dstep/slice_width+0.0001)
    ta = current_date + dstep
    tp = current_date
    tf = ta
    empty_live = (live_a.sum() == 0)
    itime = part_p['itime']
    dat = {}
    for i in range(nb_slices):
        ti = tf - slice_width
        itime -= slice_width.seconds
        # note that 0.5*(ti+tf) cannot be calculated as we are adding two dates
        coefai = (ti-tp)/dstep
        coefpi = (ta-ti)/dstep
        dat['ti'] = ti
        dat['tf'] = tf
        dat['itime'] = itime
        if empty_live:
           dat['idx_back'] = dat['xi'] = dat['yi'] = dat['pi'] = dat['ti'] = np.array([])
           dat['xf'] = dat['yf'] = dat['pf'] = dat['tf'] = []
           dat['ti'] = dat['tf'] = None
        else:
            #@@ test
#            print('i, coefs  ',i,coefai,coefpi)
#            print('pres a    ',np.min(part_a['p'][live_a]),np.max(part_a['p'][live_a]))
#            print('pres p    ',np.min(part_p['p'][live_p]),np.max(part_p['p'][live_p]))
            #@@ end test
            if i == 0:
                dat['idx_back'] = part_a['idx_back'][live_a]
                dat['xf'] = part_a['x'][live_a]
                dat['yf'] = part_a['y'][live_a]
                dat['pf'] = part_a['p'][live_a]
                dat['tempf'] = part_a['t'][live_a]
                #dat['pf'] = part_a['p'][live_a]
                #dat['tempf'] = part_a['t'][live_a]
                dat['xi'] = coefai*part_a['x'][live_a] + coefpi*part_p['x'][live_p]
                dat['yi'] = coefai*part_a['y'][live_a] + coefpi*part_p['y'][live_p]
                dat['pi'] = coefai*part_a['p'][live_a] + coefpi*part_p['p'][live_p]
                #@@ test
#                print('pi        ',np.min(dat['pi']),np.max(dat['pi']))         
                #@@ end test
                dat['tempi'] = coefai*part_a['t'][live_a] + coefpi*part_p['t'][live_p]
            elif i == nb_slices-1:
                dat['xf'] = dat['xi']
                dat['yf'] = dat['yi']
                dat['pf'] = dat['pi']
                dat['tempf'] = dat['tempi']
                #dat['pf'] = dat['pi']
                #dat['tempf'] = dat['ti']
                dat['xi'] = part_p['x'][live_p]
                dat['yi'] = part_p['y'][live_p]
                dat['pi'] = part_p['p'][live_p]
                dat['tempi'] = part_p['t'][live_p]
            else:
                dat['xf'] = dat['xi']
                dat['yf'] = dat['yi']
                dat['pf'] = dat['pi']
                dat['tempf'] = dat['tempi']
                #dat['pf'] = dat['pi']
                #dat['tempf'] = dat['ti']
                dat['xi'] = coefai*part_a['x'][live_a] + coefpi*part_p['x'][live_p]
                dat['yi'] = coefai*part_a['y'][live_a] + coefpi*part_p['y'][live_p]
                dat['pi'] = coefai*part_a['p'][live_a] + coefpi*part_p['p'][live_p]
                dat['tempi'] = coefai*part_a['t'][live_a] + coefpi*part_p['t'][live_p]
        tf -= slice_width
        yield dat

#%%
""" Function managing the exiting parcels.
Consider only exit through top or bottom """

@jit(nopython=True,cache=True)
def exiter(itime, x,y,p,t,idx_back, flag,xc,yc,pc,tc,age, ir_start, rr):
    nexits = 0
    for i in range(len(x)):
        i0 = idx_back[i]-IDX_ORGN
        if flag[i0] & I_DEAD == 0:
            nexits += 1
            xc[0,i0] = x[i]
            yc[0,i0] = y[i]
            tc[0,i0] = t[i]
            pc[0,i0] = p[i]
            age[0,i0] = ir_start[i0] - itime
            if   p[i] > highpcut - 150: excode = 1
            elif p[i] < lowpcut  + 15 : excode = 2
            else:                   excode = 7
            flag[i0] |= (excode << 13) + I_DEAD + I_CROSSED
    return nexits

@jit(nopython=True,cache=True)
def radada(itime, x,y,p,t,idx_back, flag,xc,yc,pc,tc,age, ir_start):
    nexits =  0
    for i in range(len(x)):
        i0 = idx_back[i]-IDX_ORGN
        if flag[i0] & I_DEAD == 0:
            nexits += 1
            xc[0,i0] = x[i]
            yc[0,i0] = y[i]
            tc[0,i0] = t[i]
            pc[0,i0] = p[i]
            age[0,i0] = ir_start[i0] - itime
            flag[i0] |= (8 << 13) + I_DEAD + I_CROSSED
    return nexits
#%%

""" Function finding the detrainment at the location of the parcel and doing the job 
Do not process segments with an extermity outside the ERA5 domain"""

@jit(nopython=True,cache=True)
def detrainer(itime, xi,yi,pi,ti,hyb,xf,yf, udr, idx_back,flag,ir_start,chi,passed,\
              xc,yc,pc,tc,age,source,\
              Lo1,La1,dlo,dla,detr_offset):
    nhits = [0,0,0,0,0,0]
    # loop on the kept parcels
    for i in range(len(xi)):
        i0 = idx_back[i]-IDX_ORGN
        #@@ test
#        if i0<0 or i0>=7601000:
#            print('i0',i0)
        #@@ end test
        # consider only the live parcel
        if flag[i0] & I_DEAD ==0:
            # find integer coordinates of closest location on the mesh
            # It is assumed no point outside the domain
            xig = int(math.floor((xi[i]-Lo1)/dlo+0.5))
            yig = int(math.floor((yi[i]-La1)/dlo+0.5))
            xfg = int(math.floor((xf[i]-Lo1)/dla+0.5))
            yfg = int(math.floor((yf[i]-La1)/dla+0.5))
            # This line should be useless as filtering is already done at higher level
            #if any(xig<0, xfg<0, xig>680, xfg >680, yig<0, yfg<0, yig>200, yfg>200):
            if (xig<0) or (xfg<0) or (xig>680) or (xfg >680) or (yig<0) or (yfg<0) or (yig>200) or (yfg>200):
                continue
            #@@ test
#            if xig<0 or xig>679:
#                print('xig ',xig,'i ',i,'xi ',xi[i])
#            if xfg<0 or xfg>679:
#                print('xfg ',xfg)
#            if yig<0 or yig>199:
#                print('yig ',yig)
#            if yfg<0 or yfg>199:
#                print('yfg ',yfg)
#            if hyb[i]<0 or hyb[i]>100:
#                print('hyb ',hyb[i])
            #@@ end test
            # find the meshes on the path 
            ll = line(xig,yig,xfg,yfg)
            # calculate mean detrainment on the path
            detr = 0.
            for j in range(len(ll)):
                #@@ test
#                 if ll[j][1]<0 or ll[j][1]>200:
#                     print('ll1 ',ll[j])
#                 if ll[j][1]<0 or ll[j][1]>680:
#                     print('ll0 ',ll[j])
                #@@ test    
                detr += udr[hyb[i],ll[j][1],ll[j][0]]
            detr = detr/len(ll)
            # erode the parcel
            if detr >= detr_offset:
                newchi = chi[i0] * math.exp(-3600*detr)
                xm = int(0.5*(xig+xfg))
                ym = int(0.5*(yig+yfg))
                source[ym,xm] += chi[i0] - newchi
                chi[i0] = newchi
                if passed[i0] >1:
                    if passed[i0] == 10:
                        if chi[i0] < 0.9:
                            xc[1,i0] = xi[i]
                            yc[1,i0] = yi[i]
                            pc[1,i0] = pi[i]
                            tc[1,i0] = ti[i]
                            age[1,i0] = ir_start[i0] - itime
                            flag[i0] |= I_HIT
                            passed[i0] = 9
                            nhits[1] += 1
                    if passed[i0] == 9:
                        if chi[i0] < 0.7:
                            xc[2,i0] = xi[i]
                            yc[2,i0] = yi[i]
                            pc[2,i0] = pi[i]
                            tc[2,i0] = ti[i]
                            age[2,i0] = ir_start[i0] - itime
                            passed[i0] = 7
                            nhits[2] += 1
                    if passed[i0] == 7:
                        if chi[i0] < 0.5:
                            xc[3,i0] = xi[i]
                            yc[3,i0] = yi[i]
                            pc[3,i0] = pi[i]
                            tc[3,i0] = ti[i]
                            age[3,i0] = ir_start[i0] - itime
                            passed[i0] = 5
                            nhits[3] += 1
                    if passed[i0] == 5:
                        if chi[i0] < 0.3:
                            xc[4,i0] = xi[i]
                            yc[4,i0] = yi[i]
                            pc[4,i0] = pi[i]
                            tc[4,i0] = ti[i]
                            age[4,i0] = ir_start[i0] - itime
                            passed[i0] = 3
                            nhits[4] += 1
                    if passed[i0] == 3:
                        if chi[i0] < 0.1:
                            xc[5,i0] = xi[i]
                            yc[5,i0] = yi[i]
                            pc[5,i0] = pi[i]
                            tc[5,i0] = ti[i]
                            age[5,i0] = ir_start[i0] - itime
                            passed[i0] = 1
                            nhits[5] += 1
    return nhits

#%%
""" Function related to ECMWF read """

def read_ECMWF(date):
    """ Script reading the ECMWF data.
    Not a generator as this is synchronized with the 1h part slice.
    The data are assumed valid over the 1h period that follows the timestamp.
    This is quite OK for UDR as this quantity is defined as a mean/accumuation over
    this one-hour period.
    Cloud-cover is from analysis, therefore as an instantaneous map, but varies less rapidly 
    than the UDR. """
    dat = ECMWF('STC',date)
    dat._get_var('T')
    dat.attr['dlo'] = (dat.attr['lons'][-1] - dat.attr['lons'][0]) / (dat.nlon-1)
    dat.attr['dla'] = (dat.attr['lats'][-1] - dat.attr['lats'][0]) / (dat.nlat-1) 
    if dat.attr['Lo1']>dat.attr['Lo2']:
        dat.attr['Lo1'] = dat.attr['Lo1']-360.
    #@@ test
#    print('LoLa ',dat.attr['Lo1'],dat.attr['La1'],dat.attr['dlo'],dat.attr['dla'],\
#          dat.attr['levs'][0],len(dat.attr['levs']))
    #@@ end test
    dat._get_var('UDR')
    #dat._get_var('CC')
    dat._mkp()
    dat._mkrho()
    dat.var['UDR'] /= dat.var['RHO']
    dat.close()
    return dat

#%% Function that implements the Bresenham algorithn to makes lines between points
        
@jit((int64,int64,int64,int64),nopython=True,cache=True)
def line(x0, y0, x1, y1):
    i = 0
    points = np.empty(shape=(40,2),dtype=int64)
    if (x0==x1) & (y0==y1):
        points[i,0] = x0
        points[i,1] = y0
    else:
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        dx2 = 2*dx
        dy2 = 2*dy
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1  
        x, y = x0, y0 
        if dx > dy:
            err = dx 
            while x != x1:
                points[i,0] = x
                points[i,1] = y
                i += 1
                err -= dy2
                if err < 0:
                    y += sy
                    err += dx2
                x += sx
        else:
            err = dy 
            while y != y1:
                points[i,0] = x
                points[i,1] = y
                i += 1
                err -= dx2
                if err < 0:
                    x += sx
                    err += dy2
                y += sy
        points[i,0] = x1
        points[i,1] = y1
    return points[:i+1,:]        

if __name__ == '__main__':
    main()
