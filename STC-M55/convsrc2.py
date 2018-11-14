    # -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the air sampled during the StratoClim
campaign.
This version does not use satellite data but the detrainment rates provided by the ERA5.
It is derived from convsrc1 with heavy modifications.
TO DO: another version (convsrc4) that uses the cloud cover instead of the detrainement rates.
As this postprocessing uses ERA5 data, it is mostly consistent with runs made using the ERA5 winds
and heating rates.

Created on Sat 3 February 2018

@author: Bernard Legras
"""

import socket
import numpy as np
import math
from collections import defaultdict
from numba import jit, int64
from datetime import datetime, timedelta
import os
import pickle, gzip
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
I_CROSSED = 0x2000000
I_DBORNE =  0x1000000
# ACHTUNG I_DBORNE has been set to 0x10000000 (one 0 more) in a number of earlier analysis 
# prior to 18 March 2018

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
    parser.add_argument("-d","--day",type=int,choices=1+np.arange(31),help="day")
    parser.add_argument("-a","--advect",type=str,choices=["OPZ","EAD","EAZ","EID","EIZ"],help="source of advecting winds")
    parser.add_argument("-p","--platform",type=str,choices=["M55","GLO","BAL"],help="measurement platform")
    parser.add_argument("-n","--launch_number",type=int,help="balloon launch number within a day")
    parser.add_argument("-s","--suffix",type=str,help="suffix for special cases")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
    parser.add_argument("-c","--clean0",type=bool,help="clean part_000")

    # to be updated
    if socket.gethostname() == 'graphium':
        pass
    elif 'ciclad' in socket.gethostname():
        root_dir = '/home/legras/STC/STC-M55'
        traj_dir = '/data/legras/flexout/STC/M55'
        out_dir = '/data/legras/STC'
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        root_dir = '/home/legras/STC/STC-M55'
        traj_dir = '/bdd/STRATOCLIM/flexout/M55'
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
    hmax = 732
    dstep = timedelta (hours=step)
    # time width of the parcel slice
    slice_width = timedelta(hours=1)
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)
    # defines here the offset for the detrainment (100h)
    detr_offset = 1/(100*3600.)
    # defines the domain
    domain = np.array([[-10.,160.],[0.,50.]])
    # Size of each packet of parcels
    segl = 1000  
    
    # default values of parameters
    # date of the flight
    year=2017
    month=7
    day=29
    platform = 'GLO'
    advect = 'EAZ'
    suffix =''
    launch_number=''
    quiet = False
    clean0 = False
    args = parser.parse_args()
    if args.year is not None:
        year=args.year
    if args.month is not None:
        month=args.month
    if args.day is not None:
        day=args.day
    if args.advect is not None:
        advect=args.advect
    if args.platform is not None:
        platform=args.platform
    if args.launch_number is not None:
        launch_number='-'+str(args.launch_number)
    if args.suffix is not None:
        suffix='-'+args.suffix
    if args.quiet is not None:
        if args.quiet=='y':
            quiet=True
        else:
            quiet=False
    if args.clean0 is not None:
        clean0 = args.clean0

    # Update the out_dir with the platform
    out_dir = os.path.join(out_dir,'STC-'+platform+'-DETR-OUT')

    fdate = datetime(year,month,day)

    # Manage the file that receives the print output
    if quiet:
        # Output file
        print_file = os.path.join(out_dir,'out',platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix+'.out')
        fsock = open(print_file,'w')
        sys.stdout=fsock

    # initial time to read the detrainment files
    # should be after the end of the flight
    # and a 12h or 0h boundary
    print('year',year,'month',month,'day',day)
    print('advect',advect)
    print('platform',platform)
    print('launch_number',launch_number)
    print('suffix',suffix)
    
    # Read the region mask
    mm = pickle.load(gzip.open(os.path.join(root_dir,'STCmask2-era5.pkl')))
    mask = mm['mask']

    # Directory of the backward trajectories
    ftraj = os.path.join(traj_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix)

    # Output file
    #out_file = os.path.join(out_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix+'.pkl')
    out_file2 = os.path.join(out_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix+'.hdf5')

    """ Initialization of the calculation """
    # Initialize the dictionary of the parcel dictionaries
    partStep={}   

    # Read the index file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'index_old'),quiet=False)
    print('numpart',part0['numpart'])
      
    # stamp_date not set in these runs
    # current_date actually shifted by one day / sdate
    current_date = fdate + timedelta(days=1)
    # check flag is clean
    print('check flag is clean ',((part0['flag']&I_HIT)!=0).sum(),((part0['flag']&I_DEAD)!=0).sum(),\
                                 ((part0['flag']&I_CROSSED)!=0).sum())

    # check idx_orgn
    if part0['idx_orgn'] != 0:
        print('MINCHIA, IDX_ORGN NOT 0 AS ASSUMED, CORRECTED WITH READ VALUE')
        print('VALUE ',part0['idx_orgn'])
        IDX_ORGN = part0['idx_orgn']

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
    # Dimension is that of the ERA5 field (201,681) at 0.25° resolution
    # Both latitudes and longitudes are growing
    prod0['source'] = np.zeros(shape=(201,681),dtype='float')
    # truncate eventually to 32 bits at the output stage
    # Source array that cumulates within regions as a function of time
    prod0['pl'] = np.zeros(shape=(len(mm['regcode']),int(part0['numpart']/segl)),dtype='float') 
    
     # Ininitalize the erosion 
    prod0['chi'] = np.full(part0['numpart'],1.,dtype='float')
    prod0['passed'] = np.full(part0['numpart'],10,dtype='int')
   
    # Build the interpolator to the hybrid level
    fhyb, void = tohyb()
    #vfhyb = np.vectorize(fhyb)

    # Read the part_000 file
    partStep[0] = readpart107(0,ftraj,quiet=True)
    # cleaning is necessary for runs starting in the fake restart mode
    # otherwise all parcels are thought to exit at the first step
    if clean0:
        partStep[0]['idx_back']=[]

    # number of hists and exits
    nhits = np.array([0,0,0,0,0,0])
    nexits = 0
    ndborne = 0
    nnew = 0
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

        """ IDENTIFY AND TAKE CARE OF DEADBORNE AS NON BORNE PARCELS """
        # Find and count the parcels in partpost which where not in partante, flag them as new
        if (hour <= 30) & (partpost['nact']>0):
            new[partpost['idx_back'][~kept_p]-IDX_ORGN] = True
            nnew += len(partpost['x'])-kept_p.sum()
        # When the release of parcels has ended, flag the ones that did not show up as dead
        # source their last location as the launching location
        if hour == 30:
            ndborne = np.sum(~new)
            prod0['flag_source'][~new] |= I_DBORNE + I_DEAD
            prod0['src']['x'][0,~new] = part0['x'][~new]
            prod0['src']['y'][0,~new] = part0['y'][~new]
            prod0['src']['p'][0,~new] = part0['p'][~new]
            prod0['src']['t'][0,~new] = part0['t'][~new]
            prod0['src']['age'][0,~new] = 0
            print('number of dead borne',ndborne,part0['numpart']-nnew)
            # get rid of new
            del new     

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
            nhits[0] += exits
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

        # Build generator for live parcel locations of the 1h slices
        gsp = get_slice_part(partante,partpost,live_a,live_p,current_date,dstep,slice_width)
        if verbose: print('built parcel generator for ',current_date)

        """  MAIN LOOP ON THE PARCEL TIME SLICES  """

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
            # calculate the -log surface pressure at parcel location at time ti  
            # create a 2D linear interpolar from the surface pressure field
            lsp = RegularGridInterpolator((datrean.attr['lats'],datrean.attr['lons']),\
                                          -np.log(datrean.var['SP']))
            # perform the interpolation for the location of live parcels at time ti
            datpart['lspi'] = lsp(np.transpose([datpart['yi'],datpart['xi']]))
            #@@ test
#            print('surface pressure ',np.exp(-np.min(datpart['lspi'])),np.exp(-np.max(datpart['lspi'])))
#            print('particle pressure ',np.min(datpart['pi']),np.max(datpart['pi'])) 
            #@@ end test
            # get the closest hybrid level at time ti
            # define first -log sigma = -log(p) - -log(ps)
            lsig = - np.log(datpart['pi']) - datpart['lspi']
            #@@ test
#            print('sigma ',np.exp(-np.max(lsig)),np.exp(-np.min(lsig)))
            #@@ end test
            # get the hybrid level, the rank of the first retained level is substracted to have hyb starting from 0
            # +1 because the levels are counted from 1, not 0 and +0.5 because we get the closest neighbour
            hyb = np.floor(fhyb(np.transpose([lsig,datpart['lspi']]))+1.5).astype(np.int64)-datrean.attr['levs'][0]
            #@@ test the extreme values of sigma end ps
            if np.min(lsig) < - np.log(0.95):
                print('large sigma detected ',np.exp(-np.min(lsig)))
            if np.max(datpart['lspi']) > -np.log(45000):
                print('small ps detected ',np.exp(-np.max(datpart['lspi'])))
                
            """ PROCESS THE PARCELS WHICH ARE TOO CLOSE TO GROUND
             These parcels are flagged as crossed and dead, their last location is stored in the
             index 0 of src fields.
             This test handles also the cases outside the interpolation domain as NaN produced by fhyb
             generates very large value of hyb. 
             The trajectories which are stopped here have exited the domain where winds are avilable to flexpart
             and therefore are wrong from this point. For this reason we label them from their last valid position."""
            if np.max(hyb)> 100 :
                selec = hyb>100
                nr = radada(datpart['itime'],
                        datpart['xf'][selec],datpart['yf'][selec],datpart['pf'][selec],
                        datpart['tempf'][selec],datpart['idx_back'][selec],
                        prod0['flag_source'],prod0['src']['x'],prod0['src']['y'],
                        prod0['src']['p'],prod0['src']['t'],prod0['src']['age'],
                        part0['ir_start'])
                nradada += nr
                nhits[0] += nr
           
            """ PROCESS THE (ADJOINT) DETRAINMENT """
            n1 = detrainer(datpart['itime'], 
                datpart['xi'],datpart['yi'],datpart['pi'],datpart['tempi'],hyb,
                datpart['xf'],datpart['yf'], datrean.var['UDR'], datpart['idx_back'],\
                prod0['flag_source'],part0['ir_start'], prod0['chi'],prod0['passed'],\
                prod0['src']['x'],prod0['src']['y'],prod0['src']['p'],prod0['src']['t'],\
                prod0['src']['age'],prod0['source'],prod0['pl'],\
                datrean.attr['Lo1'],datrean.attr['La1'],datrean.attr['dlo'],datrean.attr['dla'],\
                detr_offset,segl,mask)
            nhits += n1
            #@@ test
            print('return from detrainer',nhits)
            #@@ end test    
            sys.stdout.flush()

        """ End of of loop on slices """
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
            
        print('end hour ',hour,'  numact', partpost['nact'], ' nexits',nexits,' nhits',nhits)
        print('nlive', nlive,' nprist',nprist,' nouthit',nouthit,' noutprist',noutprist,' nhitpure',nhitpure)
        # check that nlive + nhits + nexits = numpart, should be true after the first day
        if partpost['nact'] != nprist + nouthit + noutprist + nhitpure + ndborne:
            print('@@@ ACHTUNG numact not equal to sum ',partpost['nact'],nprist + nouthit + noutprist + nhitpure + ndborne)

    """ End of the procedure and storage of the result """
    # clear other fields before the output as this step requires memory
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use before clean: {:4.2f} gb'.format(memoryUse))
    del partante
    del partpost
    del partStep
    del datpart
    del datrean
    del hyb
    del live_a
    del live_p
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2**30
    print('memory use after clean: {:4.2f} gb'.format(memoryUse))
    sys.stdout.flush()
    # reduction of the size of prod0 by converting float64 into float32
    prod0['chi'] = prod0['chi'].astype(np.float32)
    prod0['passed'] = prod0['passed'].astype(np.int32)
    for var in ['age','p','t','x','y']:
        prod0['src'][var] = prod0['src'][var].astype(np.float32)
    #output file
    dd.io.save(out_file2,prod0,compression='zlib')
    #pickle.dump(prod0,gzip.open(out_file,'wb'))
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
           dat['idx_back'] = dat['xi'] = dat['yi'] = dat['pi'] = dat['ti'] = []
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
                dat['xf'] = np.clip(part_a['x'][live_a],-10.,160.)
                dat['yf'] = np.clip(part_a['y'][live_a],0.,50.)
                dat['pf'] = np.clip(part_a['p'][live_a],0.,50.)
                dat['tempf'] = np.clip(part_a['t'][live_a],0.,50.)
                #dat['pf'] = part_a['p'][live_a]
                #dat['tempf'] = part_a['t'][live_a]
                dat['xi'] = np.clip(coefai*part_a['x'][live_a] + coefpi*part_p['x'][live_p],-10.,160.)
                dat['yi'] = np.clip(coefai*part_a['y'][live_a] + coefpi*part_p['y'][live_p],0.,50.)
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
                dat['xi'] = np.clip(part_p['x'][live_p],-10.,160.)
                dat['yi'] = np.clip(part_p['y'][live_p],0.,50.)
                dat['pi'] = part_p['p'][live_p]
                dat['tempi'] = part_p['t'][live_p]
            else:
                dat['xf'] = dat['xi']
                dat['yf'] = dat['yi']
                dat['pf'] = dat['pi']
                dat['tempf'] = dat['tempi']
                #dat['pf'] = dat['pi']
                #dat['tempf'] = dat['ti']
                dat['xi'] = np.clip(coefai*part_a['x'][live_a] + coefpi*part_p['x'][live_p],-10.,160.)
                dat['yi'] = np.clip(coefai*part_a['y'][live_a] + coefpi*part_p['y'][live_p],0.,50.)
                dat['pi'] = coefai*part_a['p'][live_a] + coefpi*part_p['p'][live_p]
                dat['tempi'] = coefai*part_a['t'][live_a] + coefpi*part_p['t'][live_p]
        tf -= slice_width
        yield dat

#%%
""" Function managing the exiting parcels """

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
            if   y[i] < rr[1,0] + 4.: excode = 6
            elif x[i] < rr[0,0] + 4.: excode = 3
            elif y[i] > rr[1,1] - 4.: excode = 4
            elif x[i] > rr[0,1] - 4.: excode = 5
            elif p[i] > highpcut - 150: excode = 1
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

""" Function finding the detrainment at the location of the parcel and doing the job """

@jit(nopython=True,cache=True)
def detrainer(itime, xi,yi,pi,ti,hyb,xf,yf, udr, idx_back,flag,ir_start,chi,passed,\
              xc,yc,pc,tc,age,source,pl,\
              Lo1,La1,dlo,dla,detr_offset,segl,mask):
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
                # Contribution to the source distribution
                source[ym,xm] += chi[i0] - newchi
                # Contribution to the boxed source distribution
                pl[mask[ym,xm],int(i0/segl)] += chi[i0] - newchi
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
