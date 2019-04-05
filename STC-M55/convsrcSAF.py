# -*- coding: utf-8 -*-
"""
Main code to analyse the convective sources of the air sampled during the StratoClim
campaign.
This version is based on the SAFNWC product which is availabla at the same time and same 
resolution as the satellite data.
We use the reprocessed version which is produced on atrunacted image.

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
import SAFNWCnc
import geosat

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
    parser.add_argument("-t","--step",type=int,help="step in hour between two part files")
    parser.add_argument("-f","--flight",type=str,help="flight identifier for balloons")
    parser.add_argument("-ct","--cloud_type",type=str,choices=["meanhigh","veryhigh","silviahigh"],help="cloud type filter")
    parser.add_argument("-k","--diffus",type=str,choices=['01','1','001'],help='diffusivity parameter')
    
    # to be updated
    if socket.gethostname() == 'graphium':
        pass
    elif 'ciclad' in socket.gethostname():
        #root_dir = '/home/legras/STC/STC-M55'
        main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
        traj_dir = '/data/legras/flexout/STC/M55'
        out_dir = '/data/legras/STC'
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        #root_dir = '/home/legras/STC/STC-M55'
        print('does not rune on climserv')
        exit()
    elif socket.gethostname() == 'satie':
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
    slice_width = timedelta(minutes=5)
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)
    # default values of parameters
    # date of the flight
    year=2017
    month=7
    day=29
    platform = 'GLO'
    advect = 'EAZ'
    suffix =''
    launch_number=''
    flight = 'KU005'
    quiet = False
    clean0 = False
    cloud_type = 'silviahigh'
    diffus = '01'
    args = parser.parse_args()
    if args.year is not None: year=args.year
    if args.month is not None: month=args.month
    if args.day is not None: day=args.day
    if args.advect is not None: advect=args.advect
    if args.platform is not None: platform=args.platform
    if args.launch_number is not None: launch_number='-'+str(args.launch_number)
    if args.suffix is not None: suffix='-'+args.suffix
    if args.quiet is not None:
        if args.quiet=='y': quiet=True
        else: quiet=False
    if args.clean0 is not None: clean0 = args.clean0
    if args.cloud_type is not None: cloud_type = args.cloud_type
    if args.step is not None: step = args.step
    if args.flight is not None: flight = args.flight
    if args.diffus is not None: diffus = args.diffus
    diffus = '-D' + diffus

    # Update the out_dir with the cloud type
    out_dir = os.path.join(out_dir,'STC-M55-OUT-SAF-'+cloud_type)
    fdate = datetime(year,month,day)

    # Manage the file that receives the print output
    if quiet:
        # Output file
        if platform == 'BAL':
           print_file = os.path.join(out_dir,'out',flight+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+diffus+suffix+'.out')
        else:  
           print_file = os.path.join(out_dir,'out',platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+diffus+suffix+'.out')
        fsock = open(print_file,'w')
        sys.stdout=fsock

    # initial time to read the sat files
    # should be after the end of the flight
    sdate = fdate + timedelta(days=1)
    print('year',year,'month',month,'day',day)
    print('advect',advect)
    print('platform',platform)
    print('launch_number',launch_number)
    print('suffix',suffix)
    
    # patch to fix the problem of the data hole on the evening of 2 August 2017
    if fdate == datetime(2017,8,2): sdate = datetime(2017,8,3,6)

    # Directories of the backward trajectories and name of the output file
    if platform == 'BAL':
        traj_dir = '/data/legras/flexout/STC/BAL'
        ftraj = os.path.join(traj_dir,flight+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+diffus+suffix)
        out_file2 = os.path.join(out_dir,flight+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+diffus+suffix+'.hdf5')
    else:    
        ftraj = os.path.join(traj_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+diffus+suffix)
        out_file2 = os.path.join(out_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+diffus+suffix+'.hdf5')

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
    get_sat = {'MSG1': read_sat(sdate,'MSG1',satmap.zone['MSG1']['dtRange'],satdir['MSG1'],pre=True),\
               'Hima': read_sat(sdate,'Hima',satmap.zone['Hima']['dtRange'],satdir['Hima'],pre=True)}

    # Read the index file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'index_old'),quiet=True)
    print('numpart',part0['numpart'])
    # stamp_date not set in these runs
    # current_date actually shifted by one day / sdate
    # We assume here that part time is defined from this day at 0h
    # sdate defined above should be equal or posterior to current_date
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
        if (hour <= 30) & (partpost['nact']>0):
            new[partpost['idx_back'][~kept_p]-IDX_ORGN] = True
            nnew += len(partpost['x'])-kept_p.sum()
        if hour == 30:
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
                    pm1._sat_togrid('CT')
                    pm1.attr = datsat1.attr.copy()
                    satmap.fill('MSG1',pm1,cloud_type)
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
                    pmh._sat_togrid('CT')
                    pmh.attr = datsath.attr.copy()
                    satmap.fill('Hima',pmh,cloud_type)
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

def read_sat(t0,sat,dtRange,satdir,pre=False):
    """ Generator reading the satellite data.
    The loop is infinite; sat data are called when required until the end of
    the parcel loop. """
    # get dt from satmap
    dt = dtRange
    # initial time
    current_time = t0
    namesat={'MSG1':'msg1','Hima':'himawari'}
    while True:
        # why do we need that since fname is not used, scoria of the previous non SAF version
        fname = os.path.join(satdir,current_time.strftime('%Y/%Y_%m_%d'))
        if sat=='MSG1':
            fname = os.path.join(fname,current_time.strftime('S_NWC_CTTH_MSG1_FULLAMA-VISIR_%Y%m%dT%H%M00Z.nc'))
        elif sat=='Hima':
            fname = os.path.join(fname,current_time.strftime('S_NWC_CTTH_HIMAWARI08_FULLAMA-NR_%Y%m%dT%H%M00Z.nc'))
        else:
            print('sat should be MSG1 or Hima')
            return
        try:
            dat = SAFNWCnc.SAFNWC_CTTH(current_time,namesat[sat],BBname='SAFBox')
            dat_ct = SAFNWCnc.SAFNWC_CT(current_time,namesat[sat],BBname='SAFBox')
            dat._CTTH_PRESS()
            dat_ct._CT()
            dat.var['CT'] = dat_ct.var['CT']
            # This pressure is left in hPa to allow masked with the fill_value in sat_togrid
            # The conversion to Pa is made in fill
            dat.attr['dtRange'] = dt
             # if pre, the validity interval follows the time of the satellite image
            # if not pre (default) the validity interval is before 
            if pre:
               dat.attr['lease_time'] = current_time 
               dat.attr['date'] = current_time + dtRange
            else:
               dat.attr['lease_time'] = current_time - dtRange
               dat.attr['date'] = current_time
            dat.close()
            dat_ct.close()
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

    def fill(self,zone,dat,cloud_type):
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
        # Filter according to cloud_type
        if cloud_type == 'meanhigh':
            sel = (dat.var['CT'] ==8) |  (dat.var['CT'] ==9) | (dat.var['CT'] == 12) | (dat.var['CT'] == 13)
            dat.var['CTTH_PRESS'][~sel] = np.ma.masked
        elif cloud_type == 'veryhigh':
            sel = (dat.var['CT'] ==9) | (dat.var['CT'] == 13)
            dat.var['CTTH_PRESS'][~sel] = np.ma.masked
        elif cloud_type == 'silviahigh':
            sel = (dat.var['CT'] ==9) | (dat.var['CT'] == 13) | (dat.var['CT'] == 8)
            dat.var['CTTH_PRESS'][~sel] = np.ma.masked
        # test : count the number of valid pixels
        nbValidAfterSel = len(dat.var['CTTH_PRESS'].compressed())
        print('valid pixels before & after selection',zone,nbValidBeforeSel,nbValidAfterSel)
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
