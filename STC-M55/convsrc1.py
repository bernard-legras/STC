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

@author: Bernard Legras
"""

import socket
import numpy as np
from collections import defaultdict
from numba import jit
from datetime import datetime, timedelta
import os
import deepdish as dd
import pickle, gzip
import sys
import argparse
import psutil

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

    # to be updated
    if socket.gethostname() == 'graphium':
        pass
    elif 'ciclad' in socket.gethostname():
        #root_dir = '/home/legras/STC/STC-M55'
        main_sat_dir = '/bdd/STRATOCLIM/flexpart_in'
        traj_dir = '/data/legras/flexout/STC/M55'
        out_dir = '/data/legras/STC'
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        #root_dir = '/home/legras/STC/STC-M55'
        main_sat_dir = '/bdd/STRATOCLIM/flexpart_in'
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
    slice_width = timedelta(minutes=5)
    # dtRange
    dtRange={'MSG1':timedelta(minutes=30),'Hima':timedelta(minutes=20)}
    # number of slices between two outputs
    nb_slices = int(dstep/slice_width)
    # default values of parameters
    # date of the flight
    year=2017
    month=7
    day=27
    platform = 'M55'
    advect = 'OPZ'
    suffix =''
    launch_number=''
    quiet = False
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

    # Update the out_dir with the platform
    out_dir = os.path.join(out_dir,'STC-'+platform+'-OUT')

    fdate = datetime(year,month,day)

    # Manage the file that receives the print output
    if quiet:
        # Output file
        print_file = os.path.join(out_dir,'out',platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix+'.out')
        saveout = sys.stdout
        fsock = open(print_file,'w')
        sys.stdout=fsock

    # initial time to read the sat files
    # should be after the end of the flight
    # and a 12h or 0h boundary
    sdate = fdate + timedelta(days=1)
    print('year',year,'month',month,'day',day)
    print('advect',advect)
    print('platform',platform)
    print('launch_number',launch_number)
    print('suffix',suffix)

    # Directory of the backward trajectories
    ftraj = os.path.join(traj_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix)

    # Output file
    out_file = os.path.join(out_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix+'.pkl')
    out_file2 = os.path.join(out_dir,platform+fdate.strftime('-%Y%m%d')+launch_number+'-'+advect+'-D01'+suffix+'.hdf5')

    # Directories for the satellite cloud top files
    satdir ={'MSG1':os.path.join(main_sat_dir,'StratoClim+1kmD_msg1-c'),\
             'Hima':os.path.join(main_sat_dir,'StratoClim+1kmD_himawari-d')}

    """ Initialization of the calculation """
    # Initialize the slice map to be used as a buffer for the cloudtops
    satmap = pixmap()
    satfill = {}
    datsat = {}
    # Initialize the dictionary of the parcel dictionaries
    partStep={}

    # Build the satellite field generator
    get_sat = {'MSG1': read_sat(sdate,dtRange['MSG1'],satdir['MSG1']),\
               'Hima': read_sat(sdate,dtRange['Hima'],satdir['Hima'])}

    # Read the index file that contains the initial positions
    part0 = readidx107(os.path.join(ftraj,'index_old'),quiet=True)
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
    prod0['src']['x'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['y'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['p'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['t'] = np.empty(part0['numpart'],dtype='float')
    prod0['src']['age'] = np.empty(part0['numpart'],dtype='int')
    prod0['flag_source'] = part0['flag']
    # truncate eventually to 32 bits at the output stage

    # read the part_000 file
    partStep[0] = readpart107(0,ftraj,quiet=True)

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

def read_sat(t0,dtRange,satdir):
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
        while tf > current_time - cloudtop_step/2 + dt:
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
            while current_time+timedelta(seconds=int(dat['ir_start'][idc])) \
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
        The data are read from the end to the beginning at fit the backward analysis
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
