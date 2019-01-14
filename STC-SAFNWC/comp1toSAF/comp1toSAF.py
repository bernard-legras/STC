#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare SAF and cloudtop stats.
This is based on convscr1 and consrcSAF, that is on the getsat1 and getsatSAF
which have been extracted from the former two. 

This script reads cloud top data from the CTTH files produced by the SAFNWC code
and cloudtop files for both MSG1 and Himawari. The data are retained for the satellite
images that just preceeded the time slots 0-2-4-...-20-22 h (can be changed with the parameter dstep).
The data are then projected onto the FULLAMA grid at 0.1° resolution and an output is produced
2h as a pickle file. The longitudes in the 10W-90E range are filled by MSG1 and the longitudes in 
90E-160E range are filled by Himawari. Notice that the cloudtop data for MSG1 start at 6E and not 10W.

The output contains two pixmap objects, one for the cloudtop data and one for the SAFNWC data. 

Modification to retro apply to convsrc: 
    - pre mode that corrects the shift between the interval and the data.
    - dtRange defined in a single location in the pixmap initialization

Modification to be done to convsrc
    - chage opaq policy to a policy based on the CT product

Modification not to retro apply
    - SAF: opaq is coded in a flag field rather than used to mask the data
    - cloudtop: the creation flag copied in the flag field
    
The inconvenience of this program is that it reads all the files even if only 
a subset needs to be processed. This cannot be changed easily without loosing
the main structure.

Created on Wed Dec 26 02:22:51 2018

@author: Bernard Legras
"""
import socket
import numpy as np
from collections import defaultdict
from numba import jit
from datetime import datetime, timedelta
import os
import sys
#import psutil
#import deepdish as dd
import pickle,gzip
import argparse
import SAFNWCnc
from io107 import readidx107
import geosat
import constants as cst

# misc parameters
# step in the cloudtop procedure
cloudtop_step = timedelta(hours=12)

# if True print a lot of junk
verbose = True
debug = False
debug2 = False

# Change these two values to False to process the stanfdar 2016 SAF product 
# Determines whether we use the reprocessed data or not (Std mode)
NewSAF = False
# Skip the processing of Himawari
doHima = False

#%%
"""@@@@@@@@@@@@@@@@@@@@@@@@   MAIN   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

def main():
    
    # to be updated
    if socket.gethostname() == 'graphium':
        pass
    elif 'ciclad' in socket.gethostname():
        #main_sat_dirSAF = '/data/legras/flexpart_in/SAFNWC'
        main_sat_dir1 = '/bdd/STRATOCLIM/flexpart_in'
        output_dir = '/data/legras/STC/pixmaps'
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        #main_sat_dirSAF = '/data/legras/flexpart_in/SAFNWC'
        main_sat_dir1 = '/bdd/STRATOCLIM/flexpart_in'
    else:
         print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')
         exit()
         
    parser = argparse.ArgumentParser()
    #parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m1","--month1",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-d1","--day1",type=int,choices=1+np.arange(31),help="day")
    parser.add_argument("-m2","--month2",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-d2","--day2",type=int,choices=1+np.arange(31),help="day")
    #parser.add_argument("-N","--NewSAF",type=bool,help='NewSAF')
    #parser.add_argument("-H","--doHima",type=bool,help='doHima')

    """ Parameters """
    # Starting time
    # This date must fall on a cloudtop boundary, that it 0 or 12h
    month1 = 8
    day1 = 16
    #sdate = datetime(2017,month1,day1,0)
    #sdate = datetime(2017,7,21,0)
    # Ending time (backward)
    month2 = 7
    day2 = 13
    #edate = datetime(2017,7,13,0)
    # time between two comparisons
    dstep = timedelta (hours=2)
    # dtRange
    #dtRangeSAF = {'MSG1':timedelta(minutes=15),'Hima':timedelta(minutes=20)}
    #dtRange1   = {'MSG1':timedelta(minutes=30),'Hima':timedelta(minutes=20)}
    # SAF filters
    #good = True
    #opaq = True
    
    args = parser.parse_args()
    if args.month1 is not None: month1 = args.month1
    if args.month2 is not None: month2 = args.month2
    if args.day1 is not None: day1 = args.day1
    if args.day2 is not None: day2 = args.day2
    #if args.NewSAF is not None: NewSAF = args.NewSAF
    #if args.doHima is not None: doHima = args.doHima
    
    sdate = datetime(2017,month1,day1,0)
    edate = datetime(2017,month2,day2,0)
    
    # Directories for the satellite cloud top files
    #satdirSAF ={'MSG1':os.path.join(main_sat_dirSAF,'msg1','S_NWC'),\
    #            'Hima':os.path.join(main_sat_dirSAF,'himawari','S_NWC')}
    satdir1 = {'MSG1':os.path.join(main_sat_dir1,'StratoClim+1kmD_msg1-c'),\
               'Hima':os.path.join(main_sat_dir1,'StratoClim+1kmD_himawari-d')}

    """ Initialization of the calculation """
    # Initialize the grid
    if NewSAF:
        gg = geosat.GeoGrid('FullAMA_SAFBox')
    else:
        gg = geosat.GeoGrid('FullAMA')
        output_dir += '-StdSAF'
    satmapSAF = pixmap(gg,'SAF')
    satmap1 = pixmap(gg,'1')
    satfill = {}

    # Build the satellite field generator
    get_sat1 = {'MSG1': read_sat1(sdate,satmap1.zone['MSG1']['dtRange1'],satdir1['MSG1'],pre=True),\
                'Hima': read_sat1(sdate,satmap1.zone['Hima']['dtRange1'],satdir1['Hima'],pre=True)}
    #get_satSAF = {'MSG1': read_satSAF(sdate,'MSG1',satmapSAF.zone['MSG1']['dtRangeSAF'],satdirSAF['MSG1'],pre=True),\
    #              'Hima': read_satSAF(sdate,'Hima',satmapSAF.zone['Hima']['dtRangeSAF'],satdirSAF['Hima'],pre=True)}
    get_satSAF = {'MSG1': read_satSAF(sdate,'MSG1',satmapSAF.zone['MSG1']['dtRangeSAF'],pre=True),\
                  'Hima': read_satSAF(sdate,'Hima',satmapSAF.zone['Hima']['dtRangeSAF'],pre=True)}
    current_date = sdate

    """ Main loop on the output time steps """
    while current_date > edate:
        print('searching ',current_date) 
        
        # Fill the pixmap for cloudtops
        
        while satmap1.check('MSG1',current_date) is False:
            # if not get next satellite slice
            try:
                next(satfill['MSG1'])
            # read new satellite file if the slice generator is over
            except:
                datsatm = next(get_sat1['MSG1'])
                satfill['MSG1'] = satmap1.fill1('MSG1',datsatm)
                next(satfill['MSG1'])               
        if verbose: print('check 1 MSG1 ',satmap1.check('MSG1',current_date),'##',current_date,
                      '##',satmap1.zone['MSG1']['ti'],'##',satmap1.zone['MSG1']['tf'])
        if doHima:              
            while satmap1.check('Hima',current_date) is False:
                try:
                    next(satfill['Hima'])
                except:
                    datsath = next(get_sat1['Hima'])
                    satfill['Hima'] = satmap1.fill1('Hima',datsath)
                    next(satfill['Hima'])
            if verbose: print('check 1 Hima ',satmap1.check('Hima',current_date),'##',current_date,
                          '##',satmap1.zone['Hima']['ti'],'##',satmap1.zone['Hima']['tf'])

        # Fill the pixmap for SAF
        # Check whether the present satellite image is valid
        # The while should ensure that the run synchronizes when it starts.
            
        while satmapSAF.check('MSG1',current_date) is False:
            # if not get next satellite image 
            datsatm = next(get_satSAF['MSG1'])
            # Check that the image is available
            if datsatm is not None:
                # all the data need to be read because of possible gaps 
                # and extensions to be done in such cases
                pmm = geosat.SatGrid(datsatm,gg)
                pmm._sat_togrid('CTTH_PRESS')
                #print('pm1 diag',len(datsat1.var['CTTH_PRESS'][:].compressed()),
                #                 len(pm1.var['CTTH_PRESS'][:].compressed()))
                pmm._sat_togrid('ctth_alti')
                pmm._sat_togrid('ctth_tempe')
                pmm._sat_togrid('ctth_quality')
                pmm._sat_togrid('ctth_conditions')
                pmm._sat_togrid('ctth_status_flag')
                pmm._sat_togrid('ctth_method')
                pmm._sat_togrid('CT')
                pmm._sat_togrid('ct_cumuliform')
                pmm._sat_togrid('ct_multilayer')
                pmm._sat_togrid('ct_status_flag')
                pmm._sat_togrid('CMa')
                pmm._sat_togrid('cma_status_flag')
                pmm._sat_togrid('cma_quality')
                pmm.attr = datsatm.attr.copy()
                satmapSAF.fillSAF('MSG1',pmm)
                del pmm
                #del datsat1
            else:
                # if the image is missing, extend the lease of previous image
                try:
                    satmapSAF.extend('MSG1')
                    if debug2: print ('MSG1 extension ',satmapSAF.zone['MSG1']['ti'])
                except:
                    # This handle the unlikely case where the first image is missing
                    continue
        #print('match MSG1 ',datsatm.attr['date'],datsatm.attr['lease_time'])
        if verbose: print('check SAF MSG1 ',satmapSAF.check('MSG1',current_date),'##',current_date,
                      '##',satmapSAF.zone['MSG1']['ti'],'##',satmapSAF.zone['MSG1']['tf'])
                     
        if doHima:                
            while satmapSAF.check('Hima',current_date) is False:
                # if not get next satellite image 
                datsath = next(get_satSAF['Hima'])
                # Check that the image is available
                if datsath is not None:
                    # not elegant as it duplicates the while test but saves 
                    # some processing
                    #if (current_date <= datsath.attr['lease_time']) | \
                    #   (current_date > datsath.attr['date']): continue
                    # not a good odea when some image is missing
                    pmh = geosat.SatGrid(datsath,gg)
                    pmh._sat_togrid('CTTH_PRESS')
                    #print('pmh diag',len(datsath.var['CTTH_PRESS'][:].compressed()),
                    #                 len(pmh.var['CTTH_PRESS'][:].compressed()))
                    pmh._sat_togrid('ctth_tempe')
                    pmh._sat_togrid('ctth_alti')
                    pmh._sat_togrid('ctth_quality')
                    pmh._sat_togrid('ctth_conditions')
                    pmh._sat_togrid('ctth_status_flag')
                    pmh._sat_togrid('ctth_method')
                    pmh._sat_togrid('CT')
                    pmh._sat_togrid('ct_cumuliform')
                    pmh._sat_togrid('ct_multilayer')
                    pmh._sat_togrid('ct_status_flag')
                    pmh._sat_togrid('CMa')
                    pmh._sat_togrid('cma_status_flag')
                    pmh._sat_togrid('cma_quality')
                    pmh.attr = datsath.attr
                    satmapSAF.fillSAF('Hima',pmh)
                    #del datsath
                    del pmh
                else:
                    # if the image is missing, extend the leaseof previous image
                    try:
                        satmapSAF.extend('Hima')
                        if debug2: print ('Hima extension ',satmapSAF.zone['Hima']['ti'])
                    except:
                        # This handle the unlikely case where the first image is missing
                        continue
            #print('match Hima ',datsath.attr['date'],datsath.attr['lease_time'])      
            if verbose: print('check SAF Hima ',satmapSAF.check('Hima',current_date),'##',current_date,
                          '##',satmapSAF.zone['Hima']['ti'],'##',satmapSAF.zone['Hima']['tf'])
        
        outfile = os.path.join(output_dir,current_date.strftime('%Y-%m-%dT%H%M'))
        with gzip.open(outfile,'wb') as f:
           pickle.dump([satmapSAF,satmap1],f)
        print('OUTPUT for time ',current_date,'\n')
        current_date -= dstep
        sys.stdout.flush() 

"""@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""


#%%
""" Functions related to satellite read """

def read_sat1(t0,dtRange,satdir,pre=False):
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
        if debug: print('nt ',dat['nt'],len(id))
        tf = current_time + cloudtop_step/2 - dt
        # Generate list of time intervals
        # if pre (default), the validity interval follows the time of the satellite image
        # if not pre the validity interval is before  
        dat['time'] = []
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
        if debug: print('check sorting ',np.sum(dat['numRange']),dat['numpart'])
        # iterate time
        current_time -= cloudtop_step
        yield dat
       
def read_satSAF(t0,sat,dtRange,pre=False):
    """ Generator reading the satellite data.
    The loop is infinite; sat data are called when required until the end of
    the parcel loop. """              
    # get dt from satmap
    dt = dtRange
    # initial time
    current_time = t0
    namesat={'MSG1':'msg1','Hima':'himawari'}
    while True:
        #fname = os.path.join(satdir,current_time.strftime('%Y/%Y_%m_%d'))
        #if sat=='MSG1'
            #fname = os.path.join(fname,current_time.strftime('S_NWC_CTTH_MSG1_FULLAMA-VISIR_%Y%m%dT%H%M00Z.nc'))
        #elif sat=='Hima':
            #fname = os.path.join(fname,current_time.strftime('S_NWC_CTTH_HIMAWARI08_FULLAMA-NR_%Y%m%dT%H%M00Z.nc'))
        #else:
        #    print('sat should be MSG1 or Hima')
        #    return
        try:
            if NewSAF:             
                dat = SAFNWCnc.SAFNWC_CTTH(current_time,namesat[sat],BBname='SAFBox')
                dat_ct = SAFNWCnc.SAFNWC_CT(current_time,namesat[sat],BBname='SAFBox')
                dat_cma = SAFNWCnc.SAFNWC_CMa(current_time,namesat[sat],BBname='SAFBox')
            else:
                dat = SAFNWCnc.SAFNWC_CTTH(current_time,namesat[sat])
                dat_ct = SAFNWCnc.SAFNWC_CT(current_time,namesat[sat])
                dat_cma = SAFNWCnc.SAFNWC_CMa(current_time,namesat[sat])
            dat._CTTH_PRESS()
            # This pressure i left in hPa to allow masked with the fill_value in sat_togrid
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
            dat._get_var('ctth_alti')
            dat._get_var('ctth_tempe')
            dat._get_var('ctth_status_flag')
            dat._get_var('ctth_conditions')
            dat._get_var('ctth_quality')
            dat._get_var('ctth_method')
            dat.close()
            dat_ct._CT()
            dat_ct._get_var('ct_cumuliform')
            dat_ct._get_var('ct_multilayer')
            dat_ct._get_var('ct_status_flag')
            dat.var['CT'] = dat_ct.var['CT']
            dat.var['ct_cumuliform'] = dat_ct.var['ct_cumuliform']
            dat.var['ct_multilayer'] = dat_ct.var['ct_multilayer']
            dat.var['ct_status_flag'] = dat_ct.var['ct_status_flag']
            dat_ct.close()
            dat_cma._CMa()
            dat_cma._get_var('cma_status_flag')
            dat_cma._get_var('cma_quality')
            dat.var['CMa'] = dat_cma.var['CMa']
            dat.var['cma_status_flag'] = dat_cma.var['cma_status_flag']
            dat.var['cma_quality'] = dat_cma.var['cma_quality']
            dat_cma.close()
        except FileNotFoundError:
            print('SAF file not found ',current_time,namesat[sat])
            dat = None
        current_time -= dtRange
        yield dat
        
#%%
""" Describe the pixel map that contains slice of cloudtop data """


class pixmap(geosat.GridField):

    def __init__(self,gg,typ):
        
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
        # ACHTUNG: these parameters should not be defined twice
        self.zone['MSG1']['dtRangeSAF'] = timedelta(minutes=15)
        self.zone['Hima']['dtRangeSAF'] = timedelta(minutes=20)
        self.zone['MSG1']['dtRange1'] = timedelta(minutes=30)
        self.zone['Hima']['dtRange1'] = timedelta(minutes=20)
        # define the slice
        self.ptop = np.empty(shape=self.geogrid.shapeyx,dtype=np.float)
        self.ptop.fill(cst.p0)
        if typ == 'SAF':
            self.flag_ctth = np.zeros(shape=self.geogrid.shapeyx,dtype=np.uint64)
            self.flag_ct = np.zeros(shape=self.geogrid.shapeyx,dtype=np.uint32)
            self.flag_cma = np.zeros(shape=self.geogrid.shapeyx,dtype=np.uint32)
            self.tempe = np.empty(shape=self.geogrid.shapeyx,dtype=np.float)
            self.alti = np.empty(shape=self.geogrid.shapeyx,dtype=np.float)
        elif typ == '1':
            self.flag = np.zeros(shape=self.geogrid.shapeyx,dtype=np.uint32)
        self.num  = np.zeros(shape=self.geogrid.shapeyx,dtype=np.int32)
        self.range = self.geogrid.box_range
        self.binx = self.geogrid.box_binx
        self.biny = self.geogrid.box_biny
        self.stepx = (self.range[0,1]-self.range[0,0])/self.binx
        self.stepy = (self.range[1,1]-self.range[1,0])/self.biny
        print('steps',self.stepx,self.stepy)
    
    def set_mask(self):
        # define the regional mask of the pixmap
        pass

    def erase(self,zone):
        # "ground" temperature
        T0 = 310
        # erase the data in the zone
        x1 = self.zone[zone]['xi']
        x2 = x1 + self.zone[zone]['binx']
        y1 = self.zone[zone]['yi']
        y2 = y1 + self.zone[zone]['biny']
        self.ptop[y1:y2,x1:x2].fill(cst.p0)
        self.num[y1:y2,x1:x2].fill(0)
        try:
            self.flag[y1:y2,x1:x2].fill(0)
        except:
            self.flag_ct[y1:y2,x1:x2].fill(0)
            self.flag_cma[y1:y2,x1:x2].fill(0)
            self.flag_ctth[y1:y2,x1:x2].fill(0)
            self.alti[y1:y2,x1:x2].fill(0)
            self.tempe[y1:y2,x1:x2].fill(T0)
        if debug:
             print('erase ',zone,x1,x2)

    def check(self,zone,t):
         # check that the zone is not expired
        try:
             test = (t > self.zone[zone]['ti']) and (t <= self.zone[zone]['tf'])
         # Exception for the first usage when the time keys are not defined
        except KeyError:
             test = False
        return test

    def extend(self,zone):
        self.zone[zone]['ti'] -= self.zone[zone]['dtRangeSAF']
      
    def fillSAF(self,zone,dat):
        """ Function filling the zone with new data from the satellite dictionary.
        """
        # Erase the zone
        self.erase(zone)
        # Mask outside the new data outside the zone
        if zone == 'MSG1':
            dat.var['CTTH_PRESS'][:,self.zone['Hima']['xi']:] = np.ma.masked
        elif zone == 'Hima':
            dat.var['CTTH_PRESS'][:,:self.zone['MSG1']['binx']] = np.ma.masked
        #nbValidBeforeSel = len(dat.var['CTTH_PRESS'].compressed())
#        dat.var['OPAQ'] = np.zeros(shape=dat.var['CTTH_PRESS'].shape,dtype=np.int32)
#        # Filter according to quality keeping only good retrievals with all input fields
#        if good:
#            sel = ((dat.var['ctth_quality']&0x38)==0x8) & ((dat.var['ctth_conditions']&0xFF00)==0x5500)
#            #dat.var['CTTH_PRESS'][~sel] = np.ma.masked
#            dat.var['OPAQ'][sel] += 2
#        # Filter the non opaque clouds if required
#        if opaq:
#            sel = (dat.var['ctth_status_flag']&4)==4
#            #dat.var['CTTH_PRESS'][~sel] = np.ma.masked
#            dat.var['OPAQ'][sel] += 1
        # diag for CTTH, CT, CMA
        filt = ~dat.var['CTTH_PRESS'].mask
#        print('CTTH method, conditions, status, quality ',dat.var['ctth_method'][filt].max(),\
#              dat.var['ctth_conditions'][filt].max(),dat.var['ctth_status_flag'][filt].max(),\
#              dat.var['ctth_quality'][filt].max())
#        print('CT CT cumuliform multilayer status ',dat.var['CT'][filt].max(),\
#              dat.var['ct_cumuliform'][filt].max(),dat.var['ct_multilayer'][filt].max(),\
#              dat.var['ct_status_flag'][filt].max())
#        print('CMa CMa satus quality ',dat.var['CMa'][filt].max(),\
#              dat.var['cma_status_flag'][filt].max(),dat.var['cma_quality'][filt].max())
        FLAG_CTTH = (dat.var['ctth_method'].astype(np.uint64) << 48) \
                  + (dat.var['ctth_conditions'].astype(np.uint64) << 32) \
                  + (dat.var['ctth_status_flag'].astype(np.uint64) << 16) \
                  +  dat.var['ctth_quality']
        FLAG_CT = (dat.var['CT'].astype(np.uint32) << 24) \
                + (dat.var['ct_cumuliform'].astype(np.uint32) << 16) \
                + (dat.var['ct_multilayer'].astype(np.uint32) << 8) \
                +  dat.var['ct_status_flag']
        FLAG_CMA = (dat.var['CMa'].astype(np.uint32) << 24) \
                 + (dat.var['cma_status_flag'].astype(np.uint32) << 8) \
                 + dat.var['cma_quality']
#        print('FLAGS ',FLAG_CTTH[filt].max(),FLAG_CT[filt].max(),FLAG_CMA[filt].max())
        # test : count the number of valid pixels
        #nbValidAfterSel = len(dat.var['CTTH_PRESS'].compressed())
        #print('valid pixels before & after selection',zone,nbValidBeforeSel,nbValidAfterSel)
        # Inject the non masked new data in the pixmap
        # Conversion to Pa is done here
        self.ptop[filt] = 100*dat.var['CTTH_PRESS'][filt]
        self.tempe[filt] = dat.var['ctth_tempe'][filt]
        self.alti[filt] = dat.var['ctth_alti'][filt]
        self.flag_ctth[filt] = FLAG_CTTH[filt]
        self.flag_ct[filt] = FLAG_CT[filt]
        self.flag_cma[filt] = FLAG_CMA[filt]
        # set the new expiration date 
        self.zone[zone]['tf'] = dat.attr['date']
        self.zone[zone]['ti'] = dat.attr['lease_time']
        
        if debug:
            sel = self.ptop < cst.p0
            nbact = 100*sel.sum()/(self.geogrid.box_binx*self.geogrid.box_biny)
            if verbose: print('fill ',zone,' #selec ',len(sel),' % {:4.2f}'.format(nbact),\
                ' meanP {:6.0f} minP {:5.0f}'.format(self.ptop[sel].mean(),self.ptop.min()))
        return

    def fill1(self,zone,dat):
        """ Generator filling the slice with new data from the satellite dictionary
        of cloudtop pixels, using the precalculated time ranges.
        The data are read from the end to the beginning to fit the backward analysis
        The indexRange is such that the first value point to the end of the last
        range, the second to the end of the last-1 range and so on
        """
        for i in reversed(range(dat['nt'])):
            self.erase(zone)
            self.zone[zone]['tf'] = dat['time'][i][1]
            self.zone[zone]['ti'] = dat['time'][i][0]
            if debug:
                print('fill1 ',i,dat['numRange'][i])
            if dat['numRange'][i] >0:
                selec = range(dat['indexRange'][i,0],dat['indexRange'][i,1])

                #idx = np.floor((dat[zone]['x'][selec] - self.range[0,0])/self.stepx).astype('int')
                #idy = np.floor((dat[zone]['y'][selec] - self.range[1,0])/self.stepy).astype('int')
                fillfast(self.ptop,self.num,self.flag,dat['x'][selec], \
                         dat['y'][selec],dat['p'][selec], dat['flag'][selec], \
                         self.range[0,0],self.range[1,0],self.stepx,self.stepy,self.binx,self.biny)
                if debug:
                    sel = self.ptop < cst.p0
                    nbact = 100*sel.sum()/(self.binx*self.biny)
                    if verbose: print('fill ',zone,' #selec ',len(selec),' % {:4.2f}'.format(nbact),\
                          ' meanP {:6.0f} minP {:5.0f}'.format(self.ptop[sel].mean(),self.ptop.min()),\
                          ' nmax ',self.num.max(),\
                          ' minx {:7.2f} maxx {:7.2f}'.format(dat['x'][selec].min(),dat['x'][selec].max()))
            yield i

@jit(nopython=True)
def fillfast(pSlice,numSlice,Flag,x,y,p,dflag,x0,y0,stepx,stepy,binx,biny):
    for i in range(len(x)):
        idx = min(int(np.floor((x[i]-x0)/stepx)),binx-1)
        idy = min(int(np.floor((y[i]-y0)/stepy)),biny-1)
        numSlice[idy,idx] += 1
        if p[i] < pSlice[idy,idx]:
            pSlice[idy,idx] = p[i]
            Flag[idy,idx] = dflag[i]

if __name__ == '__main__':
    main()
