#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Produces the PTOP files
Makes a composite from Himawari and msg1 SAFNWC cloud tops
in the FULLAMA domain with an image produced each time
a new image is available from Himawari or msg1.
If an image is misssing, the validity of the previous image 
is extended.
A mask is applied for non valid data from the value of ctth only.
The output contains also a flag that is defined from the CT product,
ct_multilayer, ctth_statsu_flag and ctth_quality

This script process the reprocessed  SAF product with 2016 version
algo which is available for both Himawari and MSG1 from 1st May 2017 
to 15 September 2017

CICLAD version

Caveat: does not filter out the corrupted images and therefore this filter
has to be applied at a later stage

Created on 9 January 2019 from compositMHi-NN.py

@author: Bernard Legras
"""
import SAFNWCnc
import geosat
from datetime import datetime, timedelta
import pickle,gzip
import os
import numpy as np
import argparse

delta_msg1 = 15
delta_himawari = 20
delta_time = 5

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--segment",type=str,choices=["Jun.01","Jun.11","Jun.21","Jul.01","Jul.11","Jul.21","Aug.01","Aug.11","Aug.21","Sep.01","Sep.11","Sep.21","OctBeg"],help="segment to be processed")
    seg_dates = {'Jul.01':[datetime(2017,7,1,0),datetime(2017,7,11,0)],
                 'Jul.11':[datetime(2017,7,11,0),datetime(2017,7,21,0)],
                 'Jul.21':[datetime(2017,7,21,0),datetime(2017,8,1,0)],
                 'Aug.01':[datetime(2017,8,1,0),datetime(2017,8,11,0)],
                 'Aug.11':[datetime(2017,8,11,0),datetime(2017,8,21,0)],
                 'Aug.21':[datetime(2017,8,21,0),datetime(2017,9,1,0)],
                 'Jun.01':[datetime(2017,6,1,0),datetime(2017,6,11,0)],
                 'Jun.11':[datetime(2017,6,11,0),datetime(2017,6,21,0)],
                 'Jun.21':[datetime(2017,6,21,0),datetime(2017,7,1,0)],
                 'Sep.01':[datetime(2017,9,1,0),datetime(2017,9,11,0)],
                 'Sep.11':[datetime(2017,9,11,0),datetime(2017,9,21,0)],
                 'Sep.21':[datetime(2017,9,21,0),datetime(2017,10,1,0)],
                 'OctBeg':[datetime(2017,10,1,0),datetime(2017,10,3,0)]}
    
    seg = 'Jul.11'

    args = parser.parse_args()
    if args.segment is not None:
        seg = args.segment
    
    date1 = seg_dates[seg][0]
    date2 = seg_dates[seg][1]
  
    gg=geosat.GeoGrid('FullAMA_SAFBox')

    date = date1

    output_dir = '/data/legras/STC/STC-SAFNWC-OUT'
    
    gM = get_msg1(date1,date2,gg)
    gH = get_himawari(date1,date2,gg)
   
    # date1 should correspond to a date where both images are available
    # in order to start properly
    phM = next(gM)
    phH = next(gH)
    make_new_comp = True
    
    while date < date2:
        
        # Check new msg1 to be loaded
        if phM.expiration == date:
            pp = next(gM)
            # If available image, get it
            # if not,  extend validity of previous image
            if pp is not None:
                phM = pp
                make_new_comp = True
            else:
                phM.expiration += timedelta(minutes=delta_msg1)
        if phH.expiration == date:
            pp = next(gH)
            # If available image, get it
            # if not,  extend validity of previous image
            if pp is not None:
                phH = pp
                make_new_comp = True
            else:
                phH.expiration += timedelta(minutes=delta_himawari)
        
        if make_new_comp:
            comb = phM.patch(phH,90.75,['CTTH_PRESS','CT','ctth_status_flag','ctth_quality','ct_multilayer'])
            comb.var['CTTH_PRESS'][comb.var['CTTH_PRESS']<0]=np.ma.masked
            comb.var['CTTH_PRESS'][comb.var['CTTH_PRESS']>65000]=np.ma.masked
            # compact the CT field and the flags into a single 32 bit variable
            comb.var['cloud_flag'] = (comb.var['CT'].astype(np.uint32) << 24) \
                + (comb.var['ctth_quality'].astype(np.uint32) << 16) \
                + (comb.var['ct_multilayer'].astype(np.uint32) << 8) \
                +  comb.var['ctth_status_flag'].astype(np.uint32)
            comb.var['cloud_flag'].__setmask__(comb.var['CTTH_PRESS'].mask)
            filename = date.strftime('SAFNWC-PTOP-%Y-%m-%d-%H:%M.pkl')
            day_dir = os.path.join(output_dir,date.strftime('%Y/%m/%Y-%m-%d'))
            try:
                os.mkdir(day_dir)
            except:
                pass
            # output containing the cloud top pressure and the good/opaq flag
            fullname = os.path.join(day_dir,filename)
            with gzip.open(fullname,'wb') as f:
                pickle.dump([comb.var['CTTH_PRESS'],comb.var['cloud_flag']],f,protocol=3)
            print ('new file',date,np.ma.count_masked(comb.var['CTTH_PRESS']))
            del comb
                
        date += timedelta(minutes=delta_time)
        make_new_comp = False
    
# Generator providing the msg1 data between two dates 
# The interval is [date1, date2[, including date1 but excluding date2]
def get_msg1(date1,date2,gg):
    current_date = date1
    while current_date < date2:
        try:
            ddM = SAFNWCnc.SAFNWC_CTTH(current_date,'msg1',BBname='SAFBox')
            ddM_ct = SAFNWCnc.SAFNWC_CT(current_date,'msg1',BBname='SAFBox')
            ddM._CTTH_PRESS()
            ddM_ct._CT()
            ddM.var['CT'] = ddM_ct.var['CT']
            ddM._get_var('ctth_status_flag')
            ddM._get_var('ctth_quality')
            ddM_ct._get_var('ct_multilayer')
            ddM.var['ct_multilayer'] = ddM_ct.var['ct_multilayer']
            phM=geosat.SatGrid(ddM,gg)
            phM._sat_togrid('CTTH_PRESS')
            phM._sat_togrid('CT')
            phM._sat_togrid('ctth_status_flag')
            phM._sat_togrid('ctth_quality')
            phM._sat_togrid('ct_multilayer')
            
            # If mask done here, tend to be persistent (side effect of yield?)
            #phM.var['CTTH_PRESS'][phM.var['CTTH_PRESS']>65000.]=np.ma.masked
            ddM.close()
            ddM_ct.close()
            del ddM
            del ddM_ct
            phM.expiration = current_date + timedelta(minutes=delta_msg1)
            yield phM
            del phM
        except AttributeError:
            print('msg1 AttributeError')
            yield -1
        except:
            print('MISSING msg1',current_date)
            yield None
        finally:
            current_date += timedelta(minutes=delta_msg1)
            
# Generator providing the himawari data between two dates 
# The interval is [date1, date2[, including date1 but excluding date2]
def get_himawari(date1,date2,gg):
    current_date = date1
    while current_date < date2:
        try:
            ddH = SAFNWCnc.SAFNWC_CTTH(current_date,'himawari',BBname='SAFBox')
            ddH_ct = SAFNWCnc.SAFNWC_CT(current_date,'himawari',BBname='SAFBox')
            ddH._CTTH_PRESS()
            ddH_ct._CT()
            ddH.var['CT'] = ddH_ct.var['CT']
            ddH._get_var('ctth_status_flag')
            ddH._get_var('ctth_quality')
            ddH_ct._get_var('ct_multilayer')
            ddH.var['ct_multilayer'] = ddH_ct.var['ct_multilayer']
            phH=geosat.SatGrid(ddH,gg)
            phH._sat_togrid('CTTH_PRESS')
            phH._sat_togrid('CT')
            phH._sat_togrid('ctth_status_flag')
            phH._sat_togrid('ctth_quality')
            phH._sat_togrid('ct_multilayer')
            # If mask done here, tend to be persistent (side effect of yield?)
            #phH.var['CTTH_PRESS'][phH.var['CTTH_PRESS']<0]=np.ma.masked
            ddH.close()
            ddH_ct.close()
            del ddH
            del ddH_ct
            phH.expiration = current_date + timedelta(minutes=delta_himawari)
            yield phH
            del phH
        except AttributeError:
            print('Himawari Attribute Error')
            yield -1
        except:
            print('MISSING himawari',current_date)
            yield None
        finally:
            current_date += timedelta(minutes=delta_himawari)
            
if __name__ == '__main__':
    main()
