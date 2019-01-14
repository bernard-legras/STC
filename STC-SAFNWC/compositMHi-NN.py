#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Makes a composite from Himawari and msg1 SAFNWC cloud tops
in the FULLAMMA domain with an image produced each time
a new image is available from Himawari or msg1.
If an image is misssing, the validity of the previous image 
is extended.
The output contains also the status flag. A mask is applied 
for non valid data from the value of ctth only.

This script process the standard SAF product with 2016 version
algo which is available for Himawari after 11 July 2017.

ICARE version

Created on Tue Mar 20 00:34:26 2018

@author: Bernard Legras
"""
#from SAFNWC import SAFNWC_CTTH
from SAFNWCnci import SAFNWC_CTTH
import geosati
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
    parser.add_argument("-s","--segment",type=str,choices=["Jul.11","Jul.21","Aug.1","Aug.11","Aug.21"],help="segment to be processed")
    seg_dates = {'Jul.11':[datetime(2017,7,11,0),datetime(2017,7,21,0)],
                 'Jul.21':[datetime(2017,7,21,0),datetime(2017,8,1,0)],
                 'Aug.01':[datetime(2017,8,1,0),datetime(2017,8,11,0)],
                 'Aug.11':[datetime(2017,8,11,0),datetime(2017,8,21,0)],
                 'Aug.21':[datetime(2017,8,21,0),datetime(2017,9,1,0)]}
    
    seg = 'Jul.11'

    args = parser.parse_args()
    if args.segment is not None:
        seg = args.segment
    
    date1 = seg_dates[seg][0]
    date2 = seg_dates[seg][1]
  
    gg=geosati.GeoGrid('FullAMA')

    date = date1

    output_dir = '/scratch/b.legras/STC/STC-SAFNWC-OUT'
    
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
            comb = phM.patch(phH,90.75,['CTTH_PRESS','ctth_status_flag','ctth_quality','ctth_conditions'])
            comb.var['CTTH_PRESS'][comb.var['CTTH_PRESS']<0]=np.ma.masked
            comb.var['CTTH_PRESS'][comb.var['CTTH_PRESS']>65000]=np.ma.masked
            # flag the good and the opaq clouds into a single 8bit flag (bit 0 for good, bit 1 for opaq)
            opaq = (comb.var['ctth_status_flag'] & 4) == 4
            good = ((comb.var['ctth_quality']&0x38)==0x8) & ((comb.var['ctth_conditions']&0xFF00)==0x5500)
            comb.var['cloud_flag'] = (good + 2*opaq).astype(np.int8)
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
            ddM=SAFNWC_CTTH(current_date,'msg1')
            ddM._CTTH_PRESS()
            ddM._get_var('ctth_status_flag')
            ddM._get_var('ctth_quality')
            ddM._get_var('ctth_conditions')
            phM=geosati.SatGrid(ddM,gg)
            phM._sat_togrid('CTTH_PRESS')
            phM._sat_togrid('ctth_status_flag')
            phM._sat_togrid('ctth_quality')
            phM._sat_togrid('ctth_conditions')
            # If mask done here, tend to be persistent (side effect of yield?)
            #phM.var['CTTH_PRESS'][phM.var['CTTH_PRESS']>65000.]=np.ma.masked
            ddM.close()
            del ddM
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
            ddH=SAFNWC_CTTH(current_date,'himawari')
            ddH._CTTH_PRESS()
            ddH._get_var('ctth_status_flag')
            ddH._get_var('ctth_quality')
            ddH._get_var('ctth_conditions')
            phH=geosati.SatGrid(ddH,gg)
            phH._sat_togrid('CTTH_PRESS')
            phH._sat_togrid('ctth_status_flag')
            phH._sat_togrid('ctth_quality')
            phH._sat_togrid('ctth_conditions')
            # If mask done here, tend to be persistent (side effect of yield?)
            #phH.var['CTTH_PRESS'][phH.var['CTTH_PRESS']<0]=np.ma.masked
            ddH.close()
            del ddH
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
