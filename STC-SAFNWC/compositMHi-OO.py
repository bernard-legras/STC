#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Makes a composite from Himawari and msg1 SAFNWC cloud tops

ICARE version, old MSG1 + old Himawari

Created on Tue Mar 20 00:34:26 2018

@author: Bernard Legras
"""
#from SAFNWC import SAFNWC_CTTH
from SAFNWCi import SAFNWC_CTTH
import geosati
from datetime import datetime, timedelta
import pickle,gzip
import os
import numpy as np


def main():
    date1 = datetime(2017,4,1,0)
    date2 = datetime(2017,5,1,0)
    #
    
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
                phM.expiration += timedelta(minutes=15)
        if phH.expiration == date:
            pp = next(gH)
            # If available image, get it
            # if not,  extend validity of previous image
            if pp is not None:
                phH = pp
                make_new_comp = True
            else:
                phH.expiration += timedelta(minutes=20)
        
        if make_new_comp:
            comb = phM.patch(phH,90.75,'CTTH_PRESS')
            comb.var['CTTH_PRESS'][comb.var['CTTH_PRESS']<0]=np.ma.masked
            comb.var['CTTH_PRESS'][comb.var['CTTH_PRESS']>65000]=np.ma.masked
            filename = date.strftime('SAFNWC-PTOP-%Y-%m-%d-%H:%M.pkl')
            fullname = os.path.join(output_dir,date.strftime('%Y/%m'),filename)
            with gzip.open(fullname,'wb') as f:
                pickle.dump(comb.var['CTTH_PRESS'],f,protocol=3)
            print ('new file',date,np.ma.count_masked(comb.var['CTTH_PRESS']))
            del comb
                
        date += timedelta(minutes=5)
        make_new_comp = False
    
# Generator providing the msg1 data between two dates 
# The interval is [date1, date2[, including date1 but excluding date2]
def get_msg1(date1,date2,gg):
    current_date = date1
    while current_date < date2:
        try:
            ddM=SAFNWC_CTTH(current_date,'msg1')
            ddM._CTTH_PRESS()            
            phM=geosati.SatGrid(ddM,gg)
            phM._sat_togrid('CTTH_PRESS')
            # If mask done here, tend to be persistent (side effect of yield?)
            #phM.var['CTTH_PRESS'][phM.var['CTTH_PRESS']>65000.]=np.ma.masked
            ddM.close()
            del ddM
            phM.expiration = current_date + timedelta(minutes=15)
            yield phM
            del phM
        except AttributeError:
            yield None
        except:
            print('ERROR ERROR in get_msg1')
            yield -1
        finally:
            current_date += timedelta(minutes=15)
            
# Generator providing the himawari data between two dates 
# The interval is [date1, date2[, including date1 but excluding date2]
def get_himawari(date1,date2,gg):
    current_date = date1
    while current_date < date2:
        try:
            ddH=SAFNWC_CTTH(current_date,'himawari')
            ddH._CTTH_PRESS()            
            phH=geosati.SatGrid(ddH,gg)
            phH._sat_togrid('CTTH_PRESS')
            # If mask done here, tend to be persistent (side effect of yield?)
            #phH.var['CTTH_PRESS'][phH.var['CTTH_PRESS']<0]=np.ma.masked
            ddH.close()
            del ddH
            phH.expiration = current_date + timedelta(minutes=20)
            yield phH
        except AttributeError:
            yield None
        except:
            print('ERROR ERROR in get_himawari')
            yield -1
        finally:
            current_date += timedelta(minutes=20)
            
if __name__ == '__main__':
    main()
