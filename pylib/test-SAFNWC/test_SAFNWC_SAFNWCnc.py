#!/usr/bin/env p(ython
# encoding: utf-8

'''

Routine for testing SAFNWC 
Supported products : Cloud Mask (CMa), Cloud Type (CT), Cloud Top Temperature and Height (CTTH)

S. Bucci - Created on Tue May 30 17:21:26 CET 2017

'''
#import tables
#import numpy as np
from datetime import datetime
import numpy as np
import geosat
#import os
#import sza_correc
import matplotlib.colors as mcl
import STC_cmap
import pprint
#import pylab
from SAFNWC import SAFNWC_CT, SAFNWC_CTTH
from SAFNWCnc import SAFNWC_CT_nc, SAFNWC_CTTH_nc 

d=datetime(2017,6,16,1)
#plot cloud Type classificaiton from Himawari on a FullAMA domain
# define satellites
sath='himawari'
satm='msg1'
# define the grid
gg=geosat.GeoGrid('MesoInd')
# Open CT files and read CT
ah= SAFNWC_CT(d,sath)
ah._CT()
ah.close()
am= SAFNWC_CT_nc(d,satm)
am._CT()
am.close()
# Open CTTH files and read PRESS and TEMPER
ath= SAFNWC_CTTH(d,sath)
ath._CTTH_PRESS()
ath._CTTH_TEMPER()
ath.close()
atm= SAFNWC_CTTH_nc(d,'msg1')
atm._CTTH_PRESS()
atm._CTTH_TEMPER()
atm.close()
# Merge CT and CTTH fields
ah._merge(ath)
am._merge(atm)
# Associate geogrid and geosat fields
ph_h=geosat.SatGrid(ah,gg)
ph_m=geosat.SatGrid(am,gg)
# Define bounds of color map
clims={'CT':[0.5,15.5],'CTTH_PRESS':[50,500],'CTTH_TEMPER':[190,300]}
for var_name in ['CT','CTTH_PRESS','CTTH_TEMPER']:
    # Projection onto the grid
    ph_h._sat_togrid(var_name)
    ph_m._sat_togrid(var_name)
    # Define palette
    coltab=ah.attr[var_name]['PALETTE']/255.
    if var_name == 'CT':
        cmap = mcl.ListedColormap(coltab)
    else:
        cmap = STC_cmap.mymap
    # Makes chart for the two satellites separately
    #(deactivated)
    ## ph_h.chart(var_name,cmap=cmap,clim=clims[var_name],txt=var_name+" "+sath)
    ## ph_m.chart(var_name,cmap=cmap,clim=clims[var_name],txt=var_name+" "+satm)
    # Path the two fields at longitude = 90E
    if var_name == 'CT':
        pass
    else:
        print(np.min(ph_m.var[var_name]))
        ph_h._filt(var_name,ath.attr[var_name]['FillValue'])
        print(np.min(ph_m.var[var_name]))
        ph_m._filt(var_name,atm.attr[var_name]['FillValue'])
        print(np.min(ph_m.var[var_name]))
    phm=ph_m.patch(ph_h,90,var_name)
    # Filter the content to eliminate non cloudy pixels
    if var_name == 'CT':
        pass
    #else:
        #phm._filt(var_name,am.attr[var_name]['FillValue'])
    # Plot the chart 
    phm.chart(var_name,cmap=cmap,clim=clims[var_name],txt=var_name+" "+satm+'/'+sath)
    pprint.pprint(ah.attr[var_name]['units'])
