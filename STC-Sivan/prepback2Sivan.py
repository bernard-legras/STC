#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 09:57:27 2016
Adapted for the Sivan backward runs

This script generates the part_000 file containing parcels on a given potential temperature
for a backward run. The parcels are generated on a one degree centered grid (xcent, ycent) in the
target_range which is also used for the analysis of forward runs from convection.
The parcels are launched every hour which means 8500 x 24 parcels per day in the
standard setting or 18768000 over 92 days (if JAS)
As the first grid is on the date_end (on 0h if no hour given)
and the last grid is on the date_beg (at 0h if no hour given), the total number of parcels is
18768000 + 8500 = 18776500

The output is produced as a part_000 format 107 file which can be used in a STC
backward run. The format 107 is used for convenience even if the flag field is of no use
here. flag is set to the uniform value 0x7F (=127) (Fill grid, New parcel,
Launch time from stamp date, theta level, see format 107 index)

Possible bug: the horizontal interpolation is defined in longitude between 0 and 359
as this is the longitude range in the ERA5 FULL. In interp3d_thet a 2d interpolation
is defined on this grid and then used for longitudes that go from -9.5 to 159., hence
in a extrapolation mode for the range [-9.5 -0.5]. It does not seem to cause many problems
as the parcels in this range are of negligible importance (they mostly escape immediately)
and the extrapolation does not induces a strong variation of theta (which is found to be the
nominal value +- 0.5K). Nevertheless this should be corrected by shifting the data grid in
longitude.

@author: Bernard Legras
"""

from __future__ import division, print_function

import os
import socket
from datetime import datetime, timedelta
#from dateutil.relativedelta import relativedelta
#from monthdelta import MonthDelta
import numpy as np
#import math
#import pygrib
from scipy.interpolate import RectBivariateSpline
import io107
from ECMWF_N import ECMWF
import argparse
import constants as cst
pref = 101325.
# bbox, note it differs from transit.py, y is first
target_range = np.array([[0.,50.],[-10.,160.]])
target_binx = 170; target_biny = 50
deltay = (target_range[0,1]-target_range[0,0])/target_biny
deltax = (target_range[1,1]-target_range[1,0])/target_binx
ycent = np.arange(target_range[0,0] + 0.5*deltay,target_range[0,1],deltay)
xcent = np.arange(target_range[1,0] + 0.5*deltax,target_range[1,1],deltax)

# definition zone for default parameters

#base_year = 2017
#base_month = 7
#theta_level = 380.
#R = 287.04
#Cp = 1005.7
#kappa = R/Cp
#p0 = 100000.

# time step 1 hour
tstep = timedelta(hours=1)
# interval between two ECMWF fields
ENint = timedelta(hours=3)

# %%
def interp3d(data,pt):
    """ Produces an interpolation of the temperature field to a 2d regular grid
    on a given pressure level.The grid provided by (ycent, xcent)
    might be different from that of the analysis.
    Input arguments:
        data : as read from read_ECMWF
        pt   : target pressure in Pa
        xcent : longitudes of the target grid
        ycent : latitudes of the target grid """
    # Interpolate the ECMWF grid onto the target pressure
    # pressure 'P' must be increasing, that is from top to bottom
    Tg = np.empty(shape=[data.nlat,data.nlon])
    for j in range(data.nlat):
        for i in range(data.nlon):
            Tg[j,i] = np.interp(np.log(pt),np.log(data.var['P'][:,j,i]),data.var['T'][:,j,i])
    # Then do a 2d interpolation on this surface to the target grid
    f =  RectBivariateSpline(data.attr['lats'],data.attr['lons'],Tg)
    return f(ycent,xcent)

def interp3d_thet(data,theta):
    """ Produces an interpolation of the temperature field to a 2d regular grid
    on a given potential temperature level. The grid provided by (ycent, xcent)
    might be different from that of the analysis.
    Input arguments:
        data : as read from read_ECMWF
        theta: target potential temperature in K
        ycent : longitudes of the target grid
        xcent : latitudes of the target grid """
    # Interpolate the ECMWF grid onto the target potential temperature
    # potential temperature 'PT' must be decreasing, that is from top to bottom
    Tg = np.empty(shape=[data.nlat,data.nlon])
    Pg = np.empty(shape=[data.nlat,data.nlon])
    for j in range(data.nlat):
        for i in range(data.nlon):
			# - sign because interp wants growing abscissa
            Tg[j,i] = np.interp(-theta,-data.var['PT'][:,j,i],data.var['T'][:,j,i])
            Pg[j,i] = np.exp(np.interp(-theta,-data.var['PT'][:,j,i],np.log(data.var['P'][:,j,i])))
            #f1T = interp1d(-data.var['PT'][:,j,i],data.var['T'][:,j,i],kind='cubic',assume_sorted=False)
            #f1P = interp1d(-data.var['PT'][:,j,i],np.log(data.var['P'][:,j,i]),kind='cubic',assume_sorted=False)
            #Tg[j,i] = f1T(-theta)
            #Pg[j,i] = math.exp(f1P(-theta))
            #Tg[j,i] = krogh_interpolate(-data.var['PT'][:,j,i],data.var['T'][:,j,i],-theta)
            #Pg[j,i] = krogh_interpolate(-data.var['PT'][:,j,i],data.var['P'][:,j,i],-theta)
    # Then do a 2d interpolation on this surface to the target grid
    fT =  RectBivariateSpline(data.attr['lats'],data.attr['lons'],Tg)
    fP =  RectBivariateSpline(data.attr['lats'],data.attr['lons'],Pg)
    return fT(ycent,xcent), fP(ycent,xcent)
# %%
if __name__ == '__main__':
    """ Produce the sequence of initial positions and initial times to be used in
    the backward run. Generates output as part_000 file."""

    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    #parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
    parser.add_argument("-l","--level",type=float,help="theta level")
    parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")

    # Default values
    base_year = 2017
    #base_month = 8
    #base_level = 100. # pressure in hPa
    theta_level = 380.
    quiet = False

    args = parser.parse_args()
    if args.year is not None:
        base_year = args.year
    #if args.month is not None:
    #    base_month = args.month
    if args.level is not None:
        theta_level = args.level
    if args.quiet is not None:
        if args.quiet=='y':
            quiet=True
        else:
            quiet=False

    # Derived parameters defining the run
    # Pressure
    #pt = 100*base_level

    # Dates beginning and end
    # Needs to start at 0h next month to be on an EN date, avoids complicated adjustment
    # for negligible effect.
    date_beg = datetime(year=base_year, month=7, day=1, hour=0)
    date_end = datetime(year=base_year, month=10, day=1, hour=0)
    #date_end = date_beg + timedelta(hours=2)
    # Test on gort
    if 'gort' in socket.gethostname():
        date_beg = datetime(year=2005,month=1,day=1,hour=15)
        date_end = datetime(year=2005,month=1,day=1,hour=21)

    # Define main directories
    if 'ciclad' in socket.gethostname():
            flexout = '/data/legras/flexout/STC/Sivan'
            #part_dir = os.path.join(flexout,'BACK-EAZ-'+date_beg.strftime('%b')\
            #                        +'-'+str(int(base_level))+'hPa')
            part_dir = os.path.join(flexout,'Sivan-EAD-'+date_beg.strftime('JAS-%Y-')\
                                   +str(int(theta_level))+'K')
            try:
                os.mkdir(part_dir)
            except:
                pass

    else:
        print('This script must be run on ciclad')

    # Generate the grid of points
    xg = np.tile(xcent,(target_biny,1))
    yg = np.tile(ycent,(target_binx,1)).T
    bloc_size = target_binx * target_biny
    xg = np.reshape(xg,bloc_size)
    yg = np.reshape(yg,bloc_size)

    # Generate the dictionary to be used to write part_000
    part0 = {}
    # Heading data
    part0['lhead'] = 3
    part0['outnfmt'] = 107
    part0['mode'] = 3   # modify that
    part0['stamp_date'] = date_end.year*10**10 + date_end.month*10**8 + \
        date_end.day*10**6 + date_end.hour*10**4 + date_end.minute*100
    part0['itime'] = 0
    part0['step'] = 450
    part0['idx_orgn'] = 1
    part0['nact_lastO'] = 0
    part0['nact_lastNM'] = 0
    part0['nact_lastNH'] = 0
    part0['flag'] = np.empty(0,dtype=int)
    part0['ir_start'] = np.empty(0,dtype=int)
    part0['x'] = np.empty(0,dtype=float)
    part0['y'] = np.empty(0,dtype=float)
    part0['t'] = np.empty(0,dtype=float)
    part0['p'] = np.empty(0,dtype=float)
    part0['idx_back'] = np.empty(0,dtype=int)
    # It is assumed that both date_end and date_beg are valid dates for the EN files
    date_f = date_end
    date_p = date_f - ENint
    date_current =  date_end

    # first read and first interpolation
    #data = read_ECMWF(ENdir,date_f)
    # data is read from ERA5 data in the global domain at 1° horizontal resolution
    # and longitude origin at 0°
    data = ECMWF('FULL-EA',date_f)
    data._get_var('T')
    data._mkp()   # Calculate pressure from surface pressure and hybrid coefficients
    data._mkthet()  # Calculate the potential temperature
    # ACHTUNG: TO DO : ADD HERE A SHIFT IN LONGITUDE
    # data -> data.shit2west(-179)
    if not quiet:
        print('load first EN file ',date_f)
    T_p, P_p = interp3d_thet(data,theta_level) # First interpolation to the theta level
    delta_pf = ENint.total_seconds()
    numpart = 0

    # loop on the EN time
    while date_p >= date_beg:
        T_f = T_p.copy()
        P_f = P_p.copy()
        if not quiet:
            print('load new EN file  ',date_p)
        data = ECMWF('FULL-EA',date_p)
        data._get_var('T')
        data._mkp()
        data._mkthet()
        T_p, P_p = interp3d_thet(data,theta_level)
        data.close()

        while date_current >= date_p:
            # interpol the temperature in time
            #print('date  ',date_current)
            delta_f = date_f - date_current
            delta_p = date_current - date_p
            T_current = (delta_p.total_seconds()/delta_pf)*T_f + (delta_f.total_seconds()/delta_pf)*T_p
            P_current = (delta_p.total_seconds()/delta_pf)*P_f + (delta_f.total_seconds()/delta_pf)*P_p
            # fill part0
            ir_start = - int((date_end - date_current).total_seconds())
            idx1 = numpart
            numpart += bloc_size
            part0['x'] = np.append(part0['x'],xg)
            part0['y'] = np.append(part0['y'],yg)
            part0['t'] = np.append(part0['t'],T_current)
            part0['p'] = np.append(part0['p'],P_current)
            part0['ir_start'] = np.append(part0['ir_start'],np.full(bloc_size,ir_start,dtype=int))
            part0['idx_back'] = np.append(part0['idx_back'],np.arange(idx1+1,numpart+1,dtype=int))
            part0['flag'] = np.append(part0['flag'],np.full(bloc_size,127,dtype=int))
            # theta check
            theta_check = T_current*(P_current/cst.p0)**(-cst.kappa)
            if not quiet:
                print(theta_check.min()-theta_level,theta_check.max()-theta_level)

            # increment of the date
            date_current = date_current - tstep

        date_f = date_p
        date_p = date_f - ENint

    # final size information
    part0['numpart'] = numpart
    part0['nact'] = numpart

    # write the result as part_000 file
    if not os.path.exists(part_dir):
        os.makedirs(part_dir)
    newpart0=os.path.join(part_dir,'part_000')
    io107.writeidx107(newpart0,part0)
