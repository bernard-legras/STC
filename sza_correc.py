# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 02:01:06 2016

From the Fortran code of Bob Joyce
after Joyce et al., JAMT, 2001
doi:10.1175/1520-0450(2001)040%3C0689:LASDZA%3E2.0.CO;2

This implement the zenith angle correction

According to Joyce et al., parallax correction should be done before.
Gridsat is processed by doing parallax correction after.
Does it really matter?
Anyway this code does not include parallax correction at the moment.
To be included in a next version. 

@author: Bernard Legras
"""

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from struct import unpack
import numpy as np
import os
import socket
from datetime import datetime, timedelta

# input dir where the Joyce coefficients are stored
# todo : repliace by environment variables
if socket.gethostname() == 'Graphium':
    INPUT_dir = 'C:\\cygwin64\\home\\berna\\data\\STC\\pylib\\INPUT'
elif 'ciclad' in socket.gethostname():
    INPUT_dir = '/home/legras/TRAJ/pylib'
elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
    INPUT_dir = '/home/stratocl/TRAJ/pylib/INPUT'
elif socket.gethostname() == 'grapelli':
    INPUT_dir = '/limbo/data/STC/pylib/INPUT'
elif socket.gethostname() == 'zappa':
    INPUT_dir = '/net/grapelli/limbo/data/STC/pylib/INPUT'
elif socket.gethostname() == 'gort':
    INPUT_dir = '/dkol/data/STC/pylib/INPUT'
elif 'icare' in socket.gethostname():
    INPUT_dir = '/home/b.legras/pylib/INPUT'
else:
     print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')

za_bin = 64
bt_bin = 170
# number of bins in latitude
lat_bin = 240
# min and max corrected latitudes
lat_min = -60.
lat_max = 60.
# latitude interval
lat_int = 0.5
# number of seasonal values
n_seas = 8
# maximum admitted zenith angle
sza_max = 87.0
# minimum zenith angle for correction
sza_angmin = 26.
# minimum brigthnesss temperature
btmin = 160.
# brightness temperature range (asssumed step =1K)
btrange = 170
# length of the year
year_l = 365
# Earth radius
R = 6371.
# altitude geostationary orbit
H = 35680.
# missing value
MISSING = -9999
FillValue = -31999

def sza_model():
    """ Read the Joyce coefficients from source files.
    Quite interestingly, error_regrid is stored as little-endian floats
    while latitudinal-error_annual-cycle-model is stored as big-endian floats.
    Just for fun???"""
    model = {}
    # read source files
    filename = os.path.join(INPUT_dir,'error_regrid')
    fid = open(filename,'rb')
    model['error'] = np.reshape(np.asarray(unpack('<'+str(bt_bin*za_bin)+'f',
                       fid.read(4*bt_bin*za_bin))), [bt_bin, za_bin])
    fid.close()
    filename = os.path.join(INPUT_dir,'latitudinal-error_annual-cycle-model')
    fid = open(filename,'rb')
    model['errorlat'] = np.zeros(shape=[lat_bin, n_seas])
    model['difflat']  = np.zeros(shape=[lat_bin, n_seas])
    for seas in range(n_seas):
        model['errorlat'][:,seas] =  np.asarray(unpack('>'+str(lat_bin)+'f',
                       fid.read(4*lat_bin)))                
        model['difflat'][:,seas] =  np.asarray(unpack('>'+str(lat_bin)+'f',
                       fid.read(4*lat_bin)))
    fid.close()
    # discrete times from the Fortran file - 1 to account that January 1 at
    # OUTC is 0 (not 1 like in idl or fortran)
    model['seas_day'] = np.array([57.5,104.5,148.5,206.5,268.5,326.5,395,422.5])-1
    return model

def zenith_angle(lon, lat, sublon, sublat):
    """ Calculate the zenith angle from vectors of lat and lon and the position
    of the satellite. All inputs and output in degrees. This code works whatever
    is the shape of the arguments (1D or 2D, masked array or not)"""
    # Calculate the arc angle between satellite and target point
    # Unlike Joyce code, this does not assume that the satellite is on the equator
    # (as it is not true for, e.g., Meteosat 7 which wanders about 10°
    # north and south of the equator)
    cos_psi = np.sin(np.radians(lat))*np.sin(np.radians(sublat)) \
            + np.cos(np.radians(lat))*np.cos(np.radians(sublat)) \
            *np.cos(np.radians(lon-sublon))
    sin_psi=np.sqrt(1-cos_psi**2)
    # Calculate the distance of satelitte to ground target point
    L = np.sqrt((R+H)**2 + R**2 - 2*R*(R+H)*cos_psi)
    del cos_psi
    # Calculate the zenith angle
    zeta = np.degrees(np.arcsin(np.clip(((R+H)/L)*sin_psi,0.,1.)))
    return zeta


def sel_errorlat(model, date):
    """ Calculate the vector of seasonally depend latitude correction for the
    actual date.
    Recall that the Joyce model starts on 26 Feb at 12 UTC and does not include bissextil
    years"""
    jday = (date-datetime(year=date.year,month=1,day=1)).total_seconds()/86400.
    if (jday >= 59.5) & (date.year%4==0):
        jday -= 1
    if jday < model['seas_day'][0]:
        jday += year_l
    jp = np.where(model['seas_day']>=jday)[0][0]  
    wp = (jday-model['seas_day'][jp-1])/(model['seas_day'][jp]-model['seas_day'][jp-1])
    wm = (model['seas_day'][jp]-jday)/(model['seas_day'][jp]-model['seas_day'][jp-1])
    xerrorlat = wm*model['errorlat'][:,jp-1]+wp*model['errorlat'][:,jp]
    #print('jday '+str(jday)+' jp '+str(jp)+' wm '+str(wm)+' wp '+str(wp))
    return xerrorlat


def szacorr(date, TB, lon, lat, sublon, sublat,freeze=False):
    """ Apply the Joyce correction algorithm and generate the correction to be
    added to TB.
    In addition, the input TB is set to missing when the zenit angle is larger
    than its max value. 
    This correction code should work whether TB, lon and lat are 1D or 2D 
    fields and whether they are masked array or not.
    On output TB is masked or set to MISSING for zenith angles larger that 87°, 
    and the correction in K is returned in vza. vza is a masked array with the 
    same mask as TB if inputs are masked arrays. 
   
    """
    # test that the shape of TB and lat, lon is the same
    if (TB.shape != lon.shape) | (TB.shape != lat.shape):
        raise ValueError('shapes of TB, lon lat not identical')
        return -1
    # load model
    model = sza_model()
    # calculate zenith angle
    zeta = zenith_angle(lon, lat, sublon, sublat)
    # interpolate errorlat in time
    xerrorlat = sel_errorlat(model,date)
    # calculate latitude bin
    j5x = np.clip((np.floor((lat_max-lat)/lat_int)).astype(int),0,lat_bin-1)
    # calculate zeta bins
    iz = np.clip((np.floor(zeta-sza_angmin)).astype(int),None,za_bin-1)
    # calculate T bin
    it = np.clip((np.floor(TB-btmin)).astype(int),0,btrange-1)
    vza = np.zeros(shape=it.shape)
    vza[iz>=0] = - xerrorlat[j5x[iz>=0]] * model['error'][it[iz>=0],iz[iz>=0]]
    #print('j5x '+str(j5x[iz>=0][0])+' iz '+str(iz[0])+' it '+str(it[iz>=0][0]))
    #aa =  xerrorlat[j5x[iz>=0]]
    #bb =  model['error'][it[iz>=0],iz[iz>=0]]
    #print('errorlat '+str(aa[0])+' error '+str(bb[0]))
    
    # set to missing the brightness temperatures with zenith angles larger
    # than sza_max
    # test for masked array in input
    try:
        test = np.sum(TB.mask)
    except AttributeError:
        # action for non masked array
        TB[zeta>sza_max] = MISSING
    else:
        # action for masked array
        # 1st: follow non masked array procedure for TB
        # should not happen but necessary if precalculated nearest neighbour
        # lookup table to be applied at next step on compressed() fields 
        # as next step
        # mask should be updated after interpolation
        # todo: find a clean way to handle that case
        if freeze:
            TB[zeta>sza_max] = MISSING
            vza = np.ma.array(vza, mask=TB.mask, fill_value=FillValue)
        # 2nd: normal behaviour for masked arrays     
        else:   
            TB.mask[zeta>sza_max] = True
            vza = np.ma.array(vza, mask=TB.mask, fill_value=FillValue)
            vza.mask[zeta>sza_max] = True
    return vza, zeta
  
if __name__ == '__main__':
    """ Tests for comparison with the Fortran version """
    date0=datetime(year=2013,month=1,day=1,hour=0)
    TB = np.array([230.])
    lat = np.array([59.9])
    lon = np.array([115.])
    sublat = 0.
    sublon = 140.
    #for i in range(366):
    for i in 180+np.arange(10):
        date = date0+timedelta(days=int(i))
        vza, zeta = szacorr(date,TB,lon,lat,sublon,sublat)
        TBc = TB+vza
        print('day '+str(i+1)+' zeta '+str(zeta[0])+' TB '+str(TBc[0]))
        
