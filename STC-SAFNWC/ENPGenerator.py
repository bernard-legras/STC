#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code generates the pressure ENP files required by the SAF code at high
resolution in the vertical distribution of pressure levels.

It uses the ERA5 1°x1°, 3-hourly archive on ciclad (requiring EN, VOZ, QN).

It requires additional surface data on the output grid.

What is does:
    (1) It reads the ERA5 data and interpolates them on a set of 68 pressure levels
    on a 1°x1° grid with the FULLAMA domain (0-50N, 10W-160E). Interpolation is done
    for T, Q, O3, Z
    (2) It calculates the WMO tropopause (pressure, temperature, altitude)
    (3) It generates a combined grib file using template messages for pressure levels and
    for surface data. This is done in sequence by (a) copying the surface data, (b) inserting the
    tropopause information, (c) inserting the pressure levels from top to bottom.

Notes
- The latitude order of the interpolated data and the tropopause data are reversed before writing to
the grib file. This is because the expected order from an ECMWF grib file is from north to south and
the latitudes are ordered in that way in the metadata, while the ECMWF_n package reverses this order
on input.
- The eccodes package restablishes consistency among the metadata when a parameter is modified in the
metadata or when the data is modified. There are a number of parameters which are read-only and
cannot be modified. There are however several ways to modify the date for instance. As this is not
documented and we do not know what is under the hood, it is recommended to check thoroughly the
result of surgery on grib messages.
- The ECMWF parameter ids for tropopause data are not suitable as they are on 6 digits and the SAFNWC
code does not handle that since the id is coded in a too short integer. Therefore we use 81, 82 and 83
which are still unused. Fortunately, the SAFNWC-GEO code checks only the paramId.
- The altitudes are in geopotential units, that is m**2 s**-1.
- The surface pressure is expected to be in Pascal but the tropopause pressure is expected to be in
hPa. This is at least what says the SAFNWC-GEO doc.
- The 68 retained levels are made from using the basic 37 levels pressure levels of ECMWF, reduced to
32 below 10 hPa, and to replace the levels in the range 50-650 hPa by the standard pressure levels
of the full 137 hybrid grid. These levels are rounded to closest integers as the grib messages
only support integer values for levels. The levels are in hPa.
- The vertical interpolation is wrong below the ground. This is also true in the standard ECMWF products
on pressure levels and does not seem to matter. These grid points under the ground are not flagged
by the missing value. Using such a missing value would be possibly damageable for the compression method.
- As the volume of data produced at low horizontal resolution is small, we use simple grid packing.
- The grid edition is set to 1 aq this is how the ERA5 pressure data are produced and possibly the
SAFNWC-GEO code does not read edition 2 messages.

Created on Fri Apr  9 22:03:44 2021

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
from os.path import join
from subprocess import call
import eccodes as ec
import argparse
import socket

from ECMWF_N import ECMWF
import constants as cst

# The levels
# The pressure levels of the ERA5 archive
ERA5_PL=(1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400,
350, 300, 250, 225, 200, 175, 150, 125, 100, 70, 50, 30, 20, 10)

# The ARPEGE pressure levels
ARPEGE_levs=(10 ,20 ,30 ,50 ,70 ,100 ,150 ,200 ,250, 300, 400, 500, 600, 700, 800, 850, 900, 925,
             950, 1000)

# The standard pressures on the ERA5 model levels
ERA5_ML = (0.0100018, 0.025513, 0.0388416, 0.0574703, 0.0828747, 0.116762, \
    0.161072, 0.217973, 0.289857, 0.379325, 0.489174, 0.62238, 0.782082, \
    0.971558, 1.19421, 1.45352, 1.75309, 2.09654, 2.48754, 2.9298, \
    3.42702, 3.98287, 4.60104, 5.28515, 6.03877, 6.86542, 7.76856, \
    8.75156, 9.81773, 10.9703, 12.2123, 13.5469, 14.977, 16.5054, \
    18.1348, 19.8681, 21.7076, 23.656, 25.7156, 27.8887, 30.1776, \
    32.5843, 35.1111, 37.7598, 40.5321, 43.4287, 46.4498, 49.5952, \
    52.8644, 56.2567, 59.7721, 63.4151, 67.1941, 71.1187, 75.2, 79.4495, \
    83.8814, 88.5113, 93.353, 98.4164, 103.71, 109.242, 115.02, 121.053, \
    127.349, 133.917, 140.766, 147.906, 155.344, 163.092, 171.159, \
    179.554, 188.287, 197.368, 206.808, 216.616, 226.805, 237.384, \
    248.363, 259.755, 271.57, 283.82, 296.516, 309.669, 323.291, 337.393, \
    351.989, 367.089, 382.705, 398.851, 415.539, 432.779, 450.586, \
    468.971, 487.947, 507.502, 527.569, 548.031, 568.768, 589.68, \
    610.664, 631.619, 652.442, 673.035, 693.304, 713.163, 732.533, \
    751.343, 769.533, 787.053, 803.862, 819.931, 835.236, 849.767, \
    863.519, 876.496, 888.706, 900.167, 910.897, 920.919, 930.262, \
    938.953, 947.024, 954.506, 961.431, 967.832, 973.739, 979.186, 984.2, \
    988.813, 993.053, 996.945, 1000.52, 1003.79, 1006.79, 1009.54, \
    1012.05)

# The 68 retained levels
pressures = [10 ,20 ,30 ,50, 56, 63, 71, 75, 79, 83, 88, 93, 98, 103, 109, 115, 121, \
    127, 134, 141, 148, 155, 163, 171, 179, 188, 197, 207, 217, 227, 237, \
    248, 260, 272, 284, 297, 310, 323, 337, 352, 367, 383, 399, 416, 433, 451, \
    469, 488, 508, 528, 548, 569, 590, 611, 632, 660, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]

pressPa = [100*x for x in pressures]

# Constants used to compute water vapour saturated pressure with respect to liquid water
RALPW = .6022274788e+02
RBETW = .6822400210e+04
RGAMW = .5139266694e+01
# Constants used to compute water vapour saturated pressure with respect to ice
RALPS = .3262117981e+02
RGETS = .6295421339e+04
RGAMS = .5631331575e+00

# box: domain over which the run is made
box = "ACCLIP"

def main():

    # Setting the paths according to where it runs (gort or ciclad)
    if 'gort' == socket.gethostname():
        data_dir = '/dkol/data'
        FullERA5Dir = '/dsk2/ERA5'
    elif 'spirit' in socket.gethostname():
        FullERA5Dir = '/home/legras/ERA5'
        if box == 'FullAMA': data_dir = '/data/legras/flexpart_in/STC/ERA5'
        elif box == 'ACCLIP': data_dir = '/data/legras/flexpart_in/ACCLIP/ERA5'
    else:
        print('unknown hostname for this program')

    # Setting the dimensions
    if box == 'FullAMA':
        NX = 171
        NY = 51
        lonmin = -10
        lonmax = 160
        latmin = 0
        latmax = 50
    elif box == 'ACCLIP':
        NX = 251
        NY = 61
        lonmin = -10
        lonmax = 240
        latmin = 0
        latmax = 60

    SAFNWP_dir = join(data_dir,'SAFNWP')
    MAINOUT_dir = join(SAFNWP_dir,'HVR-LHR')
    SurfData_dir = join(FullERA5Dir,'SURF')

    # Parsing the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-y","--year",type=int,help="year")
    parser.add_argument("-m1","--month1",type=int,choices=1+np.arange(12),help="start month")
    parser.add_argument("-d1","--day1",type=int,choices=1+np.arange(31),help="start day")
    parser.add_argument("-m2","--month2",type=int,choices=1+np.arange(12),help="end month")
    parser.add_argument("-d2","--day2",type=int,choices=1+np.arange(31),help="end day")

    year = 2017
    month1 = 8
    day1 = 23
    month2 = 8
    day2 = 24

    print('parsing arguments')
    args = parser.parse_args()
    if args.year is not None: year = args.year
    if args.month1 is not None: month1 = args.month1
    if args.month2 is not None: month2 = args.month2
    if args.day1 is not None: day1 = args.day1
    if args.day2 is not None: day2 = args.day2

    # To be a loop on time
    date1 = datetime(year,month1,day1,0)
    date2 = datetime(year,month2,day2,0)

    # Getting the templates
    gid = {}
    # Surface template
    ftpt = open('SurfData_template','rb')
    gid['surf'] = ec.codes_grib_new_from_file(ftpt)
    ftpt.close()
    # Column templates
    ftpt = open('ENP_template','rb')
    for var in ['Z','T','RH','O3']:
        gid[var] = ec.codes_grib_new_from_file(ftpt)
    ftpt.close()

    # Time loop
    date  = date1
    while date <= date2:
        # Building the interpolated data
        dat4 = Pinterpol(date,box=box)
        # Transform the geopotential into geopotential altitude
        try: dat4.var['Z'] *= cst.g
        except: pass
        # Calculate the relative humidity
        dat4.var['RH'] = 100 * np.array(pressPa)[:,None,None] * eq(dat4.var['Q'])/ew(dat4.var['T'])

    #%%
        # Defining output file
        file_out = join(MAINOUT_dir,date.strftime('%Y/%m'),date.strftime('ENP%y%m%d%H'))
        # Defining file hosting the surface data
        file_surf = join(SurfData_dir,date.strftime('%Y'),date.strftime('LNSP%y%m%d'))
        file_surf2 = join(SurfData_dir,date.strftime('%Y'),date.strftime('LNSP%Y%m%d.grb'))

        # Copying surface data to the output file
        call(['grib_copy','-w','hour='+str(date.hour),file_surf,file_out])
        #try: call(['grib_copy','-w','hour='+str(date.hour),file_surf,file_out])
        #except: call(['grib_copy','-w','hour='+str(date.hour),file_surf2,file_out])

        # Open the output file in append mode
        fout = open(file_out,'ab')

        # Add first the tropopause data
        dicwmo = {'pwmo':[82,'WMO tropopause pressure','hPa'],
                  'Twmo':[81,'WMO tropopause temperature','T'],
                  'zwmo':[83,'WMO tropopause altitude','m**2 s**-2']}
        dat4.d2d['pwmo'] /= 100
        dat4.d2d['zwmo'] *= cst.g
        for var in ['Twmo','pwmo','zwmo']:
            clone_id = ec.codes_clone(gid['surf'])
            #nx = ec.codes_get(gid['surf'], 'Ni')
            #ny = ec.codes_get(gid['surf'], 'Nj')
            ec.codes_set(clone_id,'Ni',NX)
            ec.codes_set(clone_id,'Nj',NY)
            ec.codes_set(clone_id,'paramId',dicwmo[var][0])
            ec.codes_set(clone_id,'dataDate',10000*date.year+100*date.month+date.day)
            ec.codes_set(clone_id,'hour',date.hour)
            ec.codes_set(clone_id,'packingType','grid_second_order')
            ec.codes_set_values(clone_id,np.reshape(dat4.d2d[var][::-1,:],NX*NY))
            ec.codes_write(clone_id,fout)
            ec.codes_release(clone_id)

        # Add now the data on pressure levels
        for ll in range(len(pressures)):
            for var in ['Z','T','RH','O3']:
                clone_id = ec.codes_clone(gid[var])
                #nx = ec.codes_get(gid[var], 'Ni')
                #ny = ec.codes_get(gid[var], 'Nj')
                ec.codes_set(clone_id,'Ni',NX)
                ec.codes_set(clone_id,'Nj',NY)
                ec.codes_set(clone_id,'lev',pressures[ll])
                ec.codes_set(clone_id,'dataDate',10000*date.year+100*date.month+date.day)
                ec.codes_set(clone_id,'hour',date.hour)
                ec.codes_set(clone_id,'packingType','grid_second_order')
                ec.codes_set_values(clone_id,np.reshape(dat4.var[var][ll,::-1,:],NX*NY))
                ec.codes_write(clone_id,fout)
                ec.codes_release(clone_id)

        # Closing the output file
        fout.close()
        print('processed ',date)
        date += timedelta(hours=3)

def Pinterpol(date,box='FullAMA'):
    """ Read the ERA5 data and interpolate to the pressure levels. """
    if box == 'FullAMA':
        lonmin = -10
        lonmax = 160
        latmin = 0
        latmax = 50
    elif box == 'ACCLIP':
        lonmin = -10
        lonmax = 240
        latmin = 0
        latmax = 60
    dat = ECMWF('FULL-EA',date,exp=['VOZ','QN'])
    dat._get_T()
    dat._get_var('O3')
    dat._get_var('Q')
    dat.close()
    dat._mkp()
    dat._mkz()
    dat2 = dat.shift2west(-20)
    dat3 = dat2.extract(lonRange=[lonmin,lonmax],latRange=[latmin,latmax],varss='All')
    dat3._WMO()
    dat4 = dat3.interpolP(pressPa,varList=['Z','T','Q','O3'])
    dat4.d2d = dat3.d2d
    return(dat4)

eq = lambda q: q/(q + (cst.R/cst.Rv)*(1-q))

def ew(T):
    if T > cst.Zero_Celsius:
        ew = np.exp(RALPW-RBETW/T-RGAMW*np.log(T))
    else:
        ew = np.exp(RALPS-RGETS/T-RGAMS*np.log(T))
    return ew
ew = np.vectorize(ew)

if __name__ == '__main__':
    main()
