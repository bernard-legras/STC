#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script is meant to manage the data from narrow extractions around stations 
(for the moment, OHP, OPAR and ACARR)
It uses data in netcdf format
Data are read on hybrid levels. These levels are not provided in the file and hence
must be et in the program.
Interpolation is provided onto geopotential or pressure levels.

Functions: Open files, read the data, extract subgrids, make charts, interpolate in time,
interpolate to pressure levels

Usage:
>> from OPZ import ECMWF
open files for a date (datetime) and a project (VOLC or STC)
>> data = ECMWF(project,date)
read a variable var
>> data._get_var(var)
example: data._get_var('CC') for cloud cover (see lists in ECMWF class)
calculate pressure field
>> data._mkp()
calculate potential temperature (requires to read T and calculate P before)
>> data._mktheta()
interpolate in time (date2 must be between the dates of data0 and data1)
>> data2 = data0.interpol-time(data1,date2)
interpolate to pressure levels
>> data1 = data.interpolP(pList,varList,latRange,lonRange)
where pList is a pressure or a list of pressures, varList is a variable or a list of variables,
latRange and lonRange as in extract

The ECMWF_pure allows to define a template object without reading files.
It can be used to modify data. It is produced as an output of the interpolation.

Created on 11/12/2022 from ECMWF_N.py

@author: Bernard Legras (bernard.legras@lmd.ipsl.fr)
@licence: CeCILL-C
"""
import numpy as np
#import math
from netCDF4 import Dataset
import os
#from cartopy import feature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.crs as ccrs
#from cartopy import feature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import socket
from scipy.interpolate import interp1d,RegularGridInterpolator
from mki2d import tohyb
import constants as cst
#import gzip,pickle
from gethyb import hyb
import satratio

#from numba import jit
#from copy import copy,deepcopy

MISSING = -999
# Physical constants
# now in the package "constants"
#R = 287.04 # or 287.053
#Cpd = 1005.7
#kappa = R/Cpd
#g = 9.80665
#pref = 101325.
#p0 = 100000.

listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']
# homemade color table
listcolors_sw=[listcolors[1],listcolors[0],listcolors[3],listcolors[2],\
               listcolors[5],listcolors[4],listcolors[7],listcolors[6],\
               listcolors[9],listcolors[8],listcolors[11],listcolors[10],\
               listcolors[13],listcolors[12],listcolors[15],listcolors[14],\
               listcolors[17],listcolors[16],listcolors[19],listcolors[18]]
mymap=colors.ListedColormap(listcolors)
mymap_sw=colors.ListedColormap(listcolors_sw)


def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

# Second order estimate of the first derivative dy/dx for non uniform spacing of x
d = lambda x,y:(1/(x[2:,:,:]-x[:-2,:,:]))\
                *((y[2:,:,:]-y[1:-1,:,:])*(x[1:-1,:,:]-x[:-2,:,:])/(x[2:,:,:]-x[1:-1,:,:])\
                -(y[:-2,:,:]-y[1:-1,:,:])*(x[2:,:,:]-x[1:-1,:,:])/(x[1:-1,:,:]-x[:-2,:,:]))
    
dvar = {'U':['u','U component of wind','m s**-1'],
         'V':['v','V component of wind','m s**-1'],
         'W':['w','Vertical velocity','Pa s**-1'],
         'T':['t','Temperature','K'],
         'LNSP':['lnsp','Logarithm of surface pressure','Log(Pa)'],
         'VO':['vo','Vorticity','s**-1'],
         'Q':['q','Specific humidity','kg kg**-1'],
         'O3':['o3','Ozone mixing ratio','kg kg**-1']
         }

((am,bm),(al,bl)) = hyb(137)

stations ={'ACARR':[76.332,10.04265],
           'OHP':[5.7133,43.9308],
           'OPAR':[55.3831,-21.0797]}

# template object produced by extraction and interpolation
class ECMWF_pure(object):
    def __init__(self):
        self.var={}
        self.attr={}
        self.warning = []

    def getxy(self,var,lev,y,x):
        """ get the interpolated value of var in x, y on the level lev """
        # Quick n' Dirty version by nearest neighbour
        jy = np.abs(self.attr['lats']-y).argmin()
        ix = np.abs(self.attr['lons']-x).argmin()
        return self.var[var][lev,jy,ix]
    
    def get_station_profile(self,method='bilin'):
        """ This method returns a dictionary of profiles with all 
        variables at the location of the station.
        It is based either on a simple bilinear interpolation or
        on a Shepard IDW method with 12 neighbours. """
        x0 = stations[self.project][0]
        y0 = stations[self.project][1]
        prof = {'project':self.project,'date':self.date}
        if method == 'bilin':
            ix = np.where(self.attr['lons']<x0)[0][-1]
            jy = np.where(self.attr['lats']<y0)[0][-1]
            dx = (x0 - self.attr['lons'][ix])/self.attr['dlo'] 
            dy = (y0 - self.attr['lats'][jy])/self.attr['dla']
            for var in self.var.keys():
                if len(self.var[var].shape) == 2:
                    prof[var] = dx*dy*self.var[var][jy+1,ix+1] \
                           + (1-dx)*dy*self.var[var][jy+1,ix] \
                           + dx*(1-dy)*self.var[var][jy,ix+1] \
                           + (1-dx)*(1-dy)*self.var[var][jy,ix]
                else:
                    prof[var] = dx*dy*self.var[var][:,jy+1,ix+1] \
                           + (1-dx)*dy*self.var[var][:,jy+1,ix] \
                           + dx*(1-dy)*self.var[var][:,jy,ix+1] \
                           + (1-dx)*(1-dy)*self.var[var][:,jy,ix]
        else:
            print('method not implemented')
        return(prof)   
          
    # INNTERPOLATION SECTION

    def interpol_time(self,other,date):
        # This code interpolate in time between two ECMWF objects with same vars
        # check the date
        if self.date < other.date:
            if (date > other.date) | (date < self.date):
                print ('error on date')
                print (self.date,date,other.date)
                return -1
        else:
            if (date < other.date) | (date > self.date):
                print ('error on date')
                print (other.date,date,self.date)
                return -1
        # calculate coefficients
        dt = (other.date-self.date).total_seconds()
        dt1 = (date-self.date).total_seconds()
        dt2 = (other.date-date).total_seconds()
        cf1 = dt2/dt
        cf2 = dt1/dt
        #print ('cf ',cf1,cf2)
        data = ECMWF_pure()
        data.date = date
        for names in ['lons','lats']:
            data.attr[names] = self.attr[names]
        data.nlon = self.nlon
        data.nlat = self.nlat
        data.nlev = self.nlev
        for var in self.var.keys():
            data.var[var] = cf1*self.var[var] + cf2*other.var[var]
        return data

    def interpolP(self,p,varList='All',latRange=None,lonRange=None):
        """ interpolate the variables to a pressure level or a set of pressure levels
            vars must be a list of variables or a single varibale
            p must be a pressure or a list of pressures in Pascal
        """
        if 'P' not in self.var.keys():
            self._mkp()
        new = ECMWF_pure()
        if varList == 'All':
            varList = list(self.var.keys())
            try:
                varList.remove('SP')
                varList.remove('P')
            except: pass # in case these variables are already removed
        elif type(varList) == str:
            varList = [varList,]
        for var in varList:
            if var not in self.var.keys():
                print(var,' not defined')
                return
        if type(p) in [float,int]:
            p = [p,]
        if 'P' not in self.var.keys():
            print('P not defined')
            return
        # first determine the boundaries of the domain
        if (latRange == []) | (latRange == None):
            nlatmin = 0
            nlatmax = self.nlat
        else:
            nlatmin = np.abs(self.attr['lats']-latRange[0]).argmin()
            nlatmax = np.abs(self.attr['lats']-latRange[1]).argmin()+1
        if (lonRange == []) | (lonRange == None):
            nlonmin = 0
            nlonmax = self.nlon
        else:
            nlonmin = np.abs(self.attr['lons']-lonRange[0]).argmin()
            nlonmax = np.abs(self.attr['lons']-lonRange[1]).argmin()+1
        new.attr['lats'] = self.attr['lats'][nlatmin:nlatmax]
        new.attr['lons'] = self.attr['lons'][nlonmin:nlonmax]
        new.attr['Lo1'] = new.attr['lons'][0]
        new.attr['Lo2'] = new.attr['lons'][-1]
        new.attr['La1'] = new.attr['lats'][0]
        new.attr['La2'] = new.attr['lats'][-1]
        new.attr['dlo'] = self.attr['dlo']
        new.attr['dla'] = self.attr['dla']
        new.project = self.project
        new.nlat = len(new.attr['lats'])
        new.nlon = len(new.attr['lons'])
        new.date = self.date
        new.nlev = len(p)
        new.attr['levtype'] = 'pressure'
        new.attr['plev'] = p
        new.attr['levs'] = p
        pmin = np.min(p)
        pmax = np.max(p)
        for var in varList:
            new.var[var] = np.empty(shape=(len(p),nlatmax-nlatmin,nlonmax-nlonmin))
            jyt = 0
            # big loop that should be paralellized or calling numa for good performance
            for jys in range(nlatmin,nlatmax):
                ixt = 0
                for ixs in range(nlonmin,nlonmax):
                    # find the range of p in the column
                    npmin = np.abs(self.var['P'][:,jys,ixs]-pmin).argmin()
                    npmax = np.abs(self.var['P'][:,jys,ixs]-pmax).argmin()+1
                    npmin = max(npmin - 3,0)
                    npmax = min(npmax + 3,self.nlev)
                    # Better version than the linear interpolation but much too slow
                    #fint = PchipInterpolator(np.log(self.var['P'][npmin:npmax,jys,ixs]),
                    #                     self.var[var][npmin:npmax,jys,ixs])
                    #new.var[var][:,jyt,ixt] = fint(np.log(p))
                    new.var[var][:,jyt,ixt] = np.interp(np.log(p),np.log(self.var['P'][npmin:npmax,jys,ixs]),self.var[var][npmin:npmax,jys,ixs])
                    ixt += 1
                jyt += 1
        return new

    def interpolZ(self,z,varList='All',latRange=None,lonRange=None):
        """ interpolate the variables to an altitude level or a set of altitude levels
            vars must be a list of variables or a single varibale
            z must be an altitude or a list of altitudes in m
        """
        if 'Z' not in self.var.keys():
            self._mkz()
        new = ECMWF_pure()
        if varList == 'All':
            varList = list(self.var.keys())
            try:
                varList.remove('SP')
                varList.remove('P')
            except: pass # in case these variables are already removed
        elif type(varList) == str:
            varList = [varList,]
        for var in varList:
            if var not in self.var.keys():
                print(var,' not defined')
                return
        if type(z) in [float,int]:
            z = [z,]
        if 'Z' not in self.var.keys():
            print('Z not defined')
            return
        # first determine the boundaries of the domain
        if (latRange == []) | (latRange == None):
            nlatmin = 0
            nlatmax = self.nlat
        else:
            nlatmin = np.abs(self.attr['lats']-latRange[0]).argmin()
            nlatmax = np.abs(self.attr['lats']-latRange[1]).argmin()+1
        if (lonRange == []) | (lonRange == None):
            nlonmin = 0
            nlonmax = self.nlon
        else:
            nlonmin = np.abs(self.attr['lons']-lonRange[0]).argmin()
            nlonmax = np.abs(self.attr['lons']-lonRange[1]).argmin()+1
        new.attr['lats'] = self.attr['lats'][nlatmin:nlatmax]
        new.attr['lons'] = self.attr['lons'][nlonmin:nlonmax]
        new.nlat = len(new.attr['lats'])
        new.nlon = len(new.attr['lons'])
        new.date = self.date
        new.attr['Lo1'] = new.attr['lons'][0]
        new.attr['Lo2'] = new.attr['lons'][-1]
        new.attr['La1'] = new.attr['lats'][0]
        new.attr['La2'] = new.attr['lats'][-1]
        new.attr['dlo'] = self.attr['dlo']
        new.attr['dla'] = self.attr['dla']
        new.project = self.project     
        new.nlev = len(z)
        new.attr['levtype'] = 'altitude'
        new.attr['levs'] = z
        zmin = np.min(z)
        zmax = np.max(z)
        for var in varList:
            new.var[var] = np.empty(shape=(len(z),nlatmax-nlatmin,nlonmax-nlonmin))
            jyt = 0
            # big loop that should be paralellized or calling numa for good performance
            for jys in range(nlatmin,nlatmax):
                ixt = 0
                for ixs in range(nlonmin,nlonmax):
                    # find the range of z in the column
                    nzmin = np.abs(self.var['Z'][:,jys,ixs]-zmin).argmin()
                    nzmax = np.abs(self.var['Z'][:,jys,ixs]-zmax).argmin()+1
                    nzmin = max(nzmin + 3,self.nlev-1)
                    nzmax = min(nzmax - 3,0)
                    # Better version than the linear interpolation but much too slow
                    #fint = PchipInterpolator(np.log(self.var['P'][npmin:npmax,jys,ixs]),
                    #                     self.var[var][npmin:npmax,jys,ixs])
                    #new.var[var][:,jyt,ixt] = fint(np.log(p))
                    fint = interp1d(-self.var['Z'][nzmax:nzmin,jys,ixs],self.var[var][nzmax:nzmin,jys,ixs])
                    new.var[var][:,jyt,ixt] =  [fint(-zz) for zz in z]
                    ixt += 1
                jyt += 1
        return new

    def interpolPT(self,pt,varList='All',latRange=None,lonRange=None):
        """ interpolate the variables to a potential temperature level or a set of
            potential tempearture levels
            vars must be a list of variables or a single varibale
            pt must be a list of potential temperatures in K
        """
        if 'PT' not in self.var.keys():
            try:
                self._mkthet()
            except:
                print('missing P or T')
                return -1
        new = ECMWF_pure()
        if varList == 'All':
            varList = list(self.var.keys())
            try:
                varList.remove('SP')
                varList.remove('P')
            except: pass # in case these variables are already removed
        elif type(varList) == str:
            varList = [varList,]
        for var in varList:
            if var not in self.var.keys():
                print(var,' not defined')
                return
        if type(pt) in [float,int]:
            pt = [pt,]
        ptrev = [-x for x in pt]
        #print(ptrev)
        # first determine the boundaries of the domain
        if (latRange == []) | (latRange == None):
            nlatmin = 0
            nlatmax = self.nlat
        else:
            nlatmin = np.abs(self.attr['lats']-latRange[0]).argmin()
            nlatmax = np.abs(self.attr['lats']-latRange[1]).argmin()+1
        if (lonRange == []) | (lonRange == None):
            nlonmin = 0
            nlonmax = self.nlon
        else:
            nlonmin = np.abs(self.attr['lons']-lonRange[0]).argmin()
            nlonmax = np.abs(self.attr['lons']-lonRange[1]).argmin()+1
        new.attr['lats'] = self.attr['lats'][nlatmin:nlatmax]
        new.attr['lons'] = self.attr['lons'][nlonmin:nlonmax]
        new.nlat = len(new.attr['lats'])
        new.nlon = len(new.attr['lons'])
        new.attr['Lo1'] = new.attr['lons'][0]
        new.attr['Lo2'] = new.attr['lons'][-1]
        new.attr['La1'] = new.attr['lats'][0]
        new.attr['La2'] = new.attr['lats'][-1]
        new.attr['dlo'] = self.attr['dlo']
        new.attr['dla'] = self.attr['dla']
        new.project = self.project     
        new.date = self.date
        new.nlev = len(pt)
        new.attr['levtype'] = 'potential temperature'
        new.attr['levs'] = pt
        new.attr['plev'] = MISSING
        thetmin = np.min(pt)
        thetmax = np.max(pt)
        for var in varList:
            new.var[var] = np.empty(shape=(len(pt),nlatmax-nlatmin,nlonmax-nlonmin))
            jyt = 0
            # big loop that should be paralellized or calling numa for good performance
            for jys in range(nlatmin,nlatmax):
                ixt = 0
                for ixs in range(nlonmin,nlonmax):
                    # find the range of pt in the column assumed ordered from top to bottom
                    npup = np.abs(self.var['PT'][:,jys,ixs]-thetmax).argmin()
                    npbot = np.abs(self.var['PT'][:,jys,ixs]-thetmin).argmin()+1
                    npup = max(npup - 1,0)
                    npbot = min(npbot + 1,self.nlev)
                    #if jyt == 20 : print(npup,npbot)
                    # Test sorting and interpolation with reverse potential temperature to ensure growth
                    # along x-axis
                    if np.any(self.var['PT'][npup+1:npbot,jys,ixs]-self.var['PT'][npup:npbot-1,jys,ixs]>0):
                        #print('sort for ',ixt,jyt)
                        pts = np.sort(-self.var['PT'][npup:npbot,jys,ixs])
                        #print(pts)
                    else:
                        pts = -self.var['PT'][npup:npbot,jys,ixs]
                    new.var[var][:,jyt,ixt] = np.interp(ptrev,pts,self.var[var][npup:npbot,jys,ixs])
                    ixt += 1
                jyt += 1
        return new

    # CHECK CAREFULLY 
    # DOES NOT PERFORM AS EXPECTED
    # DO NOT USE FOR THE MOMENT
    def interpol_part(self,p,x,y,varList='All'):
        """ Interpolate the variables to the location of particles given by [p,y,x] using trilinear method."""
        if 'P' not in self.var.keys():
            self._mkp()
        if varList == 'All':
            varList = list(self.var.keys())
            varList.remove('SP')
            #varList.remove('P')
        elif type(varList) == str:
            varList = [varList,]
        for var in varList:
            if var not in self.var.keys():
                print(var,' not defined')
                return
        # Defines output dictionary
        result = {}
        # test whether fhyb already attached to the instance
        if not hasattr(self,'fhyb'):
            self.fhyb,void = tohyb('ERA5')
        # generate the 2D interpolator of the surface pressure
        lsp = RegularGridInterpolator((self.attr['lats'],self.attr['lons']),-np.log(self.var['SP']),bounds_error=True)
        lspi = lsp(np.transpose([y,x]))
        # define -log sigma = -log(p) - -log(ps)
        lsig = - np.log(p) - lspi
        # find the non integer hybrib level with offset due to the truncation of levels
        hyb = self.fhyb(np.transpose([lsig,lspi]))-self.attr['levs'][0]+1
        lhyb = np.floor(hyb).astype(np.int64)
        # flag non valid parcels outside the domain (see discussion in convsrc2)
        # notive that hyb == 100 not flagged as invalid, but clipped below
        non_valid = np.where(hyb>100)[0]
        # clip lhyb to remove out of bounds problems due to non valid parcels
        lhyb = np.clip(lhyb,0,99)
        #@@ test the extreme values of sigma end ps
        if np.min(lsig) < - np.log(0.95):
            print('large sigma detected ',np.exp(-np.min(lsig)))
        if np.max(lspi) > -np.log(45000):
            print('small ps detected ',np.exp(-np.max(lspi)))
        # Horizontal interpolation
        # clipping to avoid boundary effects for parcels just on the edges (especially on the eastern edge)
        ix = np.clip(np.floor((x-self.attr['Lo1'])/self.attr['dlo']).astype(np.int64),0,self.nlon-2)
        jy = np.clip(np.floor((y-self.attr['La1'])/self.attr['dla']).astype(np.int64),0,self.nlat-2)
        px = ((x - self.attr['Lo1']) %  self.attr['dlo'])/self.attr['dlo']
        py = ((y - self.attr['La1']) %  self.attr['dla'])/self.attr['dla']
        for var in varList:
            if len(self.var[var].shape)==3:
                vhigh = (1-px)*(1-py)*self.var[var][lhyb,jy,ix] + (1-px)*py*self.var[var][lhyb,jy+1,ix] \
                  + px*(1-py)*self.var[var][lhyb,jy,ix+1] + px*py*self.var[var][lhyb,jy+1,ix+1]
                vlow  = (1-px)*(1-py)*self.var[var][lhyb+1,jy,ix] + (1-px)*py*self.var[var][lhyb+1,jy+1,ix] \
                  + px*(1-py)*self.var[var][lhyb+1,jy,ix+1] + px*py*self.var[var][lhyb+1,jy+1,ix+1]
                hc = hyb - lhyb
                result[var] = (1-hc)*vhigh + hc*vlow
                result[var][non_valid] = MISSING
            if len(self.var[var].shape)==2:
                result[var] = (1-px)*(1-py)*self.var[var][jy,ix] + (1-px)*py*self.var[var][jy+1,ix] \
                  + px*(1-py)*self.var[var][jy,ix+1] + px*py*self.var[var][jy+1,ix+1]
        return result

    def _CPT(self):
        """ Calculate the cold point tropopause """
        if not set(['T','P']).issubset(self.var.keys()):
            print('T or P undefined')
            return
        levbnd = [30,90]
        # TODO: find a way to avoid the big loop
        self.var['pcold'] = np.empty(shape=(self.nlat,self.nlon))
        self.var['Tcold'] = np.empty(shape=(self.nlat,self.nlon))
        if 'Z' in self.var.keys():
            self.var['zcold'] = np.empty(shape=(self.nlat,self.nlon))
        # Calculate the cold point in the discrete profile
        # TO DO: make a smoother version with vertical interpolation
        nc = np.argmin(self.var['T'][levbnd[0]:levbnd[1],...],axis=0)\
           + levbnd[0]
        for jy in range(self.nlat):
            for ix in range(self.nlon):
                self.var['pcold'][jy,ix] = self.var['P'][nc[jy,ix],jy,ix]
                self.var['Tcold'][jy,ix] = self.var['T'][nc[jy,ix],jy,ix]
        if  'Z' in self.var.keys():
            for jy in range(self.nlat):
                for ix in range(self.nlon):
                    self.var['zcold'][jy,ix] = self.var['Z'][nc[jy,ix],jy,ix]
        return

    #@jit
    def _WMO(self,highlatOffset=False):
        """ Calculate the WMO tropopause
        When highlatoffset is true the 2K/km criterion is replaced by a 3K/km
        at high latitudes latitudes above 60S or 60N """
        if not set(['T','P']).issubset(self.var.keys()):
            print('T or P undefined')
            return
        self.var['pwmo'] = np.ma.empty(shape=(self.nlat,self.nlon))
        self.var['Twmo'] = np.ma.empty(shape=(self.nlat,self.nlon))
        if 'Z' in self.var.keys():
            self.var['zwmo'] = np.empty(shape=(self.nlat,self.nlon))
            zwmo = True
        else:
            zwmo = False
        levbnd = [30,90]
        highbnd = levbnd[0]
        lowbnd =  levbnd[1]
        logp = np.log(self.var['P'])
        # dz from the hydrostatic formula dp/dz = - rho g = - p/T g /R
        # dz = dz/dp p dlogp = - 1/T R/g dlogp (units m)
        # dz is shifted b one index position / T, p
        # dz[i] is the positive logp thickness for the layer between levels i and i+1, that is
        # above level i+1
        dz = cst.R/cst.g * self.var['T'][1:,:,:] * (logp[1:,:,:]-logp[:-1,:,:])
        # calculate dT/dz = - 1/p dT/dlogp rho g = - g/R 1/T dT/dlogp
        # dTdz[i] is the vertical derivative at level [i+1]
        #lapse = - cst.g/cst.R * (1/self.var['T'][highbnd+1:lowbnd-1,...]) * \
        #               (self.var['T'][highbnd:lowbnd-2,...] - self.var['T'][highbnd+2:lowbnd,...]) / \
        #               (logp[highbnd:lowbnd-2,...]-logp[highbnd+2:lowbnd,...])
        lapse = - cst.g/cst.R * (1/self.var['T'][highbnd+1:lowbnd,...]) * \
                       (self.var['T'][highbnd:lowbnd-1,...] - self.var['T'][highbnd+1:lowbnd,...]) / \
                       (logp[highbnd:lowbnd-1,...]-logp[highbnd+1:lowbnd,...])

        for jy in range(self.nlat):
            # standard wmo criterion
            offset = - 0.002
            thicktrop = 2000
            # adaptation of the WMO offset at high latitude
            if highlatOffset & (abs(self.attr['lats'][jy]) > 60): offset = -0.003
            for ix in range(self.nlon):
                # location of lapse rate exceeding the threshod
                slope = list(np.where(lapse[:,jy,ix] > offset)[0])
                Deltaz = 0.
                found = False
                # explore slope to find the first case where the slope is maintained
                # over two km
                # This is required to avoid shallow inversion layers to be confused with
                # the tropopause
                while not found:
                    if len(slope)>0: # should be done with try but forbidden with numba
                        # candidate tropopause
                        test = slope.pop()
                    else:
                        # if all slopes have been processed without finding a suitable tropopause, mask loc
                        self.var['pwmo'][jy,ix] = np.ma.masked
                        self.var['Twmo'][jy,ix] = np.ma.masked
                        if zwmo: self.var['zwmo'][jy,ix] = np.ma.masked
                        found = True
                        break
                    # location of the basis of the interval
                    # +1 to account for the shift of the finite difference
                    lev0 = test+1+highbnd
                    lev = lev0-1
                    Deltaz = dz[lev,jy,ix]
                    # performs search above the candidate tropopause
                    search = True
                    while Deltaz < thicktrop:
                        lev -= 1
                        Deltaz += dz[lev,jy,ix]
                        # mean slope over the considered layer
                        if (self.var['T'][lev,jy,ix]-self.var['T'][lev0,jy,ix])/Deltaz < offset:
                            search = False
                            break
                    if search:
                        found = True
                        self.var['pwmo'][jy,ix] = self.var['P'][lev0,jy,ix]
                        self.var['Twmo'][jy,ix] = self.var['T'][lev0,jy,ix]
                        if zwmo: self.var['zwmo'][jy,ix] = self.var['Z'][lev0,jy,ix]
        return

    def show(self,var,lev=0,cardinal_level=True,txt=None,log=False,clim=(None,None),figsize=(5,5),
         axf=None,cmap=mymap,savfile=None,cLines=None,show=True,scale=1,aspect=1,
         horizontal=False):
        """ Chart for data fields """
        # test existence of key field
        if var in self.var.keys():
            if len(self.var[var].shape) == 3:
                # detection of the level as pressure or potential temperature with a value always
                # larger than the number of levels (may fail at the top of the model in pressure or in z coordinate (if km))
                # in this cas, force the value with cardinal_level
                if (cardinal_level==False) | (lev > self.nlev-1):
                    clev = np.abs(np.array(self.attr['levs'])-lev).argmin()
                else:
                    clev = lev
                buf = self.var[var][clev,:,:]
            else:
                buf = self.var[var]
        else:
            print ('undefined field')
            return
        if log: buf = np.log(buf)
        fs=15
        if axf is None:
            fig, ax = plt.subplots(figsize=figsize,nrows=1,ncols=1)
        else: ax = axf
        # It is assumed that the values are valid on the grid
        extent = [self.attr['Lo1']-0.5*self.attr['dlo'], 
                  self.attr['Lo2']+0.5*self.attr['dlo'],
                  self.attr['La1']-0.5*self.attr['dla'], 
                  self.attr['La2']+0.5*self.attr['dla']]
        iax = ax.imshow(scale*buf, interpolation='nearest',
                        extent=extent,origin='lower', aspect=1,
                        cmap=cmap,clim=clim)

        if cLines is not None:
                 ax.contour(buf,extent=extent,levels=cLines,
                            origin='lower')
        
        ax.scatter([stations[self.project][0],],[stations[self.project][1],] \
                    ,color='r',marker='+',s=64,linewidth=9)
        
        if txt is None:
            ax.set_title(var+' lev'+str(lev),fontsize=fs)
        else:
            ax.set_title(txt,fontsize=fs)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        # plot adjusted colorbar
        #axpos = ax.get_position()
        #pos_x = axpos.x0 + axpos.width + 0.01
        #pos_cax = fig.add_axes([pos_x,axpos.y0,0.04,axpos.height])
        if horizontal:
            cbar=plt.colorbar(iax,location='bottom')
        else:
            pos_cax = ax.inset_axes([1.02,0,0.06,1])
            cbar=plt.colorbar(iax,cax=pos_cax)
        cbar.ax.tick_params(labelsize=fs)

        if savfile is not None:
            plt.savefig(savfile,dpi=300,bbox_inches='tight')
        if show: plt.show()
        return ax

    # TO BE UPDATED
    # NON TESTED
    def chartlonz(self,var,lat,levs=(None,None),txt=None,log=False,clim=(None,None),
             cmap=mymap,savfile=None,show=True,scale=1,figsize=(11,4),axf=None,ylabel=True):
        """ Plot a lon x alt section for a given latitude
        """
        if var not in self.var.keys():
            print("UNKNOWN VAR ", var)
            return -1
        try:
            pos=np.where(self.attr['lats']>=lat)[0][0]
            #print('pos',pos)
        except:
            print('lat out of range')
            return -1
        if axf is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        else:
            ax = axf
        fs = 16
        if levs[0]==None: l1=29
        else: l1 = levs[0]
        if levs[1]==None: l2=115
        else: l2 = levs[1]
        try: dlo = self.attr['dlo']
        except: dlo =  (self.attr['lons'][-1] - self.attr['lons'][0])/(len(self.attr['lons'])-1)
        lons = np.arange(self.attr['lons'][0]-0.5*dlo,self.attr['lons'][-1]+dlo,dlo)
        #if Z variable is available, lets use it, if not use zscale
        try:
            zz1 = 0.5*(self.var['Z'][l1-1:l2+1,pos, :] + self.var['Z'][l1:l2+2,pos,:])/1000
            zz = np.empty((zz1.shape[0],zz1.shape[1]+1))
            zz[:,1:-1] = 0.5*(zz1[:,1:]+zz1[:,:-1])
            zz[:,0] = zz1[:,0]
            zz[:,-1] = zz1[:,-1]
            #print(zz.shape,len(lons))
            iax=ax.pcolormesh(lons,zz,scale*self.var[var][l1:l2+1,pos,:],
                            vmin=clim[0],vmax=clim[1],cmap=cmap)
            if ylabel: plt.ylabel('altitude (km)',fontsize=fs)
            #print('USE Z')
        except:
            iax=ax.pcolormesh(lons,self.attr['zscale_i'][l1:l2+2],self.var[var][l1:l2+1,pos, :],
                        vmin=clim[0],vmax=clim[1],cmap=cmap)

            if ylabel: plt.ylabel('baro altitude (km)',fontsize=fs)
        ax.tick_params(labelsize=fs)
        plt.xlabel('longitude',fontsize=fs)

        if txt is None:
            plt.title(var+' lat'+str(lat),fontsize=fs)
        else:
            plt.title(txt,fontsize=fs)
        #cax = fig.add_axes([0.91, 0.21, 0.03, 0.6])
        #cbar = fig.colorbar(iax,cax=cax)
        cbar = plt.colorbar(iax,ax=ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=fs)
        
        if savfile is not None:
            plt.savefig(savfile,bbox_inches='tight',dpi=300)
        if show: plt.show()
        return ax

    # TO BE UPDATED
    # NON TESTED
    def chartlatz(self,var,lon,levs=(None,None),txt=None,log=False,clim=(None,None),
             cmap=mymap,savfile=None,show=True,scale=1,figsize=(11,4),axf=None,ylabel=True):
        """ Plot a lon x alt section for a given latitude
        """
        if var not in self.var.keys():
            print("UNKNOWN VAR ", var)
            return -1
        try:
            pos=np.where(self.attr['lons']>=lon)[0][0]
            #print('pos',pos)
        except:
            print('lon out of range')
            return -1
        if axf is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        else:
            ax = axf
        fs = 16
        if levs[0]==None: l1=29
        else: l1 = levs[0]
        if levs[1]==None: l2=115
        else: l2 = levs[1]
        try: dla = self.attr['dla']
        except: dla =  (self.attr['lats'][-1] - self.attr['lats'][0])/(len(self.attr['lats'])-1)
        lats = np.arange(self.attr['lats'][0]-0.5*dla,self.attr['lats'][-1]+dla,dla)
        #if Z variable is available, lets use it, if not use zscale
        try:
            zz1 = 0.5*(self.var['Z'][l1-1:l2+1,:,pos] + self.var['Z'][l1:l2+2,:,pos])/1000
            zz = np.empty((zz1.shape[0],zz1.shape[1]+1))
            zz[:,1:-1] = 0.5*(zz1[:,1:]+zz1[:,:-1])
            zz[:,0] = zz1[:,0]
            zz[:,-1] = zz1[:,-1]
            #print(zz.shape,len(lons))
            iax=ax.pcolormesh(lats,zz,scale*self.var[var][l1:l2+1,:,pos],
                            vmin=clim[0],vmax=clim[1],cmap=cmap)
            if ylabel: plt.ylabel('altitude (km)',fontsize=fs)
            #print('USE Z')
        except:
            iax=ax.pcolormesh(lats,self.attr['zscale_i'][l1:l2+2],self.var[var][l1:l2+1,:,pos],
                        vmin=clim[0],vmax=clim[1],cmap=cmap)

            if ylabel: plt.ylabel('baro altitude (km)',fontsize=fs)
        ax.tick_params(labelsize=16)
        plt.xlabel('latitude',fontsize=fs)

        if txt is None:
            plt.title(var+' lon'+str(lon),fontsize=fs)
        else:
            plt.title(txt,fontsize=fs)
        #cax = fig.add_axes([0.91, 0.21, 0.03, 0.6])
        #cbar = fig.colorbar(iax,cax=cax)
        cbar = plt.colorbar(iax,ax=ax)
        cbar.ax.tick_params(labelsize=fs)
        if savfile is not None:
            plt.savefig(savfile,bbox_inches='tight',dpi=300)
        if show: plt.show()
        return ax
    
# standard class to read data
class ECMWF(ECMWF_pure):
# to do, raise a specific exception in case of error 
# instead of ending with a return
    
    def __init__(self,project,date,exp=[None]):
        ECMWF_pure.__init__(self)
        self.project = project
        self.date = date
        self.exp = exp
        #provision for surface data to be added as a separate file
        self.sfc = False
        
        if project not in ['OHP','OPAR', 'ACARR']:
            print('Unknown project')
            return
        if date.hour not in [0,6,12,18]:
            print('This hour is not available, only 0, 6, 12 and 18')
            return
        # index corresponding to the hour 
        self.ih = int(0.001+date.hour/6)
        
        if 'gort' == socket.gethostname():
            self.rootdir = '/dkol/data/OPZ/'+project
        elif 'spirit' in socket.gethostname():
            self.rootdir = '/proju/flexpart/flexpart_in/OPZ/'+project
        elif 'spirit' in socket.gethostname():
            self.rootdir = '/STC/ERA5'
        elif 'satie' in socket.gethostname():
            self.rootdir = '/dsk2/OPZ/'+project
        elif 'Mentat' == socket.gethostname():
            self.rootdir = 'C:\\cygwin64\\home\\berna\\data\\OPZ\\'+project

        
        if project == 'ACARR':
            self.expected_vars = ['U','V','W','T','VO','O3','Q','LNSP']
        else:
            self.expected_vars = ['U','V','T','VO','O3','LNSP']

        # opening the main file (the only one for the moment)
        try:
            self.ncid = Dataset(os.path.join(self.rootdir,date.strftime('%Y/OPZLWDA%Y%m%d.nc')))
        except:
            print(date.strftime('cannot find OPZLWDA%Y%m%d.nc'),' for ',project)
            return
            
        try:
            self.var['SP'] = np.exp(self.ncid.variables['lnsp'][self.ih,0,...])
            self.nlon = self.ncid.dimensions['longitude'].size
            self.nlat = self.ncid.dimensions['latitude'].size
            self.attr['lons'] = self.ncid.variables['longitude'][:].data          
            self.attr['lats'] = self.ncid.variables['latitude'][:].data
        except:
            print('Cannot read pressure and metadata')
            return
        
        # longitude and latitude interval
        self.attr['dlo'] = (self.attr['lons'][-1] - self.attr['lons'][0]) / (self.nlon-1)
        self.attr['dla'] = (self.attr['lats'][-1] - self.attr['lats'][0]) / (self.nlat-1)
        self.attr['Lo1'] = self.attr['lons'][0]
        self.attr['Lo2'] = self.attr['lons'][-1]
        if self.attr['Lo1']>self.attr['Lo2']:
            self.attr['Lo1'] = self.attr['Lo1']-360.
        self.attr['levtype'] = 'hybrid'
        # The list of levels is complete
        # If not, it should be changed here
        self.nlev = 137
        self.attr['levs'] = np.arange(137)
        # List of reference pressure levels
        # duplicated in _mkpscale below
        self.attr['plev'] = am + bm * cst.pref
 
        # Reverting lat order
        self.attr['lats'] = self.attr['lats'][::-1]
        self.var['SP']   = self.var['SP'][::-1,:]
        self.attr['dla'] = - self.attr['dla']
        self.attr['La1'] = self.attr['lats'][0]
        self.attr['La2'] = self.attr['lats'][-1]
        
        return
    
    def close(self):
        self.ncid.close()

# short cut for some common variables
    def _get_T(self):
        self._get_var('T')
    def _get_U(self):
        self._get_var('U')
    def _get_V(self):
        self._get_var('V')
    def _get_W(self):
        self._get_var('W')
    def _get_Q(self):
        self._get_var('Q')

# get a variable from the archive
    def _get_var(self,var):
        if (var in self.var.keys()):
            print('The variable is already read')
            return
        if (var not in self.expected_vars):
            print('This variable is not available')
            return
        try:
            self.var[var] = self.ncid.variables[dvar[var][0]][self.ih,...]
        except:
            print('Cannot read ',var)
            return
        
        # Reverting North -> South to South to North
        self.var[var] = self.var[var][:,::-1,:]
        return None

    def get_var(self,var):
        self._get_var(var)
        return self.var['var' ]
    
    # Methods to calculate new variables

    def _mkp(self):
        # Calculate the pressure field
        self.var['P'] =  np.empty(shape=(self.nlev,self.nlat,self.nlon))
        for i in range(self.nlev):
            self.var['P'][i,:,:] = am[i] \
                                 + bm[i]*self.var['SP']

    def _mkpz(self):
        # Calculate pressure field for w, which is defined on half levels 
        # starting from the exterior (check)
        # The vertical dimension of al is one unit larger than that of am
        self.var['PZ'] =  np.empty(shape=(self.nlev,self.nlat,self.nlon))
        for i in range(self.nlev):
            self.var['PZ'][i,:,:] = al[i] \
                                  + bl[i]*self.var['SP']

    def _mkpscale(self):
        # Define the standard pressure scale for this vertical grid
        self.attr['pscale'] = am + bm * 101325
        #self.attr['pscale_i'] = self.attr['ai'] + self.attr['bi'] * 101325

    def _mkzscale(self):
        # Define the standard altitude scale
        from zISA import zISA
        try:
            self.attr['zscale'] = zISA(self.attr['pscale'])
            #self.attr['zscale_i'] = zISA(self.attr['pscale_i'][1:])
            #ztop = 2*self.attr['zscale_i'][0]-self.attr['zscale_i'][1]
            #self.attr['zscale_i'] = np.concatenate(((ztop,),self.attr['zscale_i']))
        except:
            print('pscale must be generated before zscale')

    def _mkthet(self,suffix=''):
        # Calculate the potential temperature
        if not set(['T'+suffix,'P'+suffix]).issubset(self.var.keys()):
            print('T or P undefined')
            return
        self.var['PT'+suffix] = self.var['T'+suffix] * (cst.p0/self.var['P'+suffix])**cst.kappa

    def _mkrhi(self,suffix=''):
        # Calculate the relative humidity with respect to the ice
        if not set(['Q'+suffix,'T'+suffix,'P'+suffix]).issubset(self.var.keys()):
            print('T, P or Q undefined')
            return
        self.var['RHI'+suffix] = self.var['Q'+suffix]/satratio.satratio(self.var['P'+suffix],self.var['T'+suffix])

    def _mkrhl(self,suffix=''):
        # Calculate the relative humidity with respect to the iceliquid water
        if not set(['Q'+suffix,'T'+suffix,'P'+suffix]).issubset(self.var.keys()):
            print('T, P or Q undefined')
            return
        self.var['RHL'+suffix] = self.var['Q'+suffix]/satratio.liquid_satratio(self.var['P'+suffix],self.var['T'+suffix])
   
    def _mkthetei(self,suffix=''):
        # Calculate the equivalent potential temperature with respect to ice
        # Simplified version where r_T is not accounted
        if not set(['PT'+suffix,'Q'+suffix,'RHI'+suffix]).issubset(self.var.keys()):
            print('PT, RHI or Q undefined')
            return
        # views 
        r = self.var['Q'+suffix]
        T = self.var['T'+suffix]
        PT = self.var['PT'+suffix]
        self.var['PTEI'+suffix] = PT * np.exp(cst.Li(T) * r / (cst.Cp * T)) \
            * self.var['RHI'+suffix]**(-r * cst.Rv / cst.Cp)

    def _mktheteistar(self,suffix=''):
        # Calculate the equivalent saturated potential temperature with respect to ice
        # Simplified version where r_T is not accounted
        if not set(['PT'+suffix]).issubset(self.var.keys()):
            print('PT undefined')
            return
        r = satratio.satratio(self.var['P'+suffix],self.var['T'+suffix])
        T = self.var['T'+suffix]
        PT = self.var['PT'+suffix]
        self.var['PTEI*'+suffix] = PT * np.exp(cst.Li(T) * r / (cst.Cp * T)) 

    def _mkrho(self,suffix=''):
        # Calculate the dry density
        if not set(['T'+suffix,'P'+suffix]).issubset(self.var.keys()):

            return
        self.var['RHO'+suffix] = (1/cst.R) * self.var['P'+suffix] / self.var['T'+suffix]
        return

    def _mkBV(self):
        # Calculate the (dry) Brünt-Vaissale frequency
        if 'PT' not in self.var.keys():
            self._mkthet()
        if 'Z' not in self.var.keys():
            self._mkz()
        self.var['N2'] =  np.empty(shape=self.var['T'].shape)
        self.var['N2'][1:-1,...] = cst.g / self.var['PT'][1:-1,...] \
            * (self.var['PT'][2:,...] - self.var['PT'][:-2,...]) \
            / (self.var['Z'][2:,...] - self.var['Z'][:-2,...])
        # repeat N2 at top and bottom
        self.var['N2'][0,...] = self.var['N2'][1,...]
        self.var['N2'][-1,...] = self.var['N2'][-2,...]
        return
    
    def _mBVmi(self):
        """ Calculate N2m according to the formula (36) of 
        Durran & Kemp, JAS, 1982 using ice saturation"""
        # First calculate qs(T,P)
        qs = satratio.satratio(self.var['T'],self.var['P'])
        Li = cst.Li(self.var['T'])
        f1 = (1 + Li/(cst.R * self.var['T'])) \
            /(1 + cst.epsilon * Li**2 / (cst.Cp * cst.R * self.var['T']**2))
        f2 = np.empty(shape = f1.shape)
        f2[1:-1,...] = ((self.var['PT'][2:,...] - self.var['PT'][:-2,...]) \
                      / self.var['PT'][1:-1,...] \
                      + Li / (cst.R * self.var['T']) \
                      * (qs[2:,...] - qs[:-2,...])) \
                      / (self.var['Z'][2:,...] - self.var['Z'][:-2,...])
        dqsdz = (qs[2:,...] - qs[:-2,...]) \
              / (self.var['Z'][2:,...] - self.var['Z'][:-2,...])
        return( cst.g * (f1*f2 - dqsdz))
    
    def _mkrhoq(self):
        # Calculate the moist density
        if not set(['T','P','Q']).issubset(self.var.keys()):
            print('T, P or Q undefined')
            return
        pcor = 230.617*self.var['Q']*np.exp(17.5043*self.var['Q']/(241.2+self.var['Q']))
        self.var['RHOQ'] = (1/cst.R) * (self.var['P'] - pcor) / self.var['T']
        return
   
    def _checkThetProfile(self):
        # Check that the potential temperature is always increasing with height
        if 'PT' not in self.var.keys():
            print('first calculate PT')
            return
        ddd = self.var['PT'][:-1,:,:]-self.var['PT'][1:,:,:]
        for lev in range(ddd.shape[0]):
            if np.min(ddd[lev,:,:]) < 0:
                print('min level of inversion: ',lev)
                return lev
        return

    def _mkz(self,moist=False):
        """ Calculate the geopotential altitude (m) without taking moisture into account """
        if not set(['T','P']).issubset(self.var.keys()):
            print('T or P undefined')
            return
        if 'Z0' not in self.var.keys():
            try:
                nz0 = Dataset(os.path.join(self.rootdir,'GeopotZ0.nc'))
                # The variable stored is the geopotential, not z
                self.var['Z0'] = nz0.variables['z'][0,::-1,:]/cst.g
                nz0.close()
            except:
                print('Z0 not found')
                return
        if moist:
            if 'Q' not in self.var.keys():
                self._get_var('Q')
            q = 'Q'
            ff = (1 + self.var['Q']/cst.epsilon)/(1 + self.var['Q'])
        else: 
            q = ''
            ff = np.ones(shape=self.var['T'].shape)
        self.var['Z'+q] =  np.empty(shape=self.var['T'].shape)
        uu = - np.log(self.var['P'])
        uusp = -np.log(self.var['SP'])
        self.var['Z'+q][self.nlev-1,:,:] = self.var['Z0'] \
               + (cst.R/cst.g) * self.var['T'][self.nlev-1,:,:] * ff[self.nlev-1,:,:]\
               * (uu[self.nlev-1,:,:]-uusp)
        for i in range(self.nlev-2,-1,-1):
             self.var['Z'+q][i,:,:] = self.var['Z'][i+1,:,:] + 0.5*(cst.R/cst.g) \
                    * (self.var['T'][i,:,:] * ff[i,:,:] + self.var['T'][i+1,:,:] * ff[i+1,:,:])\
                    * (uu[i,:,:]-uu[i+1,:,:])
        return

    def _mkpv(self,suffix=''):
        """ Calculate the potential vorticity using the isentropic formula """
        if not set(['PT'+suffix,'P'+suffix,'U'+suffix,'V'+suffix,'VO'+suffix]).issubset(self.var.keys()):
            print('at least one of P, PT, U, V or VO is undefined')
            return
        # Calculation of the correction to the vorticity assuming that the grid is global
        dlat = np.deg2rad(self.attr['dla'])*cst.REarth
        dPTdPhi = np.empty(shape = self.var['PT'+suffix].shape)
        dPTdPhi[:,1:-1,:] = 0.5*(self.var['PT'+suffix][:,2:,:]-self.var['PT'+suffix][:,:-2,:])/dlat
        dPTdPhi[:,0,:] = (self.var['PT'+suffix][:,1,:]-self.var['PT'+suffix][:,0,:])/dlat
        dPTdPhi[:,-1,:] = (self.var['PT'+suffix][:,-1,:]-self.var['PT'+suffix][:,-2,:])/dlat
        dlon = np.deg2rad(self.attr['dlo'])*cst.REarth
        coslat = np.cos(np.deg2rad(self.attr['lats']))
        dPTdLam = np.empty(shape = self.var['PT'+suffix].shape)
        dPTdLam[...,1:-1] = 0.5*(self.var['PT'+suffix][...,2:]-self.var['PT'+suffix][...,:-2])/dlon
        dPTdLam[...,0] = (self.var['PT'+suffix][...,1]-self.var['PT'+suffix][...,0])/dlon
        dPTdLam[...,-1] = (self.var['PT'+suffix][...,-1]-self.var['PT'+suffix][...,-2])/dlon
        dPTdLam /= coslat[np.newaxis,:,np.newaxis]
        logP = np.log(self.var['P'+suffix])
        dudP = np.empty(shape = self.var['U'+suffix].shape)
        dvdP = np.empty(shape = self.var['V'+suffix].shape)
        dPTdP = np.empty(shape = self.var['PT'+suffix].shape)
        dudP[1:-1,...] = (self.var['U'+suffix][2:,...]-self.var['U'+suffix][:-2,...])/\
                             (self.var['P'+suffix][1:-1,...]*(logP[2:,...]-logP[:-2,...]))
        dvdP[1:-1,...] = (self.var['V'+suffix][2:,...]-self.var['V'+suffix][:-2,...])/\
                            (self.var['P'+suffix][1:-1,...]*(logP[2:,...]-logP[:-2,...]))
        dPTdP[1:-1,...] = (self.var['PT'+suffix][2:,...]-self.var['PT'+suffix][:-2,...])/\
                             (self.var['P'+suffix][1:-1,...]*(logP[2:,...]-logP[:-2,...]))
        # Do not loose time as PV at the ground or the top levl is not of any interest
        dudP[0,...] = dudP[1,...]
        dvdP[0,...] = dvdP[1,...]
        dPTdP[0,...] = dPTdP[1,...]
        dudP[-1,...] = dudP[-2,...]
        dvdP[-1,...] = dvdP[-2,...]
        dPTdP[-1,...] = dPTdP[-2,...]
        flat = 2*cst.Omega * np.sin(np.deg2rad(self.attr['lats']))
        self.var['PV'+suffix] = - cst.g * (dPTdP * (self.var['VO'+suffix] + flat[np.newaxis,:,np.newaxis]) \
                                       + dudP * dPTdPhi - dvdP * dPTdLam)
        del dPTdPhi, dPTdLam, logP, dudP, dvdP, dPTdP
        return

if __name__ == '__main__':
    """ test procedure """
    from datetime import datetime
    date = datetime(2022,12,9,12)
    dat = ECMWF('ACARR',date)   
    dat._get_T()
    dat._get_U()
    dat._get_V()
    dat._get_Q()
    dat._get_var('W')
    dat._get_var('VO')
    dat._get_var('O3')
    dat._mkp()
    dat._mkz()
    dat._mkz(moist=True)
    dat._mkrho()
    dat._mkrhoq()
    dat._mkthet()
    dat._checkThetProfile()
    dat._mkBV()
    dat._mkpv()
    dat._CPT()
    dat._WMO()
    dat._mkrhi()
    dat._mkthetei()
    #dat._mktheteistar()
    prof = dat.get_station_profile()
    prof['N'] = np.sqrt(prof['N2'])
    
    # Plot the profile
    fig, axs = plt.subplots(figsize=(12,12),nrows=3,ncols=3,sharey=True,squeeze=True)
    Z = prof['Z']/1000
    # Temperature
    axs[0,0].plot(prof['T'],Z)
    axs[0,0].annotate('CPT {:0.2f} km'.format(prof['zcold']/1000),(240,16))
    axs[0,0].annotate('WMO {:0.2f} km'.format(prof['zwmo']/1000),(240,21))
    # Brünt-Vaissala frequency
    axs[0,2].plot(prof['N'],Z)
    # Potential temperatures
    axs[0,1].plot(prof['PT'],Z)
    axs[0,1].plot(prof['PTEI'],Z)
    #axs[0,1].plot(prof['PTEI*'],prof['Z']/1000)
    # Vertical velocity
    axs[1,0].plot(-prof['W']/(cst.g*prof['RHOQ']),Z)
    # Horizontal velocities
    axs[1,1].plot(prof['U'],Z)
    axs[1,1].plot(prof['V'],Z)
    axs[1,1].plot(np.sqrt(prof['U']**2+prof['V']**2),Z)
    # Vertical gradient of velocity
    prof['DUDZ'] = (prof['U'][2:]-prof['U'][:-2])/(Z[2:]-Z[:-2])/1000
    prof['DVDZ'] = (prof['V'][2:]-prof['V'][:-2])/(Z[2:]-Z[:-2])/1000
    axs[1,2].plot(prof['DUDZ'],Z[1:-1])
    axs[1,2].plot(prof['DVDZ'],Z[1:-1])
    # Moisture (relative humidity / ice)
    axs[2,0].plot(prof['RHI'],Z)
    # Lait potential vorticity
    axs[2,1].plot(1.e6*prof['PV']*(600/prof['PT'])**4,Z)
    # Ozone
    axs[2,2].plot(1.e6*prof['O3'],Z)
    # Beauty
    axs[0,0].set_xlabel('Temperature (K)')
    axs[0,1].set_xlabel('Potential temperatures (K)')
    axs[0,2].set_xlabel('N (rad/s)')
    axs[1,0].set_xlabel('W (m/s)')
    axs[1,1].set_xlabel('Horizontal velocity (m/s)')
    axs[1,1].legend([r'$U$',r'$V$',r'$\sqrt{U^2+V^2}$'])
    axs[1,2].set_xlabel('Vertical shear (1/s)')
    axs[1,2].legend([r'$\frac{\partial U}{\partial Z}$',r'$\frac{\partial V}{\partial Z}$'])
    axs[2,0].set_xlabel('Relative humidity / ice')
    axs[2,1].set_xlabel('Lait PV')
    axs[2,2].set_xlabel('Ozone (ppmm)')
    for j in range(3):
        axs[j,0].set_ylabel('Altitude (km)')
        for i in range(3):
            axs[j,i].grid('True')
    fig.suptitle(date.strftime('ECMWF operational analysis at ACARR %d-%m-%Y %H:00 UTC'),fontsize=22)
    plt.show()