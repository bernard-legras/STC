#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Opens a FLEXPART type grib file from ECMWF and reads meteorological fields.
Requires pygrib, loaded from conda-forge
For both python 2 and python 3 under anaconda
ACHTUNG: in some installations of python 3, error at loading stage because
of libZ not found. Remedy: load netCDF4 in the calling program even if not used.

Functions: Open files, read the data, extract subgrids, make charts, interpolate in time,
interpolate to pressure levels

Data are read on hybrid levels. They can interpolated to pressure levels. Do it on subgrids
as it is a time consuming procedure.

Usage:
>> from ECMWF_N import ECMWF
open files for a date (datetime) and a project (VOLC or STC)
>> data = ECMWF(project,date)
read a variable var
>> data._get_var(var)
example: data._get_var('CC') for cloud cover (see lists in ECMWF class)
calculate pressure field
>> data._mkp()
calculate potential temperature (requires to read T and calculate P before)
>> data._mktheta()
extract a subgrid
>> data1 = data.extract(latRange=[lat1,lat2],lonRange=[lon1,lon2])
plot a chart for the variable var and level lev (lev is the level number)
>> data.chart(var,lev)
interpolate in time (date2 must be between the dates of data0 and data1)
>> data2 = data0.interpol-time(data1,date2)
interpolate to pressure levels
>> data1 = data.interpolP(pList,varList,latRange,lonRange)
where pList is a pressure or a list of pressures, varList is a variable or a list of variables,
latRange and lonRange as in extract

The ECMWF class is used to read the file corresponding to a date for several projects/
The relevant projects are STC and VOLC

The ECMWF_pure allows to define a template object without reading files.
It can be used to modify data. It is produced as an output of the interpolation.

Created on 21/01/2018 from ECMWF.py

@author: Bernard Legras (legras@lmd.ens.fr)
@licence: CeCILL-C
"""
from datetime import datetime, timedelta
import numpy as np
import math
import pygrib
import os
#from cartopy import feature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.crs as ccrs
#import matplotlib.pyplot as plt
from cartopy import feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import socket
from scipy.interpolate import interp1d,RegularGridInterpolator
from mki2d import tohyb
import constants as cst
import gzip,pickle
from numba import jit
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

class curtain(object):
    def __init__(self):
        self.var={}
        self.attr={}

# template object produced by extraction and interpolation
class ECMWF_pure(object):
    def __init__(self):
        self.var={}
        self.attr={}
        self.d2d={}
        self.d1d={}
        self.warning = []

    def show(self,var,lev=0,cardinal_level=True,txt=None,log=False,clim=(None,None),figsize=(11,4),
             axf=None,cmap=mymap,savfile=None,cLines=None,show=True,scale=1,aspect=1,projec=None,
             sat_H=35785831,xylim=False,polar=False):
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
        elif var in self.d2d.keys():
            buf = self.d2d[var]
        else:
            print ('undefined field')
            return
        fs=15
        # it is unclear how the trick with cm_lon works in imshow but it does
        # the web says that it is tricky to plot data accross dateline with cartopy
        # check https://stackoverflow.com/questions/47335851/issue-w-image-crossing-dateline-in-imshow-cartopy
        cm_lon =0
        ctrl_lon=(self.attr['lons'][0]+self.attr['lons'][-1])/2
        ctrl_lat=(self.attr['lats'][0]+self.attr['lats'][-1])/2
        # guess that we want to plot accross dateline
        if self.attr['lons'][-1] > 180: cm_lon=180
        projplate = ccrs.PlateCarree(central_longitude=cm_lon)
        if polar: ctrl_lat = 90*np.sign(ctrl_lat)
        if projec == 'ortho':
            proj=ccrs.Orthographic(central_longitude=ctrl_lon,central_latitude=ctrl_lat)
        elif projec == 'azimuthalequi':
            proj=ccrs.AzimuthalEquidistant(central_longitude=ctrl_lon,central_latitude=ctrl_lat)
        elif projec == 'nearside':
            proj=ccrs.NearsidePerspective(central_longitude=ctrl_lon,central_latitude=ctrl_lat,satellite_height=sat_H)
        elif projec == 'lambert':
            proj=ccrs.LambertConformal(central_longitude=ctrl_lon,central_latitude=ctrl_lat,cutoff=self.attr['lats'][0],)
        elif projec is None:
            proj=projplate
        if figsize is not None:
            fig = plt.figure(figsize=figsize)
            fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
        if axf is None: ax = plt.axes(projection = proj)
        else: ax = axf
        iax = ax.imshow(scale*buf, transform=projplate, interpolation='nearest',
                        extent=[self.attr['lons'][0]-cm_lon, self.attr['lons'][-1]-cm_lon,
                                self.attr['lats'][0], self.attr['lats'][-1]],
                        origin='lower', aspect=aspect,cmap=cmap,clim=clim)
        if xylim:
            # if cm_lon = 180, the shift must be 360, do not seek why
            x1a,y1a = proj.transform_point(self.attr['lons'][0]-2*cm_lon,self.attr['lats'][0],ccrs.Geodetic())
            x1b,y1b = proj.transform_point(self.attr['lons'][0]-2*cm_lon,self.attr['lats'][-1],ccrs.Geodetic())
            x1c,y1c = proj.transform_point(ctrl_lon-2*cm_lon,self.attr['lats'][0],ccrs.Geodetic())
            x1 = min(x1a,x1b,x1c)
            y1 = min(y1a,y1b,y1c)
            x2a,y2a = proj.transform_point(self.attr['lons'][-1]-2*cm_lon,self.attr['lats'][-1],ccrs.Geodetic())
            x2b,y2b = proj.transform_point(self.attr['lons'][-1]-2*cm_lon,self.attr['lats'][0],ccrs.Geodetic())
            x2c,y2c = proj.transform_point(ctrl_lon-2*cm_lon,self.attr['lats'][-1],ccrs.Geodetic())
            x2 = max(x2a,x2b,x2c)
            y2 = max(y2a,y2b,y2c)
            #print(x1a,y1a,x1b,y1b,x2a,y2a,x2b,y2b)
            ax.set_xlim(x1,x2)
            ax.set_ylim(y1,y2)
        if polar:
            yy = np.empty(len(self.attr['lons']))
            if ctrl_lat>0: yy.fill(self.attr['lats'][0])
            else: yy.fill(self.attr['lats'][-1])
            x1,y1,_ = (proj.transform_points(ccrs.Geodetic(),self.attr['lons']-2*cm_lon,yy)).T
            del yy
            ax.set_xlim(np.min(x1),np.max(x1))
            ax.set_ylim(np.min(y1),np.max(y1))
        xlocs = None
        if cm_lon == 180:
                interx = 30
                # next multiple of interx on the east of the western longitude boundary
                minx = self.attr['lons'][0] + interx - self.attr['lons'][0]%interx
                xlocs = list(np.arange(minx,181,interx))+list(np.arange(interx-180,self.attr['lons'][-1]-360,interx))
        if projec in ['ortho','azimuthalequi','nearside','lambert']:
            gl = ax.gridlines(draw_labels=True)
        elif projec in [None,'mercator','plate']:
            gl = ax.gridlines(draw_labels=True, xlocs=xlocs,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        if cLines is not None:
                 ax.contour(self.var[var][lev,:,:],transform=projplate,extent=(self.attr['lons'][0]-cm_lon,
                               self.attr['lons'][-1]-cm_lon,self.attr['lats'][0],self.attr['lats'][-1]),levels=cLines,origin='lower')
        ax.add_feature(feature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='50m',
                facecolor='none'))
        ax.coastlines('50m')
            #ax.add_feature(feature.BORDERS)
            # The grid adjusts automatically with the following lines
            # If crossing the dateline, superimposition of labels there
            # can be suppressed by specifying xlocs

        #ax.set_xlim(self.attr['lons'][0],self.attr['lons'][-1])
        #ax.set_ylim(self.attr['lats'][0],self.attr['lats'][-1])
        if projec == None:
            gl.top_labels = False
            gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': fs}
        gl.ylabel_style = {'size': fs}
        #gl.xlabel_style = {'color': 'red', 'weight': 'bold'}

        if txt is None:
            plt.title(var+' lev'+str(lev),fontsize=fs)
        else:
            plt.title(txt,fontsize=fs)
        # plot adjusted colorbar
        #axpos = ax.get_position()
        #pos_x = axpos.x0 + axpos.width + 0.01
        #pos_cax = fig.add_axes([pos_x,axpos.y0,0.04,axpos.height])
        pos_cax = ax.inset_axes([1.02,0,0.06,1])
        cbar=plt.colorbar(iax,cax=pos_cax)
        cbar.ax.tick_params(labelsize=fs)

        if savfile is not None:
            plt.savefig(savfile,dpi=300,bbox_inches='tight')
        if show: plt.show()
        ax.cm_lon = cm_lon
        return ax

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

    def shift2zero(self):
        # generate an object with the same content as self but with the origin
        # of longitude shifted to 0
        new = ECMWF_pure()
        # find the location of 0 in the initial grid
        ll0 = np.abs(self.attr['lons']).argmin()
        if ll0==0:
            print ('This dataset is already gridded from zero longitude')
            return None

        new.attr['lons'] = np.concatenate((self.attr['lons'][ll0:],self.attr['lons'][:ll0]+360))
        new.attr['lats'] = self.attr['lats']
        new.nlon = len(new.attr['lons'])
        new.nlat = len(new.attr['lats'])
        new.nlev = self.nlev
        new.attr['levtype'] = self.attr['levtype']
        new.attr['levs'] = self.attr['levs']
        new.attr['plev'] = self.attr['plev']
        new.date = self.date
        new.project = self.project
        new.attr['La1'] = self.attr['La1']
        new.attr['Lo1'] = 0
        try:
            new.attr['dla'] = self.attr['dla']
            new.attr['dlo'] = self.attr['dlo']
        except:
            new.attr['dla'] = (self.attr['lats'][-1]-self.attr['lats'][0])/(new.nlat-1)
            new.attr['dlo'] = (self.attr['lons'][-1]-self.attr['lons'][0])/(new.nlon-1)
        try:
            new.attr['pscale'] = self.attr['pscale']
            new.attr['zscale'] = self.attr['zscale']
            new.attr['zscale_i'] = self.attr['zscale_i']
        except: pass
        for var in self.var.keys():
            ndim = len(self.var[var].shape)
            new.var[var] = np.concatenate((self.var[var][...,ll0:],self.var[var][...,:ll0]),axis=ndim-1)
        return new

    def shift2west(self,lon0):
        # generate an object with the same content as self but with the origin
        # of longitude shifted to lon0
        # It is assumed that the first longitude is presently 0
        if self.attr['lons'][0] != 0:
            print('This operation can only be made on a grid strating at 0 longitude')
            return None
        new = ECMWF_pure()
        # find the location of lon0 in the initial grid
        # first make sure it is in the grid
        if lon0 < 0 : lon0 += 360
        ll0 = np.abs(self.attr['lons']-lon0).argmin()
        new.attr['lons'] = np.concatenate((self.attr['lons'][ll0:]-360,self.attr['lons'][:ll0]))
        new.attr['lats'] = self.attr['lats']
        new.nlon = len(new.attr['lons'])
        new.nlat = len(new.attr['lats'])
        new.nlev = self.nlev
        new.attr['levtype'] = self.attr['levtype']
        new.attr['levs'] = self.attr['levs']
        new.attr['plev'] = self.attr['plev']
        new.date = self.date
        new.project = self.project
        new.attr['La1'] = self.attr['La1']
        new.attr['Lo1'] = 0
        try:
            new.attr['dla'] = self.attr['dla']
            new.attr['dlo'] = self.attr['dlo']
        except:
            new.attr['dla'] = (self.attr['lats'][-1]-self.attr['lats'][0])/(new.nlat-1)
            new.attr['dlo'] = (self.attr['lons'][-1]-self.attr['lons'][0])/(new.nlon-1)
        try:
            new.attr['pscale'] = self.attr['pscale']
            new.attr['zscale'] = self.attr['zscale']
            new.attr['zscale_i'] = self.attr['zscale_i']
        except: pass
        for var in self.var.keys():
            ndim = len(self.var[var].shape)
            new.var[var] = np.concatenate((self.var[var][...,ll0:],self.var[var][...,:ll0]),axis=ndim-1)
        return new

    def extract(self,latRange=None,lonRange=None,varss=None,vard=None,levs=None,copy=False):
        """ Extract all variables on a reduced grid
           varss is a list of the 3D variables
           vard is a list of the 2D variables
           The reduced grid is defined by latRange and lonRange
           The 3D extracted fields can be extracted on a range of levels
        """
        # first determine the boundaries of the domain
        new = ECMWF_pure()
        new.date = self.date
        if (latRange == []) | (latRange == None):
            new.attr['lats'] = self.attr['lats']
            nlatmin = 0
            nlatmax = len(self.attr['lats'])
        else:
            nlatmin = np.argmax(self.attr['lats']>latRange[0])-1
            nlatmax = np.argmax(self.attr['lats']>latRange[1])
        if (lonRange == []) | (lonRange == None):
            new.attr['lons'] = self.attr['lons']
            nlonmin = 0
            nlonmax = len(self.attr['lons'])
        else:
            nlonmin = np.argmax(self.attr['lons']>lonRange[0])-1
            nlonmax = np.argmax(self.attr['lons']>=lonRange[1])+1
        if levs == None:
            kup = 0
            kbot = self.nlev
        else:
            kup = max(levs[0],0)
            kbot = min(levs[1],self.nlev)
        new.attr['lats'] = self.attr['lats'][nlatmin:nlatmax]
        new.attr['lons'] = self.attr['lons'][nlonmin:nlonmax]
        new.nlat = len(new.attr['lats'])
        new.nlon = len(new.attr['lons'])
        new.nlev = kbot - kup
        new.attr['levtype'] = self.attr['levtype']
        new.attr['levs'] = self.attr['levs'][kup:kbot]
        try: new.attr['plev'] = self.attr['plev'][kup:kbot]
        except: pass
        new.project = self.project
        new.date = self.date
        new.attr['La1'] = self.attr['lats'][nlatmin]
        new.attr['Lo1'] = self.attr['lons'][nlonmin]
        try:
            new.attr['dla'] = self.attr['dla']
            new.attr['dlo'] = self.attr['dlo']
        except:
            new.attr['dla'] = (self.attr['lats'][-1]-self.attr['lats'][0])/(new.nlat-1)
            new.attr['dlo'] = (self.attr['lons'][-1]-self.attr['lons'][0])/(new.nlon-1)
        try:
            new.attr['pscale'] = self.attr['pscale']
            new.attr['zscale'] = self.attr['zscale']
            new.attr['zscale_i'] = self.attr['zscale_i']
        except: pass
        new.globalGrid = False
        # extraction
        if varss is None:
            list_vars = []
        elif varss == 'All':
            list_vars = self.var.keys()
        else:
            list_vars = varss
        for var in list_vars:
            if len(self.var[var].shape) == 3:
                new.var[var] = self.var[var][kup:kbot,nlatmin:nlatmax,nlonmin:nlonmax]
            else:
                new.var[var] = self.var[var][nlatmin:nlatmax,nlonmin:nlonmax]
        if vard is None:
            list_vars = []
        elif vard == 'All':
            list_vars = self.d2d.keys()
        else:
            list_vars = vard
        for var in list_vars:
            new.d2d[var] = self.d2d[var][nlatmin:nlatmax,nlonmin:nlonmax]
        # make a copy if the source is not intended to be kept
        if copy:
            for var in new.var.keys():
                new.var[var] = new.var[var].copy()
            for var in new.d2d.keys():
                new.d2d[var] = new.d2d[var].copy()
            for var in new.attr.keys():
                try: new.attr[var] = new.attr[var].copy()
                except AttributeError: pass
        return new

    def zonal(self,vars=None,vard=None):
        """ Calculate the zonal averages of the variables contained in self"""
        new = ECMWF_pure()
        new.attr['lats'] = self.attr['lats']
        new.nlat = len(new.attr['lats'])
        new.nlev = self.nlev
        if vars is None:
            list_vars = []
        elif vars == 'All':
            list_vars = self.var.keys()
        else:
            list_vars = vars
        for var in list_vars:
            lx = len(self.var[var].shape)-1
            new.var[var] = np.mean(self.var[var],axis=lx)
        if vard is None:
            list_vars = []
        elif vard == 'All':
            list_vars = self.d2d.keys()
        else:
            list_vars = vard
        for var in list_vars:
            new.d2d[var] = np.mean(self.d2d[var],axis=1)
        return new

    def getxy(self,var,lev,y,x):
        """ get the interpolated value of var in x, y on the level lev """
        # Quick n' Dirty version by nearest neighbour
        jy = np.abs(self.attr['lats']-y).argmin()
        ix = np.abs(self.attr['lons']-x).argmin()
        return self.var[var][lev,jy,ix]

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
        try:
            for var in self.d2d.keys():
               data.d2d[var] = cf1*self.d2d[var] + cf2*other.d2d[var]
        except:
            pass
        return data

    def interpolP(self,p,varList='All',latRange=None,lonRange=None):
        """ interpolate the variables to a pressure level or a set of pressure levels
            vars must be a list of variables or a single varibale
            p must be a list or np.array of pressures in Pascal
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
            z must be a list of altitudes inm
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

    def interpol_part(self,p,x,y,varList='All'):
        """ Interpolate the variables to the location of particles given by [p,y,x] using trilinear method."""
        if 'P' not in self.var.keys():
            self._mkp()
        if varList == 'All':
            varList = list(self.var.keys())
            varList.remove('SP')
            varList.remove('P')
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
            self.fhyb,void = tohyb()
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
            vhigh = (1-px)*(1-py)*self.var[var][lhyb,jy,ix] + (1-px)*py*self.var[var][lhyb,jy+1,ix] \
                  + px*(1-py)*self.var[var][lhyb,jy,ix+1] + px*py*self.var[var][lhyb,jy+1,ix+1]
            vlow  = (1-px)*(1-py)*self.var[var][lhyb+1,jy,ix] + (1-px)*py*self.var[var][lhyb+1,jy+1,ix] \
                  + px*(1-py)*self.var[var][lhyb+1,jy,ix+1] + px*py*self.var[var][lhyb+1,jy+1,ix+1]
            hc = hyb - lhyb
            result[var] = (1-hc)*vhigh + hc*vlow
            result[var][non_valid] = MISSING
        return result

    def _CPT(self):
        """ Calculate the cold point tropopause """
        if not set(['T','P']).issubset(self.var.keys()):
            print('T or P undefined')
            return
        levbnd = {'FULL-EA':[30,90],'FULL-EI':[15,43],'STC':[10,90]}
        # TODO: find a way to avoid the big loop
        self.d2d['pcold'] = np.empty(shape=(self.nlat,self.nlon))
        self.d2d['Tcold'] = np.empty(shape=(self.nlat,self.nlon))
        if 'Z' in self.var.keys():
            self.d2d['zcold'] = np.empty(shape=(self.nlat,self.nlon))
        # Calculate the cold point in the discrete profile
        # TO DO: make a smoother version with vertical interpolation
        nc = np.argmin(self.var['T'][levbnd[self.project][0]:levbnd[self.project][1],...],axis=0)\
           + levbnd[self.project][0]
        for jy in range(self.nlat):
            for ix in range(self.nlon):
                self.d2d['pcold'][jy,ix] = self.var['P'][nc[jy,ix],jy,ix]
                self.d2d['Tcold'][jy,ix] = self.var['T'][nc[jy,ix],jy,ix]
        if  'Z' in self.var.keys():
            for jy in range(self.nlat):
                for ix in range(self.nlon):
                    self.d2d['zcold'][jy,ix] = self.var['Z'][nc[jy,ix],jy,ix]
        return

    def _lzrh(self):
        """ Calculate the clear sky LZRH. Translated from LzrnN.m
        Not yet validated. Use with caution.
        Add the calculation of the all sky LZRH
        """
        if not set(['P','ASSWR','ASLWR','CSSWR','CSLWR']).issubset(self.var.keys()):
            print('P or heating rate missing')
            return
        # bound for the test of positive heating in the stratosphere
        levbnd = {'FULL-EA':31,'FULL-EI':15,'STC':11}
        # Restrict the column (top at 60 hPa)
        ntop1 = 60 - self.attr['levs'][0]
        nbot1 = self.attr['levs'][-1] - 30
        self.d2d['plzrh'] = np.ma.empty(shape=(self.nlat,self.nlon))
        self.d2d['ptlzrh'] = np.ma.empty(shape=(self.nlat,self.nlon))
        self.d2d['aslzrh'] = np.ma.empty(shape=(self.nlat,self.nlon))
        uniq = np.empty(shape=(self.nlat,self.nlon))
        for jy in range(self.nlat):
            for ix in range(self.nlon):
                #print(jy,ix)
                # Clear sky heating in the column
                csh = self.var['CSSWR'][ntop1:nbot1,jy,ix] + self.var['CSLWR'][ntop1:nbot1,jy,ix]
                # Define location of positive heating
                lp = (csh>0).astype(np.int64)
                # Find position of vertical crossings of zero heating
                pos = np.where(lp[1:]-lp[:-1]==-1)[0]
                # Number of crossings
                uniq[jy,ix] = len(pos)
                # Select column with at least one crossing and make sure there is heating
                # in the strato higher than 66 hpa level to eliminate tropopause folds in the subtropics
                #
                if (len(pos)==0) | (np.max(csh[:levbnd[self.project]])<0) :
                    self.d2d['plzrh'][jy,ix] = np.ma.masked
                    self.d2d['ptlzrh'][jy,ix] = np.ma.masked
                    self.d2d['aslzrh'][jy,ix] = np.ma.masked
                    continue
                # Calculate pressure at each crossing
                pzk = []
                px = []
                for k in range(len(pos)):
                    id = pos[k]
                    #@@ Temporary test
                    if (csh[id]*csh[id+1]>0) | (csh[id]<0):
                        print('ERROR: localization')
                        return
                    px.append(csh[id]/(csh[id]-csh[id+1]))
                    # Interopolate to get intersection
                    pzk.append(math.exp(px[k]*math.log(self.var['P'][id+1,jy,ix])+(1-px[k])*math.log(self.var['P'][id,jy,ix])))
                # Now we test the best continuity for pressure
                # Collect neighbour values of the LZRH pressure already calculated
                # and find intersections in the current column which is closest to these values

                if len(pos)==1:
                    k=0
                else:
                    if jy>0 & ix>0:
                        if ix<self.nlon-1:
                            pzt = list(self.d2d['plzrh'][jy-1,ix-1:ix+2])
                            pzt.append(self.d2d['plzrh'][jy,ix-1])
                        else:
                            pzt = list(self.d2d['plzrh'][jy-1,ix-1:ix+1])
                            pzt.append(self.d2d['plzrh'][jy,ix-1])
                    elif jy==0 & ix>0:
                        pzt = [self.d2d['plzrh'][jy,ix-1],]
                    elif jy>0 & ix==0:
                        pzt = list(self.d2d['plzrh'][jy-1,ix:ix+2])
                    else:
                        pzt = []
                    if len(pzt)>0:
                        avpzt = sum(pzt)/len(pzt)
                        k = np.argmin((np.array(pzk)-avpzt)**2)
                        # This patch corrects some rare artefacts near the tropopause
                        # (empirical)
                        if pzk[k] < 9500:
                            print(jy,ix,len(pzk),avpzt,len(pzt),pzt,type(pzt))
                            idx = np.argsort((np.array(pzk)-avpzt)**2)
                            if pzk[idx[1]]>12000:
                                k=idx[1]
                    else:
                        k=0
                try:
                    self.d2d['plzrh'][jy,ix] = pzk[k]
                except:
                    print('ERROR')
                    print(self.d2d['plzrh'].shape,jy,ix)
                    print(len(pzk),k)
                # test
                if jy==166 & ix==0:
                    print(k,pzk[k])
                # Calculate potential temperature and all sky heating at the LZRH
                tz = px[k]*self.var['T'][pos[k]+1,jy,ix] + (1-px[k])*self.var['T'][pos[k],jy,ix]
                self.d2d['ptlzrh'][jy,ix] = tz * (cst.p0/pzk[k])**cst.kappa
                self.d2d['aslzrh'][jy,ix] = px[k]*self.var['ASSWR'][pos[k]+1,jy,ix] + (1-px[k])*self.var['ASSWR'][pos[k],jy,ix] \
                                          + px[k]*self.var['ASLWR'][pos[k]+1,jy,ix] + (1-px[k])*self.var['ASLWR'][pos[k],jy,ix]
        return

    #@jit
    def _WMO(self,highlatOffset=False):
        """ Calculate the WMO tropopause
        When highlatoffset is true the 2K/km criterion is replaced by a 3K/km
        at high latitudes latitudes above 60S or 60N """
        if not set(['T','P']).issubset(self.var.keys()):
            print('T or P undefined')
            return
        self.d2d['pwmo'] = np.ma.empty(shape=(self.nlat,self.nlon))
        self.d2d['Twmo'] = np.ma.empty(shape=(self.nlat,self.nlon))
        if 'Z' in self.var.keys():
            self.d2d['zwmo'] = np.empty(shape=(self.nlat,self.nlon))
            zwmo = True
        else:
            zwmo = False
        levbnd = {'FULL-EA':[30,90],'FULL-EI':[15,43],'STC':[10,85]}
        highbnd = levbnd[self.project][0]
        lowbnd =  levbnd[self.project][1]
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
                        self.d2d['pwmo'][jy,ix] = np.ma.masked
                        self.d2d['Twmo'][jy,ix] = np.ma.masked
                        if zwmo: self.d2d['zwmo'][jy,ix] = np.ma.masked
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
                        self.d2d['pwmo'][jy,ix] = self.var['P'][lev0,jy,ix]
                        self.d2d['Twmo'][jy,ix] = self.var['T'][lev0,jy,ix]
                        if zwmo: self.d2d['zwmo'][jy,ix] = self.var['Z'][lev0,jy,ix]
        return

    def interpol_track(self,p,x,y,varList='All'):
        """ Writing in progress. For the moment, this is a copy of interpol_part.
        Calculate the distance to the cold point and to the LZRH .
        We use RegularGridInterpolator that allows interpolation to a list of
        points unlike interp2d."""
        if 'P' not in self.var.keys():
            self._mkp()
        if varList == 'All':
            varList = list(self.var.keys())
            varList.remove('SP')
            varList.remove('P')
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
            self.fhyb,void = tohyb()
        # generate the 2D interpolator of the surface log pressure
        lsp = RegularGridInterpolator((self.attr['lats'],self.attr['lons']),-np.log(self.var['SP']),bounds_error=True)
        lspi = lsp(np.transpose([y,x]))
        # define -log sigma = -log(p) - -log(ps)
        lsig = - np.log(p) - lspi
        # find the non integer hybrib level
        hyb = self.fhyb(np.transpose([lsig,lspi]))-self.attr['levs'][0]+1
        lhyb = np.floor(hyb).astype(np.int64)
        #@@ test the extreme values of sigma end ps
        if np.min(lsig) < - np.log(0.95):
            print('large sigma detected ',np.exp(-np.min(lsig)))
        if np.max(lspi) > -np.log(45000):
            print('small ps detected ',np.exp(-np.max(lspi)))
        # Horizontal interpolation
        ix = np.floor((x-self.attr['Lo1'])/self.attr['dlo']).astype(np.int64)
        jy = np.floor((y-self.attr['La1'])/self.attr['dla']).astype(np.int64)
        px = ((x - self.attr['Lo1']) %  self.attr['dlo'])/self.attr['dlo']
        py = ((y - self.attr['La1']) %  self.attr['dla'])/self.attr['dla']
        for var in varList:
            vhigh = (1-px)*(1-py)*self.var[var][lhyb,jy,ix] + (1-px)*py*self.var[var][lhyb,jy+1,ix] \
                  + px*(1-py)*self.var[var][lhyb,jy,ix+1] + px*py*self.var[var][lhyb,jy+1,ix+1]
            vlow  = (1-px)*(1-py)*self.var[var][lhyb+1,jy,ix] + (1-px)*py*self.var[var][lhyb+1,jy+1,ix] \
                  + px*(1-py)*self.var[var][lhyb+1,jy,ix+1] + px*py*self.var[var][lhyb+1,jy+1,ix+1]
            hc = hyb % 1
            result[var] = (1-hc)*vhigh + hc*vlow
            del hc; del vhigh; del vlow
        return result

    def interpol_orbit(self,x,y,varList='All',var2='All'):
        """ Generate an interpolation to an orbit curtain in 2d.
        Quick n' dirty method  using floor. Should be improved for better accuracy"""
        if varList == 'All':
            varList = list(self.var.keys())
        if var2=='All':
            var2 = list(self.d2d.keys())
        sect = curtain()
        sect.attr['lons'] = x
        sect.attr['lats'] = y
        sect.n = len(x)
        # extract within a bounding box
        #bb = self.extract(latRange=(y.min(),y.max()),lonRange=(x.min(),x.max()),vard=['SP',],varss=varList)
        # Horizontal interpolation
        ix = np.floor((x-self.attr['Lo1'])/self.attr['dlo']).astype(np.int)
        # trick to handle periodicity in longitude and profiles in the cutting strip
        ix1 = (ix+1) % self.nlon
        jy = np.floor((y-self.attr['La1'])/self.attr['dla']).astype(np.int)
        px = ((x - self.attr['Lo1']) %  self.attr['dlo'])/self.attr['dlo']
        py = ((y - self.attr['La1']) %  self.attr['dla'])/self.attr['dla']
        for var in var2:
            sect.var[var] = (1-px)*(1-py)*self.d2d[var][jy,ix] + (1-px)*py*self.d2d[var][jy+1,ix] \
                  + px*(1-py)*self.d2d[var][jy,ix1] + px*py*self.d2d[var][jy+1,ix1]
        for var in varList:
            if len(self.var[var].shape)==2:
                sect.var[var] = (1-px)*(1-py)*self.var[var][jy,ix] + (1-px)*py*self.var[var][jy+1,ix] \
                    + px*(1-py)*self.var[var][jy,ix1] + px*py*self.var[var][jy+1,ix1]
            else:
                sect.var[var] = np.empty(shape=(self.nlev,len(x)))
                for lev in range(self.nlev):
                    sect.var[var][lev,:] = (1-px)*(1-py)*self.var[var][lev,jy,ix] + (1-px)*py*self.var[var][lev,jy+1,ix] \
                        + px*(1-py)*self.var[var][lev,jy,ix1] + px*py*self.var[var][lev,jy+1,ix1]
        return sect

# standard class to read data
class ECMWF(ECMWF_pure):
    # to do: raise exception in case of an error

    def __init__(self,project,date,step=0,exp=[None]):
        ECMWF_pure.__init__(self)
        self.project = project
        self.date = date
        self.exp = exp
        #SP_expected = False
        self.EN_expected = False
        self.DI_expected = False
        self.WT_expected = False
        self.VD_expected = False
        self.DE_expected = False
        self.x4I_expected = False
        self.VOZ_expected = False
        self.QN_expected = False
        self.F12_expected = False
        self.CAMS_expected = False
        self.CF12_expected = False
        # set when some auxilliary files are hemispheric (one choice for all concerned)
        self.hemis = None
        self.offd = 100 # offset to be set for ERA-I below, 100 for ERA5
        if self.project=='VOLC':
            if 'satie' in socket.gethostname():
                self.rootdir = '/dsk2/ERA5/VOLC'
            elif 'ens' in socket.gethostname():
                self.rootdir = '/net/satie/dsk2/ERA5/VOLC'
            else:
                print('unknown hostname for this dataset')
                return
            self.globalGrid = False
            self.EN_expected = True
            self.DI_expected = True
            self.WT_expected = True

        elif project=='STC':
            if 'gort' == socket.gethostname():
                self.rootdir = '/dkol/dc6/samba/STC/ERA5/STC'
            elif 'ciclad' in socket.gethostname():
                self.rootdir = '/proju/flexpart/flexpart_in/STC/ERA5'
            elif 'climserv' in socket.gethostname():
                self.rootdir = '/proju/flexpart/flexpart_in/STC/ERA5'
            elif 'camelot' in socket.gethostname():
                self.rootdir = '/proju/flexpart/flexpart_in/STC/ERA5'
            elif 'satie' in socket.gethostname():
                self.rootdir = '/dsk2/ERA5/STC'
            elif socket.gethostname() in ['grapelli','coltrane','zappa','couperin','puccini','lalo']:
                self.rootdir = '/net/satie/dsk2/ERA5/STC'
            else:
                print('unknown hostname for this dataset')
                return
            self.globalGrid = False
            self.EN_expected = True
            self.DI_expected = True
            self.WT_expected = True
            self.VD_expected = True

        elif project=='FULL-EI':
            if 'gort' == socket.gethostname():
                self.rootdir = '/dkol/data/NIRgrid'
            elif 'ciclad' in socket.gethostname():
                self.rootdir = '/proju/flexpart/flexpart_in/NIRgrid'
            elif 'climserv' in socket.gethostname():
                self.rootdir = '/proju/flexpart/flexpart_in/NIRgrid'
            elif 'camelot' in socket.gethostname():
                self.rootdir = '/proju/flexpart/flexpart_in/NIRgrid'
            elif 'satie' in socket.gethostname():
                self.rootdir = '/data/NIRgrid'
            else:
                print('unknown hostname for this dataset')
                return
            self.globalGrid = True
            self.EN_expected = True
            self.DI_expected = True
            self.DE_expected = True
            self.offd = 300

        elif project=='FULL-EA':
            if 'gort' == socket.gethostname():
                self.rootdir = '/dkol/data/ERA5'
            elif 'ciclad' in socket.gethostname():
                self.rootdir = '/data/legras/flexpart_in/ERA5'
            elif 'climserv' in socket.gethostname():
                self.rootdir = '/data/legras/flexpart_in/ERA5'
            elif 'satie' in socket.gethostname():
                self.rootdir = '/data/ERA5'
            elif 'Graphium' == socket.gethostname():
                self.rootdir = 'C:\\cygwin64\\home\\berna\\data\\ERA5'
            else:
                print('unknown hostname for this dataset')
                return
            self.globalGrid = True
            self.EN_expected = True
            if (self.exp == 'VOZ') | ('VOZ' in self.exp): self.VOZ_expected = True
            if (self.exp == 'DI') | ('DI' in self.exp): self.DI_expected = True
            if (self.exp == 'QN') | ('QN' in self.exp): self.QN_expected = True
            if (self.exp == 'F12') | ('F12' in self.exp):
                if self.date.hour not in [6,18]:
                    print('non valid time for the 12h forecast')
                    return
                self.F12_expected = True

        elif project=='OPZ':
            if 'gort' == socket.gethostname():
                self.rootdir = '/dkol/data/OPZ'
            elif 'ciclad' in socket.gethostname():
                self.rootdir = '/data/legras/flexpart_in/OPZ'
            elif 'climserv' in socket.gethostname():
                self.rootdir = '/data/legras/flexpart_in/OPZ'
            elif 'satie' in socket.gethostname():
                self.rootdir = '/data/OPZ'
            else:
                print('unknown hostname for this dataset')
                return
            self.globalGrid = True
            self.EN_expected = True
            if (self.exp == 'x4I') | ('x4I' in self.exp): self.x4I_expected = True
            if (self.exp == 'CAMS') | ('CAMS' in self.exp): self.CAMS_expected = True
            if (self.exp == 'F12') | ('F12' in self.exp):
                if self.date.hour not in [6,18]:
                    print('non valid time for the 12h forecast')
                    return
                self.F12_expected = True
            if (self.exp == 'CF12') | ('CF12' in self.exp):
                if self.date.hour not in [0,12]:
                    print('non valid time for the 12h forecast')
                    return
                self.CF12_expected = True

        elif project=='OPZFCST':
            if 'gort' == socket.gethostname():
                self.rootdir = '/dkol/data/OPZ'
            elif 'ciclad' in socket.gethostname():
                self.rootdir = '/data/legras/flexpart_in/OPZ'
            elif 'climserv' in socket.gethostname():
                self.rootdir = '/data/legras/flexpart_in/OPZ'
            elif 'satie' in socket.gethostname():
                self.rootdir = '/data/OPZ'
            else:
                print('unknown hostname for this dataset')
                return
            self.globalGrid = True
            self.EN_expected = True
        else:
            print('Non implemented project')
            return

        #if SP_expected:
        #    self.fname = 'SP'+date.strftime('%y%m%d%H')
        #    path1 = 'SP-true/grib'
        if self.EN_expected:
            if project in ['STC','VOLC']:
                self.fname = date.strftime('EN%y%m%d%H')
                path1 = 'EN-true/grib'
                self.ENvar = {'U':['u','U component of wind','m s**-1'],
                     'V':['v','V component of wind','m s**-1'],
                     'W':['w','Vertical velocity','Pa s**-1'],
                     'T':['t','Temperature','K'],
                     'LNSP':['lnsp','Logarithm of surface pressure','Log(Pa)']
                     }
            elif project == 'FULL-EA':
                self.fname = date.strftime('ERA5EN%Y%m%d.grb')
                path1 ='EN-true'
                self.ENvar = {'U':['u','U component of wind','m s**-1'],
                     'V':['v','V component of wind','m s**-1'],
                     'W':['w','Vertical velocity','Pa s**-1'],
                     'T':['t','Temperature','K'],
                     'LNSP':['lnsp','Logarithm of surface pressure','Log(Pa)']}
            elif project == 'OPZ':
                self.fname = date.strftime('OPZLWDA%Y%m%d.grb')
                path1 ='EN-true'
                self.ENvar = {'U':['u','U component of wind','m s**-1'],
                     'V':['v','V component of wind','m s**-1'],
                     'W':['w','Vertical velocity','Pa s**-1'],
                     'T':['t','Temperature','K'],
                     'LNSP':['lnsp','Logarithm of surface pressure','Log(Pa)'],
                     'VO':['vo','Vorticity','s**-1'],
                     'D':['d','Divergence','s**-1'],
                     'Q':['q','Specific humidity','kg kg**-1'],
                     'O3':['o3','Ozone mixing ratio','kg kg**-1'],
                     'WE':['etadot','Eta-coordinate vertical velocity','s**-1']
                     }
            elif project == 'OPZFCST':
                self.fname = date.strftime('OPZFCST%Y%m%d_SH.grb')
                path1 ='EN-true'
                self.ENvar = {'T':['t','Temperature','K'],
                     'LNSP':['lnsp','Logarithm of surface pressure','Log(Pa)'],
                     'VO':['vo','Vorticity','s**-1'],
                     'O3':['o3','Ozone mixing ratio','kg kg**-1']}
            elif project == 'FULL-EI':
                self.fname = date.strftime('EN%y%m%d%H')
                path1 = 'EN-true'
                self.ENvar = {'U':['u','U component of wind','m s**-1'],
                     'V':['v','V component of wind','m s**-1'],
                     'W':['w','Vertical velocity','Pa s**-1'],
                     'T':['t','Temperature','K'],
                     'SP':['sp','Surface pressure','Pa'],
                     'U10':['10u','10 metre U wind component','m s**-1'],
                     'V10':['10v','10 metre V wind component','m s**-1'],
                     'T2':['2t','2 metre temperature','K'],
                     'TD2':['2d','2 metre dewpoint temperature','K']
                     }
                # A number of other surface variables are available

        if self.DI_expected:
            if project == 'FULL-EI':
                # for ERA-I: tendencies over 3-hour intervals following file date, provided as temperature increments
                self.DIvar = {'ASSWR':['srta','Mean temperature tendency due to short-wave radiation','K'],
                     'ASLWR':['trta','Mean temperature tendency due to long-wave radiation','K'],
                     'CSSWR':['srtca','Mean temperature tendency due to short-wave radiation, clear sky','K'],
                     'CSLWR':['trtca','Mean temperature tendency due to long-wave radiation, clear sky','K'],
                     'PHR':['ttpha','Mean temperature tendency due to physics','K']}
                self.dname = 'DI'+date.strftime('%y%m%d%H')
            elif project == 'FULL-EA':
                # for ERA5: tendencies over 1-hour intervals following file date, provided as temperature increments
                # notice that native names differ from that of ERA-I
                self.DIvar = {'ASSWR':['mttswr','Mean temperature tendency due to short-wave radiation','K s**-1'],
                     'ASLWR':['mttlwr','Mean temperature tendency due to long-wave radiation','K s**-1'],
                     'CSSWR':['mttswrcs','Mean temperature tendency due to short-wave radiation, clear sky','K s**-1'],
                     'CSLWR':['mttlwrcs','Mean temperature tendency due to long-wave radiation, clear sky','K s**-1'],
                     'PHR':['mttpm','Mean temperature tendency due to parametrerizations','K s**-1'],}
                self.dname = date.strftime('ERA5DI%Y%m%d.grb')
            else:
                # for ERA5: tendencies over 1-hour intervals following file date
                self.DIvar = {'ASSWR':['mttswr','Mean temperature tendency due to short-wave radiation','K s**-1'],
                     'ASLWR':['mttlwr','Mean temperature tendency due to long-wave radiation','K s**-1'],
                     'CSSWR':['mttswrcs','Mean temperature tendency due to short-wave radiation, clear sky','K s**-1'],
                     'CSLWR':['mttlwrcs','Mean temperature tendency due to long-wave radiation, clear sky','K s**-1'],
                     'PHR':['mttpm','Mean temperature tendency due to parametrisations','K s**-1'],
                     'UMF':['mumf','Mean updraught mass flux','kg m**-2 s**-1'],
                     'DMF':['mdmf','Mean downdraught mass flux','kg m**-2 s**-1'],
                     'UDR':['mudr','Mean updraught detrainment rate','kg m**-3 s**-1'],
                     'DDR':['mddr','Mean downdraught detrainement rate','kg m**-3 s**-1']}
                self.dname = 'DI'+date.strftime('%y%m%d%H')
        if self.WT_expected:
            # for ERA5
            self.WTvar = {'CRWC':['crwc','Specific rain water content','kg kg**-1'],
                     'CSWC':['cswc','Specific snow water content','kg kg**-1'],
                     'Q':['q','Specific humidity','kg kg**-1'],
                     'QL':['clwc','Specific cloud liquid water content','kg kg**-1'],
                     'QI':['ciwc','Specific cloud ice water content','kg kg**-1'],
                     'CC':['cc','Fraction of cloud cover','0-1']}
            self.wname = 'WT'+date.strftime('%y%m%d%H')
        if self.VD_expected:
            # for ERA5
            self.VDvar = {'DIV':['d','Divergence','s**-1'],
                          'VO':['vo','Vorticity','s**-1'],
                          #'Z':['xxx','Geopotential','m'],
                          'WE':['etadot','Eta dot','s**-1*']}
            self.vname = 'VD'+date.strftime('%y%m%d%H')
        if self.DE_expected:
            # for ERA-I
            self.DEvar = {'UDR':['udra','Mean updraught detrainment rate','kg m**-3 s**-1']}
            self.dename = 'DE'+date.strftime('%y%m%d%H')
        if self.x4I_expected:
            self.x4Ivar = {'VOD':['vodiff',"Vorticity increment","s**-1"],
                           'O3D':['o3diff','Ozone mixing ratio increment','kg kg**-1'],
                           'DD':['ddiff','Divergence increment','s**-1'],
                           'QD':['qdiff','Specidic humidity increment','kg kg**-1'],
                           'TD':['tdiff','Temperature increment','K'],
                           'LNSPD':['lnspdiff','Surface log pressure increment','']}
            self.x4iname = date.strftime('4ILWDA%Y%m%d.grb')
        if self.VOZ_expected:
            self.VOZvar = {'VO':['vo','Vorticity','s**-1'],
                           'O3':['o3','Ozone mixing ratio','kg kg**-1']}
            self.vozname = date.strftime('ERA5VOZ%Y%m%d.grb')
        if self.QN_expected:
            self.QNvar = {'Q':['q','Specific humidity','kg kg**-1'],
                          'QI':['ciwc','Specific cloud ice water content','kg kg**-1'],
                          'Ql':['clwc','Specific cloud ice water content','kg kg**-1'],
                          'QR':['crwc','Specific rain water content:kg kg**-1'],
                          'QS':['cswc','Specific snow water content:kg kg**-1']}
            self.qnname = date.strftime('ERA5QN%Y%m%d.grb')
        if self.F12_expected:
            self.F12var = {'UF':['u','12h forecast U component of wind','m s**-1'],
                     'VF':['v','12h forecast V component of wind','m s**-1'],
                     'TF':['t','12h forecast temperature','K'],
                     'VOF':['vo',"12h forecast vorticity","s**-1"],
                     'O3F':['o3',"12h forecast ozone mixing ratio",'kg kg**-1']}
             #self.f12name defined in the opening
        if self.CAMS_expected:
            self.CAMSvar = {'CO':['co','Carbon monoxide','kg kg**-1'],
                'O3G':['go3','GEMS Ozone','kg kg**-1'],
                'O3S':['o3s','Stratospheric ozone:kg kg**-1'],
                'O3C':['o3','Ozone mass mixing ratio','kg kg**-1']}
            self.camsname = date.strftime('OPZCAMS-%Y%m%d_SH.grb')
        if self.CF12_expected:
            self.CF12var = {'COF':['co','12h forecast carbon monoxide','kg kg**-1'],
                'O3GF':['go3','12h forecast GEMS Ozone','kg kg**-1'],
                'O3SF':['o3s','12h forecast stratospheric ozone:kg kg**-1'],
                'O3CF':['o3','12h ozone mass mixing ratio','kg kg**-1']}
            #self.cf12name defined in the opening
        # opening first the DI file as it might be needed to get pressure for hours not multiple of 3
        self.DI_open = False
        if self.DI_expected:
            try:
                self.drb = pygrib.open(os.path.join(self.rootdir,date.strftime('DI-true/grib/%Y/%m'),self.dname))
                self.DI_open = True
            except:
                try:
                    self.drb = pygrib.open(os.path.join(self.rootdir,date.strftime('DI-true/%Y'),self.dname))
                    self.DI_open = True
                except:
                    print('cannot open '+os.path.join(self.rootdir,date.strftime('DI-true/grib/%Y/%m'),self.dname))
        # opening first the CAMS file as it might be needed to read pressure at 0 or 12
        self.CAMS_open = False
        if self.CAMS_expected:
            try:
                self.crb = pygrib.open(os.path.join(self.rootdir,date.strftime('MC-true/%Y'),self.camsname))
                self.CAMS_open = True
                self.hemis = 'SH' # future: to be dynamically determined
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('MC-true/%Y'),self.camsname))
        # opening the main EN file
        try:
            self.grb = pygrib.open(os.path.join(self.rootdir,path1,date.strftime('%Y/%m'),self.fname))
        except:
            try:
                self.grb = pygrib.open(os.path.join(self.rootdir,path1,date.strftime('%Y'),self.fname))
            except:
                print('cannot open '+os.path.join(self.rootdir,path1,date.strftime('%Y/%m'),self.fname))
                # We do not need to open EN if we only want CAMS at 0 and 12 to calculate assimilation increment
                if self.project=='OPZ' & date.hour in [0,12]:
                    pass
                else:
                    return

        # Define searched valid date and time, and step
        validityDate = int(self.date.strftime('%Y%m%d'))
        validityTime = int(self.date.strftime('%H%M'))
        self.EN_open = True
        step
        try:
            sp = self.grb.select(name='Logarithm of surface pressure',validityTime=validityTime)[0]
            logp = True
        except:
            try:
                if self.project=='OPZFCST':
                    sp = self.grb.select(name='Logarithm of surface pressure',validityTime=validityTime,step=step)[0]
                    logp = True
                elif (self.project == 'FULL-EA') & self.DI_open:
                    sp = self.drb.select(name='Logarithm of surface pressure',validityTime=validityTime)[0]
                    logp = True
                elif (self.project == 'OPZ') & self.CAMS_open:
                    sp = self.crb.select(name='Logarithm of surface pressure',validityTime=validityTime)[0]
                    logp = True
                else:
                    sp = self.grb.select(name='Surface pressure',validityTime=validityTime)[0]
                    logp = False
            except:
                print('no surface pressure in '+self.fname)
                self.grb.close()
                return
        # Check date matching (should be OK)
        if (sp['validityDate'] != validityDate) & (self.project!='OPZFCST'):
            print('WARNING: dates do not match')
            print('called date    ',self.date.strftime('%Y%m%d %H%M'))
            print('date from file ',sp['validityDate'],sp['validityTime'])
            self.grb.close()
            return
        # Get general info from this message
        # ACHTUNG if CAMS OPZ at 0h and 12h the sp is at a hemispheric resolution,
        # if 6h and 18h, the sp is from EN file at full resolution
        self.attr['Date'] = sp['dataDate']
        self.attr['Time'] = sp['dataTime']
        self.attr['valDate'] = sp['validityDate']
        self.attr['valTime'] = sp['validityTime']
        self.attr['step'] = sp['step']
        self.nlon = sp['Ni']
        self.nlat = sp['Nj']
        self.attr['lons'] = sp['distinctLongitudes']
        self.attr['lats'] = sp['distinctLatitudes']
        # in ERA-Interim, it was necessary to divide by 1000
        if project == 'FULL-EI':
          self.attr['Lo1'] = sp['longitudeOfFirstGridPoint']/1000  # in degree
          self.attr['Lo2'] = sp['longitudeOfLastGridPoint']/1000  # in degree
          # reversing last and first for latitudes
          self.attr['La1'] = sp['latitudeOfLastGridPoint']/1000  # in degree
          self.attr['La2'] = sp['latitudeOfFirstGridPoint']/1000 # in degree
        else:
            self.attr['Lo1'] = sp['longitudeOfFirstGridPoint']/1000000  # in degree
            self.attr['Lo2'] = sp['longitudeOfLastGridPoint']/1000000  # in degree
            # reversing last and first for latitudes
            self.attr['La1'] = sp['latitudeOfLastGridPoint']/1000000  # in degree
            self.attr['La2'] = sp['latitudeOfFirstGridPoint']/1000000 # in degree
        # longitude and latitude interval
        self.attr['dlo'] = (self.attr['lons'][-1] - self.attr['lons'][0]) / (self.nlon-1)
        self.attr['dla'] = (self.attr['lats'][-1] - self.attr['lats'][0]) / (self.nlat-1)
        if self.attr['Lo1']>self.attr['Lo2']:
            self.attr['Lo1'] = self.attr['Lo1']-360.
        if sp['PVPresent']==1 :
            pv = sp['pv']
            self.attr['ai'] = pv[0:int(pv.size/2)]
            self.attr['bi'] = pv[int(pv.size/2):]
            self.attr['am'] = (self.attr['ai'][1:] + self.attr['ai'][0:-1])/2
            self.attr['bm'] = (self.attr['bi'][1:] + self.attr['bi'][0:-1])/2
        else:
            # todo: provide a fix, eg for JRA-55
            print('missing PV not implemented')
        self.attr['levtype'] = 'hybrid'
        # Set the surface pressure from sp field
        if logp:
            self.var['SP'] = np.exp(sp['values'])
        else:
            self.var['SP'] = sp['values']
        # Reverting lat order
        self.attr['lats'] = self.attr['lats'][::-1]
        self.var['SP']   = self.var['SP'][::-1,:]
        self.attr['dla'] = -  self.attr['dla']
        # Set the hemis parameter for OPZ in the northern hemisphere
        if (self.project == 'OPZ') & (self.attr['La1']==0):
            print('set NH')
            self.hemis = 'NH'
        # Opening of the other files
        self.WT_open = False
        self.VD_open = False
        self.DE_open = False
        self.x4I_open = False
        self.VOZ_open = False
        self.QN_open = False
        self.F12_open = False
        self.CF12_open = False
        if self.DI_expected & ~self.DI_open:
            try:
                self.drb = pygrib.open(os.path.join(self.rootdir,date.strftime('DI-true/grib/%Y/%m'),self.dname))
                self.DI_open = True
            except:
                try:
                    self.drb = pygrib.open(os.path.join(self.rootdir,date.strftime('DI-true/%Y'),self.dname))
                    self.DI_open = True
                except:
                    print('cannot open '+os.path.join(self.rootdir,date.strftime('DI-true/grib/%Y/%m'),self.dname))
        if self.CAMS_expected & ~self.CAMS_open:
            try:
                self.crb = pygrib.open(os.path.join(self.rootdir,date.strftime('MC-true/%Y'),self.camsname))
                self.CAMS_open = True
                self.hemis = 'SH' # future: to be dynamically determined
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('MC-true/%Y'),self.camsname))
        if self.WT_expected:
            try:
                self.wrb = pygrib.open(os.path.join(self.rootdir,date.strftime('WT-true/grib/%Y/%m'),self.wname))
                self.WT_open = True
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('WT-true/grib/%Y/%m'),self.wname))
        if self.VD_expected:
            try:
                self.vrb = pygrib.open(os.path.join(self.rootdir,date.strftime('VD-true/grib/%Y/%m'),self.vname))
                self.VD_open = True
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('VD-true/grib/%Y/%m'),self.vname))
        if self.DE_expected:
            try:
                self.derb = pygrib.open(os.path.join(self.rootdir,date.strftime('DE-true/%Y'),self.dename))
                self.DE_open = True
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('DE-true/%Y'),self.dename))
        if self.x4I_expected:
            try:
                self.x4Irb = pygrib.open(os.path.join(self.rootdir,date.strftime('EN-true/%Y'),self.x4iname))
                self.x4I_open = True
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('EN-true/%Y'),self.x4iname))
        if self.VOZ_expected:
            try:
                self.vozrb = pygrib.open(os.path.join(self.rootdir,date.strftime('VO3-true/%Y'),self.vozname))
                self.VOZ_open = True
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('VO3-true/%Y'),self.vozname))
        if self.QN_expected:
            try:
                self.qnrb = pygrib.open(os.path.join(self.rootdir,date.strftime('QN-true/%Y'),self.qnname))
                self.QN_open = True
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('QN-true/%Y'),self.qnname))
        if self.F12_expected:
            if self.date.hour == 6: datef = date - timedelta(days=1)
            else: datef = date
            if project == 'OPZ':
                self.f12name = datef.strftime('OPZLWDA-FCST12-%Y%m%d_SH.grb')
            else:
                self.f12name = datef.strftime('ERA5FCST12%Y%m%d.grb')
            try:
                if project == 'OPZ':
                    self.f12rb = pygrib.open(os.path.join(self.rootdir,datef.strftime('EN-true/%Y'),self.f12name))
                    # contrary to initial intention these data are on full grid
                    #self.hemis = 'SH'
                else:
                    self.f12rb = pygrib.open(os.path.join(self.rootdir,datef.strftime('FCST12-true/%Y'),self.f12name))
                self.F12_open = True
                sp = self.f12rb.select(name='Logarithm of surface pressure',validityTime=self.attr['valTime'])[0]
                self.var['SPF'] = np.exp(sp['values'])
                self.var['SPF']   = self.var['SPF'][::-1,:]
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('FCST12-true/%Y'),self.f12name))
        if self.CF12_expected:
            if self.date.hour == 0: datef = date - timedelta(days=1)
            else: datef = date
            self.cf12name = datef.strftime('OPZCAMS-FCST12-%Y%m%d_SH.grb')
            try:
                self.cf12rb = pygrib.open(os.path.join(self.rootdir,datef.strftime('MC-true/%Y'),self.cf12name))
                self.CF12_open = True
                sp = self.cf12rb.select(name='Logarithm of surface pressure',validityTime=self.attr['valTime'])[0]
                self.var['SPF'] = np.exp(sp['values'])
                self.var['SPF']   = self.var['SPF'][::-1,:]
            except:
                print('cannot open '+os.path.join(self.rootdir,date.strftime('FCST12-true/%Y'),self.f12name))

    def close(self):
        self.grb.close()
        try: self.drb.close()
        except: pass
        try: self.wrb.close()
        except: pass
        try: self.derb.close()
        except: pass
        try: self.vrb.close()
        except: pass

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
    def _get_var(self,var,step=None):
        if (var in self.var.keys()) & (self.project != 'OPZFCST'):
            return
        get = False
        try:
            if self.EN_open:
                if var in self.ENvar.keys():
                    if self.project=='OPZFCST':
                        if step is not None:
                            stepi = step
                            self.attr['step'] = step
                        else: stepi = self.attr['step']
                        TT = self.grb.select(shortName=self.ENvar[var][0],validityTime=self.attr['valTime'],step=stepi)
                    else:
                        TT = self.grb.select(shortName=self.ENvar[var][0],validityTime=self.attr['valTime'])
                    get = True
            if ~get & self.DI_open:
                if var in self.DIvar.keys():
                    TT = self.drb.select(shortName=self.DIvar[var][0],validityTime=(self.attr['valTime']+self.offd) % 2400)
                    get = True
            if ~get & self.WT_open:
                if var in self.WTvar.keys():
                    TT = self.wrb.select(shortName=self.WTvar[var][0],validityTime=self.attr['valTime'])
                    get = True
            if ~get & self.VD_open:
                if var in self.VDvar.keys():
                    TT = self.vrb.select(shortName=self.VDvar[var][0],validityTime=self.attr['valTime'])
                    get = True
            if ~get & self.DE_open:
                if var in self.DEvar.keys():
                    # Shift of validity time due to ERA-I convention
                    TT = self.derb.select(shortName=self.DEvar[var][0],validityTime=(self.attr['valTime']+self.offd) % 2400)
                    get = True
            if ~get & self.x4I_open:
                if var in self.x4Ivar.keys():
                    # sum of the four partial increments
                    TT = {}
                    print('TT created')
                    for i in range(4):
                        TT[i] = self.x4Irb.select(shortName=self.x4Ivar[var][0],validityTime=self.attr['valTime']+300,iterationNumber=i)
                        print('exit',i)
                    get = True
            if ~get & self.VOZ_open:
                if var in self.VOZvar.keys():
                    TT = self.vozrb.select(shortName=self.VOZvar[var][0],validityTime=self.attr['valTime'])
                    get = True
            if ~get & self.QN_open:
                if var in self.QNvar.keys():
                    TT = self.qnrb.select(shortName=self.QNvar[var][0],validityTime=self.attr['valTime'])
                    get = True
            if ~get & self.F12_open:
                if var in self.F12var.keys():
                    TT = self.f12rb.select(shortName=self.F12var[var][0],validityTime=self.attr['valTime'])
                    get = True
            if ~get & self.CAMS_open:
                if var in self.CAMSvar.keys():
                    TT = self.crb.select(shortName=self.CAMSvar[var][0],validityTime=self.attr['valTime'])
                    get = True
            if ~get & self.CF12_open:
                if var in self.CF12var.keys():
                    TT = self.cf12rb.select(shortName=self.CF12var[var][0],validityTime=self.attr['valTime'])
                    get = True
            if get == False:
                    print(var+' not found')
                    return
        except:
            print(var+' not found or read error')
            return
        # Process each message corresponding to a level
        # Special case of the 4D var increment first
        # ugly programming
        x4ISpecial = False
        if self.x4I_open:
            if var in self.x4Ivar.keys(): x4ISpecial = True
        if x4ISpecial:
            self.var[var] = np.zeros(shape=[self.nlev,self.nlat,self.nlon])
            for i in range(4):
               for l in range(len(TT[i])):
                   self.var[var][l,:,:] += TT[i][l]['values']
            readlev = False
        # Common case
        else:
            self.TT = TT
            if 'levs' not in self.attr.keys():
                readlev=True
                self.nlev = len(TT)
                self.attr['levs'] = np.full(len(TT),MISSING,dtype=int)
                self.attr['plev'] = np.full(len(TT),MISSING,dtype=float)
            else:
                readlev=False
                if (len(TT) !=  self.nlev) & (var != 'LNSP'):
                    print('new record inconsistent with previous ones ',var,len(TT))
            self.var[var] = np.empty(shape=[self.nlev,self.nlat,self.nlon])
            #print(np.isfortran(self.var[var]))
            # assuming levels are stored from top to bottom
            offset = TT[0]['level']
            for l in range(len(TT)):
                self.var[var][TT[l]['level']-offset,:,:] = TT[l]['values']
                if readlev:
                    try:
                        lev = TT[l]['lev']
                    except:
                        pass
                    try:
                        lev = TT[l]['level']
                    except:
                        pass
                    self.attr['levs'][l] = lev
                    self.attr['plev'][l] = self.attr['am'][lev-1] + self.attr['bm'][lev-1]*cst.pref
            del self.TT

#       # Check the vertical ordering of the file
        if readlev:
            if not strictly_increasing(self.attr['levs']):
                self.warning.append('NOT STRICTLY INCREASING LEVELS')
        # Reverting North -> South to South to North
        #print(np.isfortran(self.var[var]))
        self.var[var] = self.var[var][:,::-1,:]
        #print(np.isfortran(self.var[var]))
        return None

    def get_var(self,var):
        self._get_var(var)
        return self.var['var' ]

    def _mkp(self,suffix=''):
        # Calculate the pressure field
        self.var['P'+suffix] =  np.empty(shape=(self.nlev,self.nlat,self.nlon))
        for i in range(self.nlev):
            lev = self.attr['levs'][i]
            self.var['P'+suffix][i,:,:] = self.attr['am'][lev-1] \
                                 + self.attr['bm'][lev-1]*self.var['SP'+suffix]

    def _mkpz(self):
        # Calculate pressure field for w (check)
        self.var['PZ'] =  np.empty(shape=(self.nlev,self.nlat,self.nlon))
        for i in range(self.nlev):
            lev = self.attr['levs'][i]
            self.var['PZ'][i,:,:] = self.attr['ai'][lev-1] \
                                  + self.attr['bi'][lev-1]*self.var['SP']

    def _mkpscale(self):
        # Define the standard pressure scale for this vertical grid
        # could also be defined as as d1d field
        self.attr['pscale'] = self.attr['am'] + self.attr['bm'] * 101325
        self.attr['pscale_i'] = self.attr['ai'] + self.attr['bi'] * 101325

    def _mkzscale(self):
        # Define the standard altitude scale
        from zISA import zISA
        try:
            self.attr['zscale'] = zISA(self.attr['pscale'])
            self.attr['zscale_i'] = zISA(self.attr['pscale_i'][1:])
            ztop = 2*self.attr['zscale_i'][0]-self.attr['zscale_i'][1]
            self.attr['zscale_i'] = np.concatenate(((ztop,),self.attr['zscale_i']))
        except:
            print('pscale must be generated before zscale')

    def _mkthet(self,suffix=''):
        # Calculate the potential temperature
        if not set(['T'+suffix,'P'+suffix]).issubset(self.var.keys()):
            print('T or P undefined')
            return
        self.var['PT'+suffix] = self.var['T'+suffix] * (cst.p0/self.var['P'+suffix])**cst.kappa

    def _mkrho(self,suffix=''):
        # Calculate the dry density
        if not set(['T'+suffix,'P'+suffix]).issubset(self.var.keys()):
            print('T or P undefined')
            return
        self.var['RHO'+suffix] = (1/cst.R) * self.var['P'+suffix] / self.var['T'+suffix]

    def _mkrhoq(self):
        # Calculate the moist density
        if not set(['T','P','Q']).issubset(self.var.keys()):
            print('T, P or Q undefined')
            return
        pcor = 230.617*self.var['Q']*np.exp(17.5043*self.var['Q']/(241.2+self.var['Q']))
        self.var['RHO'] = (1/cst.R) * (self.var['P'] - pcor) / self.var['T']

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

    def _mkz(self,suffix=''):
        """ Calculate the geopotential altitude (m) without taking moisture into account """
        if not set(['T','P']).issubset(self.var.keys()):
            print('T or P undefined')
            return
        try:
            with gzip.open(os.path.join(self.rootdir,'EN-true','Z0_'+self.project+'.pkl')) as f:
                Z0 = pickle.load(f)
        except:
            print('Cannot read ground geopotential altitude')
        try:
            self.var['Z0'] = Z0.var['Z0']
        except:
            self.var['Z0'] = Z0
        # processing special cases
        if self.hemis == 'SH':
            self.var['Z0'] = self.var['Z0'][0:self.nlat,:]
        if self.hemis == 'NH':
            # we assume here the ground geopotential needs to be rotated and truncated
            # works here for the OPZ case (Californian fire)
            print('truncate and rotate Z0')
            ll0 = 181
            self.var['Z0'] = np.concatenate((self.var['Z0'][self.nlat-1:,ll0:],self.var['Z0'][self.nlat-1:,:ll0]),axis=1)
        self.var['Z'+suffix] =  np.empty(shape=self.var['T'+suffix].shape)
        uu = - np.log(self.var['P'+suffix])
        uusp = -np.log(self.var['SP'+suffix])
        self.var['Z'+suffix][self.nlev-1,:,:] = self.var['Z0'] \
               + (cst.R/cst.g) * self.var['T'+suffix][self.nlev-1,:,:] \
               * (uu[self.nlev-1,:,:]-uusp)
        for i in range(self.nlev-2,-1,-1):
             self.var['Z'+suffix][i,:,:] = self.var['Z'+suffix][i+1,:,:] + 0.5*(cst.R/cst.g) \
                    * (self.var['T'+suffix][i,:,:] + self.var['T'+suffix][i+1,:,:])\
                    * (uu[i,:,:]-uu[i+1,:,:])

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
        if self.globalGrid:
            dPTdLam[:,1:-1,1:-1] = 0.5*(self.var['PT'+suffix][:,1:-1,2:]-self.var['PT'+suffix][:,1:-1,:-2])/dlon
            dPTdLam[:,1:-1,0] = 0.5*(self.var['PT'+suffix][:,1:-1,1]-self.var['PT'+suffix][:,1:-1,-1])/dlon
            dPTdLam[:,1:-1,-1] = 0.5*(self.var['PT'+suffix][:,1:-1,0]-self.var['PT'+suffix][:,1:-1,-2])/dlon
            dPTdLam[:,1:-1,:] /= coslat[np.newaxis,1:-1,np.newaxis]
            dPTdLam[:,0,:] = dPTdLam[:,1,:]
            dPTdLam[:,-1,:] = dPTdLam[:,-2,:]
        else: # it is assumed that the non global grid does not contain poles
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

    def _mkinc(self,varss='All' ):
        """ Calculate the increment on the listed variables or on all the variables in the forecast
        Beware that the analysis and forecast should be on same grid. If not make a step to meet this
        criterion before """
        if (self.F12_open == False) & (self.CF12_open == False):
            print ('The increment can only be calculated if a forecast file is open')
            return
        if varss == 'All': listvar = [v for v in ['U','V','O3','VO','T'] if v in self.var]
        else: listvar = varss
        for var in listvar:
            if any(x not in self.var for x in [var,var+'F']):
                print('one of the two fields',var,var+'F','not available')
                continue
            if self.var[var].shape != self.var[var+'F'].shape :
                print('shape of analysis and forecast do not correspond for var ',var)
            else:
                self.var[var+'D'] = self.var[var]-self.var[var+'F']
        return

if __name__ == '__main__':
    date = datetime(2017,8,11,18)
    dat = ECMWF('STC',date)
    dat._get_T()
    dat._get_U()
    dat._get_var('ASLWR')
    #dat._get_var('CC')
    dat._mkp()
    dat._mkthet()
    dat._checkThetProfile()
