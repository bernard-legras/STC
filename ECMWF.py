#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Opens a FLEXPART type grib file from ECMWF and reads temperature field.
The    

Created on Mon Apr 24 02:43:04 2017

@author: Bernard
"""
from __future__ import division, print_function
#from __future__ import unicode_literals
from datetime import datetime
import numpy as np
import pygrib
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
pref = 101325.
MISSING = -999

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

class ECMWF_pure(object):
    def __init__(self,date):
        self.var={}
        self.attr={}
        self.date = date
        self.warning = []
        
    def chart(self,var,lev=0,txt=''):
        # test existence of key field
        if var not in self.var.keys():
            print ('undefined field')
            return
        
        fig=plt.figure(figsize=[11,4])
        m = Basemap(projection='cyl', llcrnrlat=self.attr['lat'][0], 
                    urcrnrlat=self.attr['lat'][-1],
                    llcrnrlon=self.attr['lon'][0],
                    urcrnrlon=self.attr['lon'][-1], resolution='c')
        m.drawcoastlines(color='k')
        m.drawcountries(color='k')
        if self.attr['lon'][-1] - self.attr['lon'][0] <= 50.:
            spacex = 5.
        else:
            spacex = 10.
        if self.attr['lat'][-1] - self.attr['lat'][0] <= 50.:
            spacey = 5.
        else:
            spacey = 10.
        meridians = np.arange(self.attr['lon'][0], self.attr['lon'][-1], spacex)
        parallels = np.arange(self.attr['lat'][0], self.attr['lat'][-1], spacey)
        m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=8)
        m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=8)
        if len(self.var[var].shape) == 3:
            iax = plt.imshow(self.var[var][lev,:,:], interpolation='nearest', 
                     extent=[self.attr['lon'][0], self.attr['lon'][-1],
                             self.attr['lat'][0], self.attr['lat'][-1]],
                     origin='lower', aspect=1.,cmap='jet')
        else:
            iax = plt.imshow(self.var[var], interpolation='nearest',
                     extent=[self.attr['lon'][0], self.attr['lon'][-1],
                             self.attr['lat'][0], self.attr['lat'][-1]],
                     origin='lower', aspect=1.,cmap='jet')
        cax = fig.add_axes([0.91, 0.15, 0.03, 0.7])
        fig.colorbar(iax, cax=cax)
        plt.title(txt)
        plt.show()
        return None
                
class ECMWF(ECMWF_pure):
    # to do: raise exception in case of an error
    
    def __init__(self,path,date):
        ECMWF_pure.__init__(self,date)
        self.fname = 'EN'+date.strftime('%y%m%d%H')
        print(path,self.fname)
        try:
            self.grb = pygrib.open(os.path.join(path,self.fname))
        except:
            print('cannot open '+self.fname)
        
        try:
            sp = self.grb.select(name='Surface pressure')[0]
            logp = False
        except:
            try:
                sp = self.grb.select(name='Logarithm of surface pressure')[0]
                logp = True
            except:
                print('no surface pressure in '+self.fname)
                self.grb.close()
        # Check time matching  (made from validity date and time)
        vd = sp['validityDate']
        vt = sp['validityTime']
        day = vd % 100
        vd //=100
        month = vd % 100
        vd //=100
        minute = vt % 100
        vt //=100
        dateread = datetime(year=vd,month=month,day=day,hour=vt,minute=minute)
        if dateread != self.date:
            print('WARNING: dates do not match')
            print('called date    '+self.date.strftime('%Y-%m-%d %H:%M'))
            print('date from file '+dateread.strftime('%Y-%m-%d %H:%M'))
            # self.grb.close()   
        # Get general info from this message
        self.attr['Date'] = sp['dataDate']
        self.attr['Time'] = sp['dataTime']
        self.attr['valDate'] = sp['validityDate']
        self.attr['valTime'] = sp['validityTime']         
        self.attr['nlon'] = sp['Ni']
        self.attr['nlat'] = sp['Nj']
        self.attr['lon'] = sp['distinctLongitudes']
        self.attr['lat'] = sp['distinctLatitudes']
        self.attr['Lo1'] = sp['longitudeOfFirstGridPoint']  # in millidegrees
        self.attr['Lo2'] = sp['longitudeOfLastGridPoint']  # in millidegrees
        self.attr['La1'] = sp['latitudeOfFirstGridPoint']  # in millidegrees
        self.attr['La2'] = sp['latitudeOfLastGridPoint'] # in millidegrees
        if sp['PVPresent']==1 :
            pv = sp['pv']
            self.attr['ai'] = pv[0:int(pv.size/2)]
            self.attr['bi'] = pv[int(pv.size/2):]
            self.attr['am'] = (self.attr['ai'][1:] + self.attr['ai'][0:-1])/2
            self.attr['bm'] = (self.attr['bi'][1:] + self.attr['bi'][0:-1])/2
        else: 
            # todo: provide a fix, eg for JRA-55 
            print('missing PV not implemented')          
        # Read the surface pressure
        if logp:
            self.var['SP'] = np.exp(sp['values'])
        else:
            self.var['SP'] = sp['values']
        #  Reverting lat order
        self.attr['lat'] = self.attr['lat'][::-1]
        self.var['SP']   = self.var['SP'][::-1,:]
        
    def close(self):
        self.grb.close()
    
    def _get_T(self,pressure=True):
        # Now proceed to read the temperature field
        self._get_var('T',pressure=pressure)
        
    def _get_U(self,pressure=True):
        self._get_var('U',pressure=pressure)
        
    def _get_V(self,pressure=True):
        self._get_var('V',pressure=pressure)
        
    def _get_W(self,pressure=True):
        self._get_var('W',pressure=pressure)
        
    def _get_Q(self,pressure=True):
        self._get_var('Q',pressure=pressure)     
    
    def _get_var(self,var,pressure=True):
        longname = {'T':'Temperature','U':'U component of wind','V':'V component of wind',
                    'w':'Vertical velocity','Q':'Specific humidity'} 
        try:
            TT = self.grb.select(name=longname[var])
        except: 
            print('no '+var+' in '+self.fname )
        # Process each message corresponding to a level
        self.attr['nlev'] = len(TT)
        self.attr['levs'] = np.full(len(TT),MISSING,dtype=int)
        self.attr['plev'] = np.full(len(TT),MISSING,dtype=int)
        self.var[var] = np.empty(shape=[self.attr['nlev'],self.attr['nlat'],self.attr['nlon']])
        #print(np.isfortran(self.var[var]))
        if pressure:
            if var=='W':
                self.var['PZ'] = np.empty(shape=[self.attr['nlev'],self.attr['nlat'],self.attr['nlon']])
            else:    
                self.var['P'] = np.empty(shape=[self.attr['nlev'],self.attr['nlat'],self.attr['nlon']])
        for i in range(len(TT)):
            try:
                lev = TT[i]['lev']
            except:
                pass
            try:
                lev = TT[i]['level']
            except:
                pass
            self.attr['levs'][i] = lev
            self.var[var][i,:,:] = TT[i]['values']
            #print(np.isfortran(self.var[var]))
            # calculate pressure, reordering in lat unecessary since selfvar['SP'] already reordered
            if pressure:
                if var == 'W':
                    self.var['PZ'][i,:,:] = self.attr['ai'][lev-1] + self.attr['bi'][lev-1]*self.var['SP']
                else:    
                    self.var['P'][i,:,:] = self.attr['am'][lev-1] + self.attr['bm'][lev-1]*self.var['SP']
            self.attr['plev'][i] = self.attr['am'][lev-1] + self.attr['bm'][lev-1]*pref
        # Check the vertical ordering of the file
        if not strictly_increasing(self.attr['levs']):
            self.warning.append('NOT STRICTLY INCREASING LEVELS')
        # Reverting North -> South to South to North
        #print(np.isfortran(self.var[var]))
        self.var[var] = self.var[var][:,::-1,:]
        #print(np.isfortran(self.var[var]))
        return None

    def interpol(self,other,date):
        # check the date
        if self.date < other.date:
            if (date > other.date) | (date < self.date):
                print ('error on date')
                return -1
        else:
            if (date < other.date) | (date > self.date):
                print ('error on date')
                return -1
        # calculate coefficients
        dt = (other.date-self.date).total_seconds()
        dt1 = (date-self.date).total_seconds()
        dt2 = (other.date-date).total_seconds()
        cf1 = dt2/dt
        cf2 = dt1/dt
        print ('cf ',cf1,cf2)
        data = ECMWF_pure(date)
        for names in ['lon','lat','nlon','nlat']:
            data.attr[names] = self.attr[names]
        for var in self.var.keys():
            data.var[var] = cf1*self.var[var] + cf2*other.var[var] 
        return data              
