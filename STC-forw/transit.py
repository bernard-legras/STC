# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 01:15:15 2016
Exploit the trajectories to investigate the transit properties from sources to 
target levels

The methods of transit class are:
    
__init__: initializes the source and target dictionaries
For all quantities defined in these dictionary, the horizontal grid is either in the source
or the target space. The vertical levels always refer to the target space. Therefore this tool
cannot be used to study the vertical distribution of sources, exccept by the mean theta or z and
the mean displacements.
Variables diagnosed in transit
'hist': count of the parcels,
'age': mean age,
'z': mean altitude,
'dz': mean altitude displacement,
'dx': mean longitude displacement,
'dy': mean latitude displacement,
'dx2': mean squarred longitude displacement,
'dy2': mean squarred latitude displacement,
'thet': mean potential temperature,
'dthet': mean potential temperature displacement
'r_v': water mixing ratio (if water_path argument is true)
The prefix tot in __init__ means that these quantities are first accumulated and 
then divided by the count registry (in hist)
The accumulation is made in "update" and the division is made in "complete".
The grid domain is either fullAMA (for source an target space) or global (for target space only).
The grids are 1° resolution in [10W, 160E, 0, 50N] and [179W, 181E, 90S, 90N]
made of 170x50 and 360x180 meshes. The parcels are localized within the mesh.
A centered grid is also defined with the same dimensions.
The vertical levels are defined either in termes of barometric altitude or in terms of potential temperature.
The current definition is such that in pressure, there are 21 levels of 500m each from 10 to 20 km, and
in potential temperature 20 levels of 5K each from 325 to 420 K.

update: Performs the accumulation.
The accumulation consist in counting all the parcel provided with dat in both the source and the 
target grid (_s and _t suffixes).
A weight is applied due to the pixelisation of the sources in the PTOP files that serve to initialize
the forward runs (see prepforw5Box and STC-SAFNWC/). The weight is 0.1°x0.1°xcos(lat source). See also notes 
on impact scaling.

complete: Calculate the means from the accumulations and the counts.

merge: Merge two accumulations

chart: makes a density plot of a field at a given level, in addition can plot another field
as contours with the option of choosing and labelling contour from its cumulative sum starting
from the maximum value

vect: shows the horizontal displacement as a vector field.

chartv: shows a vertical longitude x altitude section

Modified 30 March 2018 to process the FORW and FORWN runs.
Modified 16 January 2019 to accommodate the new calculations based on meanhigh selection
and doing a further selection on veryhigh (_vh index) or silviahigh (_sh index) in the analysis.
Vect changed to cartopy but untested

This code needs a correction near line 208 for future usages if the resolution of the 
analysis grid is no longer 1°

Fixing data error:
A patch is applied to filter out the contribution of 30 August at 11:00 in the statistics.
No filter for spurious high lat contributions at the moment.

@author: Bernard Legras
"""

#import os
#from struct import unpack
import numpy as np
#import io107
#import matplotlib
#matplotlib.use('Agg') # to avoid requiring X output
#matplotlib.use('Qt4Agg') # anaconda spyder2 on Windows 10 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from cartopy import feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from scipy.ndimage import gaussian_filter
#import gzip, pickle
from os.path import join
import sys
from constants import kappa
from zISA import z1,z2,z3

# Discretization of age histograms
#minage = 0.
#maxage = 20. 
#binage = 50

# No more pre-calculated area field. Weight is included in the calculation below
# area_pix: area of the pixel in the projected SAFNWC map (in degree^2)
area_pix = 0.1*0.1

# Color list with 20 colors      
listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']            
mymap=colors.ListedColormap(listcolors)

p00=100000.

bdget = lambda ll,ss : np.argmin(np.abs(ll-ss))

class transit(object):
    """ Class defined to contain the transit properties of forward runs """

    def __init__(self,water_path=False,target='FullAMA',vert='baro'):
        # Horizontal source range with 1° bins
        # ACHTUNG ACHTUNG: is this is changed, the weight calculation made near line 200
        # might be wrong
        self.source = {}
        self.source['range'] = np.array([[-10.,160.],[0.,50.]])
        self.source['binx'] = 170
        self.source['biny'] = 50
        self.source['xedge'] = np.arange(self.source['range'][0,0],self.source['range'][0,1]+0.001,
                   (self.source['range'][0,1]-self.source['range'][0,0])/self.source['binx'])
        self.source['yedge'] = np.arange(self.source['range'][1,0],self.source['range'][1,1]+0.001,
                   (self.source['range'][1,1]-self.source['range'][1,0])/self.source['biny'])
        self.source['xcent'] = 0.5*(self.source['xedge'][0:-1]+self.source['xedge'][1:])
        self.source['ycent'] = 0.5*(self.source['yedge'][0:-1]+self.source['yedge'][1:])
        # Horizontal target range with 1° bins
        if target =='global':
            self.target = {}
            self.target['type'] = 'global'
            self.target['range'] = np.array([[-179.,181.],[-90.,90.]])
            self.target['binx'] = 360 
            self.target['biny'] = 180
            self.target['xedge'] = np.arange(self.target['range'][0,0],self.target['range'][0,1]+0.001,
                       (self.target['range'][0,1]-self.target['range'][0,0])/self.target['binx'])
            self.target['yedge'] = np.arange(self.target['range'][1,0],self.target['range'][1,1]+0.001,
                       (self.target['range'][1,1]-self.target['range'][1,0])/self.target['biny'])
            self.target['xcent'] = 0.5*(self.target['xedge'][0:-1]+self.target['xedge'][1:])
            self.target['ycent'] = 0.5*(self.target['yedge'][0:-1]+self.target['yedge'][1:])
        else:
            self.target = self.source
            self.target['type'] = 'FullAMA'
        # Vertical discretization
        if vert == 'baro':
            self.vertype = 'baro'
            # Vertical discretization of the barometric altitude (in km)
#            self.minv = 10.
#            self.maxv = 20. 
#            self.binv = 20
            # New definition centered on full levels
            self.minv = 9.75
            self.maxv = 20.25
            self.binv = 21
        elif  'thet' in vert:
            self.vertype = 'theta'
#            self.minv = 325
#            self.maxv = 425
#            self.binv = 20
            # New definition centered on full levels
            self.minv = 322.5
            self.maxv = 422.5
            self.binv = 20
        else:
            print('unknown vertical discretization \n IMMEDIATE ABORT')
            exit()
        self.deltav = (self.maxv-self.minv)/self.binv
        self.vcent = np.arange(self.minv+0.5*self.deltav,self.maxv,self.deltav)
        self.vedge = np.arange(self.minv,self.maxv+0.01,(self.maxv-self.minv)/self.binv) 
        """ Initialization of the transit dictionary """
        self.transit={}
        for var in ['hist','totage','totz','totdz','totdx','totdy','totdx2','totdy2','totdthet','totthet']:
            self.transit[var+'_s'] = np.zeros(shape=[self.binv,self.source['biny'],self.source['binx']])
            self.transit[var+'_t'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
            self.transit[var+'_s_vh'] = np.zeros(shape=[self.binv,self.source['biny'],self.source['binx']])
            self.transit[var+'_t_vh'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
            self.transit[var+'_s_sh'] = np.zeros(shape=[self.binv,self.source['biny'],self.source['binx']])
            self.transit[var+'_t_sh'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
    
        self.water_path = water_path
        if water_path:
            self.transit['rv_t'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
            self.transit['rv_t_vh'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
            self.transit['rv_t_sh'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
            self.transit[var+'_s_vh'] = np.zeros(shape=[self.binv,self.source['biny'],self.source['binx']])
            self.transit[var+'_t_sh'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
        self.transit['count']=0
        self.transit['count_vh']=0
        self.transit['count_sh']=0
        return
    
    def update(self,dat):
        """ Main calculation routine to be applied to each data dictionary read 
        from part files to update the statistics of transit."""
        # initialization of temporary fields
        print('enter transit_update')
        #H_s=H_t=np.zeros(shape=[binv,source_binx,source_biny])
        
        # Calculation of the barometric altitude and the potential temperature
        altbaro=z2(0.01*dat['p'])
        altbaro = np.empty(len(dat['p']))
        id1 = dat['p']>22632.
        altbaro[~id1]=z2(0.01*dat['p'][~id1])
        altbaro[id1]=z1(0.01*dat['p'][id1])
        id1=dat['p']<5474.88
        altbaro[id1]=z3(0.01*dat['p'][id1])
        thet = dat['t']*(p00/dat['p'])**kappa
        
        # Indexing of the vertical location with vedge
        print('digitizing')
        if self.vertype == 'baro':
            idv=np.digitize(altbaro,self.vedge)-1
        elif self.vertype == 'theta':
            idv=np.digitize(thet,self.vedge)-1
        # here we select the target level j
        for j in range(self.binv):
            #print(j,'.',end="")
            sys.stdout.write(str(j)+'.')
            sys.stdout.flush()
            # Select the parcels located within a target layer
            selec=(idv==j)
            lsel=np.sum(selec)
            if lsel==0:
                continue
            # Make a veryhigh boolean slice among active parcels from layer j 
            vh = dat['veryhigh'][selec]
            sh = dat['silviahigh'][selec]
            print('transit',j,np.sum(vh),np.sum(sh))
            # Extract the source and target location and the motion
            xx_t=dat['x'][selec];  yy_t=dat['y'][selec]
            x0_s=dat['x0'][selec]; y0_s=dat['y0'][selec]
            # Extract water vapour mixing ratio
            if self.water_path:
                rv_t = dat['rv_t'][selec]
            # age
            age = dat['age'][selec]
            alt_t = altbaro[selec]
            alt_s = dat['alt0'][selec]
            thet_t = thet[selec]
            thet_s = dat['thet0'][selec]
            dx = xx_t-x0_s 
            dy = yy_t-y0_s
            dz = alt_t - alt_s
            dx2 = dx**2 
            dy2 = dy**2
            dthet = thet_t-thet_s
            # index the latitude and longitude within the transit grid
            # for source and target locations 
            # TO DO target grid for idxt and idyt
            idxt=np.digitize(xx_t,self.target['xedge'])-1
            idyt=np.digitize(yy_t,self.target['yedge'])-1
            idx0=np.digitize(x0_s,self.source['xedge'])-1
            idy0=np.digitize(y0_s,self.source['yedge'])-1
            # todo : diagnose here number of clipped parcels
            # clipping for non managed roundoff effects
            # This procedure may hide errors.
            idxt=np.clip(idxt,0,self.target['binx']-1)
            idyt=np.clip(idyt,0,self.target['biny']-1)
            idx0=np.clip(idx0,0,self.source['binx']-1)
            idy0=np.clip(idy0,0,self.source['biny']-1)
            # Determination of the weight according to the location
            # ACHTUNG ACHTUNG
            # This is OK by chance because the source domain starts from 0 and
            # the resolution of the grid is 1° but it will turn wrong if the origin
            # or the resolution is changed
            ww = area_pix * np.cos(np.deg2rad(idy0))
            # should be instead
            # ww = area_pix * np.cos(np.deg2rad(y0_s))
            
            # Histograms in the source and target space
            # notice, however, that the cumul is done by levels in the target space
            # the sum of H_t and H_s per level should be the same
            H_t,_,_=np.histogram2d(yy_t,xx_t,weights=ww,bins=[self.target['biny'],self.target['binx']],
                         range=np.flip(self.target['range'],0))
            H_s,_,_=np.histogram2d(y0_s,x0_s,weights=ww,bins=[self.source['biny'],self.source['binx']],
                         range=np.flip(self.source['range'],0))            
            H_t_vh,_,_=np.histogram2d(yy_t[vh],xx_t[vh],weights=ww[vh],bins=[self.target['biny'],self.target['binx']],
                         range=np.flip(self.target['range'],0))
            H_s_vh,_,_=np.histogram2d(y0_s[vh],x0_s[vh],weights=ww[vh],bins=[self.source['biny'],self.source['binx']],
                         range=np.flip(self.source['range'],0))
            H_t_sh,_,_=np.histogram2d(yy_t[sh],xx_t[sh],weights=ww[sh],bins=[self.target['biny'],self.target['binx']],
                         range=np.flip(self.target['range'],0))
            H_s_sh,_,_=np.histogram2d(y0_s[sh],x0_s[sh],weights=ww[sh],bins=[self.source['biny'],self.source['binx']],
                         range=np.flip(self.source['range'],0))
            
            self.transit['count'] += len(yy_t)
            self.transit['count_vh'] += np.sum(vh)
            # should add
            # self.transit['count_sh'] += np.sum(vh)
            
            np.add.at(self.transit['totage_t'][j,:,:],(idyt,idxt),age*ww)
            np.add.at(self.transit['totage_s'][j,:,:],(idy0,idx0),age*ww)
            np.add.at(self.transit['totdz_t'][j,:,:],(idyt,idxt),dz*ww)
            np.add.at(self.transit['totdz_s'][j,:,:],(idy0,idx0),dz*ww)
            np.add.at(self.transit['totz_t'][j,:,:],(idyt,idxt),alt_t*ww)
            np.add.at(self.transit['totz_s'][j,:,:],(idy0,idx0),alt_s*ww)
            np.add.at(self.transit['totdx_t'][j,:,:],(idyt,idxt),dx*ww)
            np.add.at(self.transit['totdx_s'][j,:,:],(idy0,idx0),dx*ww)
            np.add.at(self.transit['totdy_t'][j,:,:],(idyt,idxt),dy*ww)
            np.add.at(self.transit['totdy_s'][j,:,:],(idy0,idx0),dy*ww)
            np.add.at(self.transit['totdx2_t'][j,:,:],(idyt,idxt),dx2*ww)
            np.add.at(self.transit['totdx2_s'][j,:,:],(idy0,idx0),dx2*ww)
            np.add.at(self.transit['totdy2_t'][j,:,:],(idyt,idxt),dy2*ww)
            np.add.at(self.transit['totdy2_s'][j,:,:],(idy0,idx0),dy2*ww)
            np.add.at(self.transit['totthet_t'][j,:,:],(idyt,idxt),thet_t*ww)
            np.add.at(self.transit['totthet_s'][j,:,:],(idy0,idx0),thet_s*ww)
            np.add.at(self.transit['totdthet_t'][j,:,:],(idyt,idxt),dthet*ww)
            np.add.at(self.transit['totdthet_s'][j,:,:],(idy0,idx0),dthet*ww)
            
            np.add.at(self.transit['totage_t_vh'][j,:,:],(idyt[vh],idxt[vh]),age[vh]*ww[vh])
            np.add.at(self.transit['totage_s_vh'][j,:,:],(idy0[vh],idx0[vh]),age[vh]*ww[vh])
            np.add.at(self.transit['totdz_t_vh'][j,:,:],(idyt[vh],idxt[vh]),dz[vh]*ww[vh])
            np.add.at(self.transit['totdz_s_vh'][j,:,:],(idy0[vh],idx0[vh]),dz[vh]*ww[vh])
            np.add.at(self.transit['totz_t_vh'][j,:,:],(idyt[vh],idxt[vh]),alt_t[vh]*ww[vh])
            np.add.at(self.transit['totz_s_vh'][j,:,:],(idy0[vh],idx0[vh]),alt_s[vh]*ww[vh])
            np.add.at(self.transit['totdx_t_vh'][j,:,:],(idyt[vh],idxt[vh]),dx[vh]*ww[vh])
            np.add.at(self.transit['totdx_s_vh'][j,:,:],(idy0[vh],idx0[vh]),dx[vh]*ww[vh])
            np.add.at(self.transit['totdy_t_vh'][j,:,:],(idyt[vh],idxt[vh]),dy[vh]*ww[vh])
            np.add.at(self.transit['totdy_s_vh'][j,:,:],(idy0[vh],idx0[vh]),dy[vh]*ww[vh])
            np.add.at(self.transit['totdx2_t_vh'][j,:,:],(idyt[vh],idxt[vh]),dx2[vh]*ww[vh])
            np.add.at(self.transit['totdx2_s_vh'][j,:,:],(idy0[vh],idx0[vh]),dx2[vh]*ww[vh])
            np.add.at(self.transit['totdy2_t_vh'][j,:,:],(idyt[vh],idxt[vh]),dy2[vh]*ww[vh])
            np.add.at(self.transit['totdy2_s_vh'][j,:,:],(idy0[vh],idx0[vh]),dy2[vh]*ww[vh])
            np.add.at(self.transit['totthet_t_vh'][j,:,:],(idyt[vh],idxt[vh]),thet_t[vh]*ww[vh])
            np.add.at(self.transit['totthet_s_vh'][j,:,:],(idy0[vh],idx0[vh]),thet_s[vh]*ww[vh])
            np.add.at(self.transit['totdthet_t_vh'][j,:,:],(idyt[vh],idxt[vh]),dthet[vh]*ww[vh])
            np.add.at(self.transit['totdthet_s_vh'][j,:,:],(idy0[vh],idx0[vh]),dthet[vh]*ww[vh])
            
            np.add.at(self.transit['totage_t_sh'][j,:,:],(idyt[sh],idxt[sh]),age[sh]*ww[sh])
            np.add.at(self.transit['totage_s_sh'][j,:,:],(idy0[sh],idx0[sh]),age[sh]*ww[sh])
            np.add.at(self.transit['totdz_t_sh'][j,:,:],(idyt[sh],idxt[sh]),dz[sh]*ww[sh])
            np.add.at(self.transit['totdz_s_sh'][j,:,:],(idy0[sh],idx0[sh]),dz[sh]*ww[sh])
            np.add.at(self.transit['totz_t_sh'][j,:,:],(idyt[sh],idxt[sh]),alt_t[sh]*ww[sh])
            np.add.at(self.transit['totz_s_sh'][j,:,:],(idy0[sh],idx0[sh]),alt_s[sh]*ww[sh])
            np.add.at(self.transit['totdx_t_sh'][j,:,:],(idyt[sh],idxt[sh]),dx[sh]*ww[sh])
            np.add.at(self.transit['totdx_s_sh'][j,:,:],(idy0[sh],idx0[sh]),dx[sh]*ww[sh])
            np.add.at(self.transit['totdy_t_sh'][j,:,:],(idyt[sh],idxt[sh]),dy[sh]*ww[sh])
            np.add.at(self.transit['totdy_s_sh'][j,:,:],(idy0[sh],idx0[sh]),dy[sh]*ww[sh])
            np.add.at(self.transit['totdx2_t_sh'][j,:,:],(idyt[sh],idxt[sh]),dx2[sh]*ww[sh])
            np.add.at(self.transit['totdx2_s_sh'][j,:,:],(idy0[sh],idx0[sh]),dx2[sh]*ww[sh])
            np.add.at(self.transit['totdy2_t_sh'][j,:,:],(idyt[sh],idxt[sh]),dy2[sh]*ww[sh])
            np.add.at(self.transit['totdy2_s_sh'][j,:,:],(idy0[sh],idx0[sh]),dy2[sh]*ww[sh])
            np.add.at(self.transit['totthet_t_sh'][j,:,:],(idyt[sh],idxt[sh]),thet_t[sh]*ww[sh])
            np.add.at(self.transit['totthet_s_sh'][j,:,:],(idy0[sh],idx0[sh]),thet_s[sh]*ww[sh])
            np.add.at(self.transit['totdthet_t_sh'][j,:,:],(idyt[sh],idxt[sh]),dthet[sh]*ww[sh])
            np.add.at(self.transit['totdthet_s_sh'][j,:,:],(idy0[sh],idx0[sh]),dthet[sh]*ww[sh])
            
            if self.water_path:
                np.add.at(self.transit['rv_t'][j,:,:],(idyt,idxt),rv_t*ww)
                np.add.at(self.transit['rv_t_vh'][j,:,:],(idyt[vh],idxt[vh]),rv_t[vh]*ww[vh]) 
                np.add.at(self.transit['rv_t_sh'][j,:,:],(idyt[sh],idxt[sh]),rv_t[sh]*ww[sh]) 
            self.transit['hist_s'][j,:,:]+=H_s
            self.transit['hist_t'][j,:,:]+=H_t
            self.transit['hist_s_vh'][j,:,:]+=H_s_vh
            self.transit['hist_t_vh'][j,:,:]+=H_t_vh
            self.transit['hist_s_sh'][j,:,:]+=H_s_sh
            self.transit['hist_t_sh'][j,:,:]+=H_t_sh
        #print("")       
        return

    def complete(self):
        H_s=self.transit['hist_s']+0. # to avoid identification and unwished mod*
        H_t=self.transit['hist_t']+0.
        self.transit['Hsum_s']=np.sum(H_s,axis=(1,2))
        self.transit['Hsum_t']=np.sum(H_t,axis=(1,2))
        Hsum_s=self.transit['Hsum_s']+0. # to avoid identification and unwished mod
        Hsum_t=self.transit['Hsum_t']+0.
        Hsum_s[Hsum_s==0]=1
        Hsum_t[Hsum_t==0]=1
        self.transit['Hnorm_s']=H_s/Hsum_s[:,None,None]
        self.transit['Hnorm_t']=H_t/Hsum_t[:,None,None]
        H_s[H_s==0]=1
        H_t[H_t==0]=1
        for var in ['age','dz','z','dx','dy','dx2','dy2','thet','dthet']:
            self.transit['m'+var+'_s'] = self.transit['tot'+var+'_s']/H_s
            self.transit['m'+var+'_t'] = self.transit['tot'+var+'_t']/H_t
        if self.water_path: self.transit['mrv_t'] = self.transit['rv_t']/H_t
        try:
            H_s_vh=self.transit['hist_s_vh']+0.
            H_t_vh=self.transit['hist_t_vh']+0.
            self.transit['Hsum_s_vh']=np.sum(H_s_vh,axis=(1,2))
            self.transit['Hsum_t_vh']=np.sum(H_t_vh,axis=(1,2))
            Hsum_s_vh=self.transit['Hsum_s_vh']+0.
            Hsum_t_vh=self.transit['Hsum_t_vh']+0.
            Hsum_s_vh[Hsum_s_vh==0]=1
            Hsum_t_vh[Hsum_t_vh==0]=1
            self.transit['Hnorm_s_vh']=H_s_vh/Hsum_s_vh[:,None,None]
            self.transit['Hnorm_t_vh']=H_t_vh/Hsum_t_vh[:,None,None]
            H_s_vh[H_s_vh==0]=1
            H_t_vh[H_t_vh==0]=1
            for var in ['age','dz','z','dx','dy','dx2','dy2','thet','dthet']:
                self.transit['m'+var+'_s_vh'] = self.transit['tot'+var+'_s_vh']/H_s_vh
                self.transit['m'+var+'_t_vh'] = self.transit['tot'+var+'_t_vh']/H_t_vh
            if self.water_path: self.transit['mrv_t_vh'] = self.transit['rv_t_vh']/H_t_vh
        except: pass
        try:
            H_s_sh=self.transit['hist_s_sh']+0.
            H_t_sh=self.transit['hist_t_sh']+0.                
            self.transit['Hsum_s_sh']=np.sum(H_s_sh,axis=(1,2))
            self.transit['Hsum_t_sh']=np.sum(H_t_sh,axis=(1,2))             
            Hsum_s_sh=self.transit['Hsum_s_sh']+0.
            Hsum_t_sh=self.transit['Hsum_t_sh']+0.       
            Hsum_s_sh[Hsum_s_sh==0]=1
            Hsum_t_sh[Hsum_t_sh==0]=1       
            self.transit['Hnorm_s_sh']=H_s_sh/Hsum_s_sh[:,None,None]
            self.transit['Hnorm_t_sh']=H_t_sh/Hsum_t_sh[:,None,None]               
            H_s_sh[H_s_sh==0]=1
            H_t_sh[H_t_sh==0]=1       
            for var in ['age','dz','z','dx','dy','dx2','dy2','thet','dthet']:    
                self.transit['m'+var+'_s_sh'] = self.transit['tot'+var+'_s_sh']/H_s_sh
                self.transit['m'+var+'_t_sh'] = self.transit['tot'+var+'_t_sh']/H_t_sh
            if self.water_path: self.transit['mrv_t_sh'] = self.transit['rv_t_sh']/H_t_sh
        except: pass
        return
    
    def merge(self,other):
        """ Merge two classes into a single one """
        self.transit['hist_s'] += other.transit['hist_s']
        self.transit['hist_t'] += other.transit['hist_t']
        try :
            self.transit['hist_s_vh'] += other.transit['hist_s_vh']
            self.transit['hist_t_vh'] += other.transit['hist_t_vh']
        except: pass
        try:
            self.transit['hist_t_sh'] += other.transit['hist_t_sh']
            self.transit['hist_s_sh'] += other.transit['hist_s_sh']
        except: pass
        for var in ['age','dz','z','dx','dy','dx2','dy2','thet','dthet']:
            self.transit['tot'+var+'_s'] += other.transit['tot'+var+'_s']
            self.transit['tot'+var+'_t'] += other.transit['tot'+var+'_t']
            try:
                self.transit['tot'+var+'_s_vh'] += other.transit['tot'+var+'_s_vh']
                self.transit['tot'+var+'_t_vh'] += other.transit['tot'+var+'_t_vh']
            except: pass
            try: 
                self.transit['tot'+var+'_s_sh'] += other.transit['tot'+var+'_s_sh']
                self.transit['tot'+var+'_t_sh'] += other.transit['tot'+var+'_t_sh']
            except: pass
        if self.water_path:
            self.transit['rv_t'] += other.transit['rv_t']
            try: self.transit['rv_t_vh'] += other.transit['rv_t_vh']
            except: pass
            try: self.transit['rv_t_sh'] += other.transit['rv_t_sh']
            except: pass 
                
    def chart(self,field,lev,vmin=0,vmax=0,back_field=None,txt="",fgp="",cumsum=False,show=True,dpi=300):
        """ Plots a 2d array field with colormap between min and max.
        This array is extracted from a 3D field with level lev.
        Optionnaly, a background field is added as contours. This field is smoothed by default.
        The limit values for the main field are given in vmin and vmax. 
        TODO: add a parameter to fix the levels of the background field (levels argument of contour).
        """
        # Set TeX to write the labels and annotations
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        if field not in self.transit.keys():
            print("UNKNOWN FIELD ", field)
            return -1
        elif lev > self.binv-1:
            print('lev TOO LARGE ',lev)
            return -1
        try:
            n1=plt.get_fignums()[-1]+1
            fig=plt.figure(plt.get_fignums()[-1]+1,figsize=[13,6])
        except:
            n1=1
            fig=plt.figure(n1+1,figsize=[13,6])
        ax = fig.add_subplot(111)
        if '_t' in field:
            irange = self.target['range']
            xedge = self.target['xedge']
            yedge = self.target['yedge']
            xcent = self.target['xcent']
            ycent = self.target['ycent']
        elif '_s' in field:
            irange = self.source['range']
            xedge = self.source['xedge']
            yedge = self.source['yedge']
            xcent = self.source['xcent']
            ycent = self.source['ycent']
                       
        fs = 18
        # it is unclear how the trick with cm_lon works in imshow but it does
        # the web says that it is tricky to plot data accross dateline with cartopy
        # check https://stackoverflow.com/questions/47335851/issue-w-image-crossing-dateline-in-imshow-cartopy
        cm_lon =0
        # guess that we want to plot accross dateline (to be improved) 
        if irange[0,1]> 200: cm_lon = 180
        proj = ccrs.PlateCarree(central_longitude=cm_lon)
        fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
        ax = plt.axes(projection = proj)
        if vmin==0:    
            vmin=np.min(self.transit[field][lev,:,:])
        if vmax==0:    
            vmax=np.max(self.transit[field][lev,:,:])
        try:
            bounds=np.arange(vmin,vmax*(1+0.0001),(vmax-vmin)/mymap.N)
        except:
            print('ERROR in bounds vmin vmax',vmin,vmax)
        norm=colors.BoundaryNorm(bounds,mymap.N)
        #print('chart')
        #print(field,self.transit[field].shape,vmin,vmax,range)
        #iax=ax.imshow(self.transit[field][lev,:,:],interpolation='nearest',extent=irange.flatten(),
        #               clim=[vmin,vmax],origin='lower',cmap=mymap,norm=norm,aspect=1.)
        # Plot of the background field if needed. This field is smoothed before plotting by default in order
        # to improv the appearance
        if back_field is not None:            
            if cumsum:
                # Here we plot levels for areas that contain a percentage of the sum
                h,edges = np.histogram(self.transit[back_field][lev,:,:],bins=50)
                cc = 0.5*(edges[:-1]+edges[1:])
                ee = np.cumsum((cc*h/np.sum(cc*h))[::-1])[::-1]         
                nl = [np.argmin(np.abs(ee-x)) for x in [0.9,0.7,0.5,0.3,0.1]]
                CS=ax.contour(xcent,ycent,gaussian_filter(self.transit[back_field][lev,:,:],2),
                              transform=proj,levels=edges[nl],linewidths=3)
                strs = ['90%','70%','50%','30%','10%']                
                fmt={}
                for l,s in zip(CS.levels,strs):
                    fmt[l] = s
                #plt.clabel(CS,CS.levels[::2],inline=True,fmt=fmt,fontsize=12)
                plt.clabel(CS,inline=True,fmt=fmt,fontsize=fs)
            else:
                CS=ax.contour(xcent,ycent,gaussian_filter(self.transit[back_field][lev,:,:],2),
                              transform=proj,linewidths=3)                
                plt.clabel(CS)
        # Plot of the main field
        iax=ax.pcolormesh(xedge,yedge,self.transit[field][lev,:,:],transform=proj,
                       clim=[vmin,vmax],cmap=mymap,norm=norm)
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
        xlocs = None
        if cm_lon == 180: xlocs = [0,30,60,90,120,150,180,-150,-120,-90,-60,-30]
        gl = ax.gridlines(draw_labels=True, xlocs=xlocs,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': fs}
        gl.ylabel_style = {'size': fs}
        #gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
        # Eliminate white borders
        plt.xlim(xedge[0],xedge[-1])
        plt.ylim(yedge[0],yedge[-1])
        plt.title(txt,fontsize=fs)
        # plot adjusted colorbar
        axpos = ax.get_position()
        pos_x = axpos.x0 + axpos.x0 + axpos.width + 0.01
        pos_cax = fig.add_axes([pos_x,axpos.y0,0.04,axpos.height])
        cbar=fig.colorbar(iax,cax=pos_cax)
        cbar.ax.tick_params(labelsize=fs)
        if len(fgp)>0:
            plt.savefig(join('figs','chart-'+fgp+'.png'),dpi=dpi,bbox_inches='tight')
        if show:
            plt.show()
    
    def vect(self,lev,thresh=0.00025,type='source',txt="",fgp="",show=True):
        """ Plots a 2d array field with colormap between min and max"""
        try:
            n1=plt.get_fignums()[-1]+1
            fig = plt.figure(plt.get_fignums()[-1]+1,figsize=[13,6])
        except:
            n1=1
            fig = plt.figure(n1+1,figsize=[13,6])
        # Getting and sparsing the data    
        if 'source' in type:
            suf = '_s'
            xcent_e = self.source['xcent'][1:-1:5]
            ycent_e = self.source['ycent'][1:-1:5]
        else:
            suf = '_t'          
            xcent_e = self.target['xcent'][1:-1:5]
            ycent_e = self.target['ycent'][1:-1:5]
        dx_e=np.ma.array(self.transit['mdx'+suf][lev,1:-1:5,1:-1:5])
        dy_e=np.ma.array(self.transit['mdy'+suf][lev,1:-1:5,1:-1:5])
        H_e=self.transit['Hnorm'+suf][lev,1:-1:5,1:-1:5]
        dx_e[H_e<thresh] = np.ma.masked
        dy_e[H_e<thresh] = np.ma.masked
        if 'inv' in type:
            dx_e = -dx_e
            dy_e = -dy_e
          
        fs = 15
        # it is unclear how the trick with cm_lon works in imshow but it does
        # the web says that it is tricky to plot data accross dateline with cartopy
        # check https://stackoverflow.com/questions/47335851/issue-w-image-crossing-dateline-in-imshow-cartopy
        cm_lon =0
        # guess that we want to plot accross dateline 
        if xcent_e[-1] > 180: cm_lon = 180
        proj = ccrs.PlateCarree(central_longitude=cm_lon)
        fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
        ax = plt.axes(projection = proj)
        
        #bounds=np.arange(vmin,vmax*(1+0.0001),(vmax-vmin)/mymap.N)
        #norm=colors.BoundaryNorm(bounds,mymap.N)
        #iax=plt.imshow(field.T,interpolation='nearest',extent=source_range.flatten(),
        #               clim=[vmin,vmax],origin='lower',cmap=mymap,norm=norm,aspect=1.)
        # Sparsing the data  
        ax.quiver(xcent_e,ycent_e,dx_e,dy_e,H_e,angles='xy',scale_units='xy',scale=1)    
        plt.title(txt,fontsize=fs)
        ax.add_feature(feature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none'))
        ax.coastlines('50m')
        xlocs = None
        if cm_lon == 180: xlocs = [0,30,60,90,120,150,180,-150,-120,-90,-60,-30]
        gl = ax.gridlines(draw_labels=True, xlocs=xlocs,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        #gl.xformatter = LONGITUDE_FORMATTER
        #gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': fs}
        gl.ylabel_style = {'size': fs}
        if len(fgp)>0:
            plt.savefig(join('figs','vect-'+fgp+'.png'),bbox_inches='tight')
        if show:
            plt.show()    
   
    def chartv(self,field,lat,txt="",vmin=0,vmax=0,fgp="",show=True):
        """ Plot a longitude x altitude section"""
        if field not in self.transit.keys():
            print("UNKNOWN FIELD ", field)
            return -1
        try:
            pos=np.where(self.target['ycent']>=lat)[0][0]
        except:
            print('lat out of range') 
            return
        try:
            n1=plt.get_fignums()[-1]+1
        except:
            n1=1
        fig=plt.figure(n1+1,figsize=[13,6])
        ax = fig.add_subplot(111)
        if vmin==0:    
            vmin=np.min(self.transit[field])
        if vmax==0:    
            vmax=np.max(self.transit[field])
        bounds=np.arange(vmin,vmax*(1+0.0001),(vmax-vmin)/mymap.N)
        norm=colors.BoundaryNorm(bounds,mymap.N)
#        if self.vertype == 'baro':
#            aspect = 6
#        elif self.vertype == 'theta':
#            aspect = 0.6
        #iax=ax.imshow(self.transit[field][:,pos, :],interpolation='nearest',
        #               extent=np.append(self.target['range'][0],[self.vcent[0],self.vcent[-1]]),
        #               clim=[vmin,vmax],origin='lower',cmap=mymap,norm=norm,aspect=aspect)
        iax=ax.pcolormesh(self.target['xedge'],self.vedge,self.transit[field][:,pos, :],
                        clim=[vmin,vmax],cmap=mymap,norm=norm)
        ax.tick_params(labelsize=16)
        plt.xlabel('longitude',fontsize=18)
        if self.vertype == 'baro':
            plt.ylabel('baro altitude (km)',fontsize=18)
        elif self.vertype == 'theta':
            plt.ylabel('pot temperature (K)',fontsize=18)
        plt.title(txt+" lat="+str(lat),fontsize=18)
        cax = fig.add_axes([0.91, 0.21, 0.03, 0.6])
        cbar = fig.colorbar(iax,cax=cax)
        cbar.ax.tick_params(labelsize=18)
        if len(fgp)>0:
            plt.savefig(join('figs','chartv-'+fgp+'.png'),bbox_inches='tight')
        #fig.suptitle("txt"+" lat="+str(lat))
        if show:
            plt.show()
