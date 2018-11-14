# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 01:15:15 2016
Exploit the trajectories to investigate the transit properties from sources to 
target levels
Modified 30 March 2018 to process the FORW and FORWN runs.

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
from mpl_toolkits.basemap import Basemap
from scipy.ndimage import gaussian_filter
#import gzip, pickle
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
            self.transit[var+'_s_go'] = np.zeros(shape=[self.binv,self.source['biny'],self.source['binx']])
            self.transit[var+'_t_go'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])    
    
        self.water_path = water_path
        if water_path:
            self.transit['rv_t'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
            self.transit['rv_t_go'] = np.zeros(shape=[self.binv,self.target['biny'],self.target['binx']])
        self.transit['count']=0
        self.transit['count_go']=0
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
        for j in range(self.binv):
            #print(j,'.',end="")
            sys.stdout.write(str(j)+'.')
            sys.stdout.flush()
            # Select the parcels located within a target layer
            selec=(idv==j)
            lsel=np.sum(selec)
            if lsel==0:
                continue
            # Get the good-opaq among them
            go = dat['go'][selec]
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
            ww = area_pix * np.cos(np.deg2rad(idy0))
            # Histograms in the source and target space
            H_t,_,_=np.histogram2d(yy_t,xx_t,weights=ww,bins=[self.target['biny'],self.target['binx']],
                         range=np.flip(self.target['range'],0))
            H_s,_,_=np.histogram2d(y0_s,x0_s,weights=ww,bins=[self.source['biny'],self.source['binx']],
                         range=np.flip(self.source['range'],0))            
            H_t_go,_,_=np.histogram2d(yy_t[go],xx_t[go],weights=ww[go],bins=[self.target['biny'],self.target['binx']],
                         range=np.flip(self.target['range'],0))
            H_s_go,_,_=np.histogram2d(y0_s[go],x0_s[go],weights=ww[go],bins=[self.source['biny'],self.source['binx']],
                         range=np.flip(self.source['range'],0))
            
            self.transit['count'] += len(yy_t)
            self.transit['count_go'] += np.sum(go)
            
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
            np.add.at(self.transit['totage_t_go'][j,:,:],(idyt[go],idxt[go]),age[go]*ww[go])
            np.add.at(self.transit['totage_s_go'][j,:,:],(idy0[go],idx0[go]),age[go]*ww[go])
            np.add.at(self.transit['totdz_t_go'][j,:,:],(idyt[go],idxt[go]),dz[go]*ww[go])
            np.add.at(self.transit['totdz_s_go'][j,:,:],(idy0[go],idx0[go]),dz[go]*ww[go])
            np.add.at(self.transit['totz_t_go'][j,:,:],(idyt[go],idxt[go]),alt_t[go]*ww[go])
            np.add.at(self.transit['totz_s_go'][j,:,:],(idy0[go],idx0[go]),alt_s[go]*ww[go])
            np.add.at(self.transit['totdx_t_go'][j,:,:],(idyt[go],idxt[go]),dx[go]*ww[go])
            np.add.at(self.transit['totdx_s_go'][j,:,:],(idy0[go],idx0[go]),dx[go]*ww[go])
            np.add.at(self.transit['totdy_t_go'][j,:,:],(idyt[go],idxt[go]),dy[go]*ww[go])
            np.add.at(self.transit['totdy_s_go'][j,:,:],(idy0[go],idx0[go]),dy[go]*ww[go])
            np.add.at(self.transit['totdx2_t_go'][j,:,:],(idyt[go],idxt[go]),dx2[go]*ww[go])
            np.add.at(self.transit['totdx2_s_go'][j,:,:],(idy0[go],idx0[go]),dx2[go]*ww[go])
            np.add.at(self.transit['totdy2_t_go'][j,:,:],(idyt[go],idxt[go]),dy2[go]*ww[go])
            np.add.at(self.transit['totdy2_s_go'][j,:,:],(idy0[go],idx0[go]),dy2[go]*ww[go])
            np.add.at(self.transit['totthet_t_go'][j,:,:],(idyt[go],idxt[go]),thet_t[go]*ww[go])
            np.add.at(self.transit['totthet_s_go'][j,:,:],(idy0[go],idx0[go]),thet_s[go]*ww[go])
            np.add.at(self.transit['totdthet_t_go'][j,:,:],(idyt[go],idxt[go]),dthet[go]*ww[go])
            np.add.at(self.transit['totdthet_s_go'][j,:,:],(idy0[go],idx0[go]),dthet[go]*ww[go])
            if self.water_path:
                np.add.at(self.transit['rv_t'][j,:,:],(idyt,idxt),rv_t*ww)
                np.add.at(self.transit['rv_t_go'][j,:,:],(idyt[go],idxt[go]),rv_t[go]*ww[go])                 
            self.transit['hist_s'][j,:,:]+=H_s
            self.transit['hist_t'][j,:,:]+=H_t
            self.transit['hist_s_go'][j,:,:]+=H_s_go
            self.transit['hist_t_go'][j,:,:]+=H_t_go
        #print("")        
        return

    def complete(self):
        H_s=self.transit['hist_s']+0. # to avoid identification and unwished mod*
        H_t=self.transit['hist_t']+0.
        H_s_go=self.transit['hist_s_go']+0.
        H_t_go=self.transit['hist_t_go']+0.
        self.transit['Hsum_s']=np.sum(H_s,axis=(1,2))
        self.transit['Hsum_t']=np.sum(H_t,axis=(1,2))
        self.transit['Hsum_s_go']=np.sum(H_s_go,axis=(1,2))
        self.transit['Hsum_t_go']=np.sum(H_t_go,axis=(1,2))
        Hsum_s=self.transit['Hsum_s']+0. # to avoid identification and unwished mod
        Hsum_t=self.transit['Hsum_t']+0.
        Hsum_s_go=self.transit['Hsum_s_go']+0.
        Hsum_t_go=self.transit['Hsum_t_go']+0.
        Hsum_s[Hsum_s==0]=1
        Hsum_t[Hsum_t==0]=1
        Hsum_s_go[Hsum_s_go==0]=1
        Hsum_t_go[Hsum_t_go==0]=1
        self.transit['Hnorm_s']=H_s/Hsum_s[:,None,None]
        self.transit['Hnorm_t']=H_t/Hsum_t[:,None,None]
        self.transit['Hnorm_s_go']=H_s_go/Hsum_s_go[:,None,None]
        self.transit['Hnorm_t_go']=H_t_go/Hsum_t_go[:,None,None]
        H_s[H_s==0]=1
        H_t[H_t==0]=1
        H_s_go[H_s_go==0]=1
        H_t_go[H_t_go==0]=1
        for var in ['age','dz','z','dx','dy','dx2','dy2','thet','dthet']:
            self.transit['m'+var+'_s'] = self.transit['tot'+var+'_s']/H_s
            self.transit['m'+var+'_t'] = self.transit['tot'+var+'_t']/H_t
            self.transit['m'+var+'_s_go'] = self.transit['tot'+var+'_s_go']/H_s_go
            self.transit['m'+var+'_t_go'] = self.transit['tot'+var+'_t_go']/H_t_go
        if self.water_path:
            self.transit['mrv_t']=self.transit['rv_t']/H_t
            self.transit['mrv_t_go']=self.transit['rv_t_go']/H_t_go
        return
                
    def chart(self,field,lev,vmin=0,vmax=0,back_field=None,txt="",fgp="",cumsum=False,show=True):
        """ Plots a 2d array field with colormap between min and max.
        This array is extracted from a 3D field with level lev.
        Optionnaly, a background field is added as contours. This field is smoothed by default.
        The limit values for the main field are given in vmin and vmax. 
        TODO: add a parameter to fix the levels of the background field (levels argument of contour).
        """
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
        if '_s' in field:
            irange = self.source['range']
            xedge = self.source['xedge']
            yedge = self.source['yedge']
            xcent = self.source['xcent']
            ycent = self.source['ycent']
        elif '_t' in field:
            irange = self.target['range']
            xedge = self.target['xedge']
            yedge = self.target['yedge']
            xcent = self.target['xcent']
            ycent = self.target['ycent']
        m = Basemap(projection='cyl',llcrnrlat=irange[1,0],urcrnrlat=irange[1,1],
            llcrnrlon=irange[0,0],urcrnrlon=irange[0,1],resolution='c')
        m.drawcoastlines(color='w'); m.drawcountries(color='k')
        if '_s' in field:
            meridians = np.arange(-10.,180.,20.)
            parallels = np.arange(0.,50.,10.)
        else:
            meridians = np.arange(-180.,180.,30.)
            parallels = np.arange(-80.,80.,20.)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
        if vmin==0:    
            vmin=np.min(self.transit[field][lev,:,:])
        if vmax==0:    
            vmax=np.max(self.transit[field][lev,:,:])
        bounds=np.arange(vmin,vmax*(1+0.0001),(vmax-vmin)/mymap.N)
        norm=colors.BoundaryNorm(bounds,mymap.N)
        #print('chart')
        #print(field,self.transit[field].shape,vmin,vmax,range)
        #iax=ax.imshow(self.transit[field][lev,:,:],interpolation='nearest',extent=irange.flatten(),
        #               clim=[vmin,vmax],origin='lower',cmap=mymap,norm=norm,aspect=1.)
        # Plot of the background field if needed. This field is smoothed before plotting by default in order
        # to improv the appearance
        if back_field is not None:            
            if cumsum:
                h,edges = np.histogram(self.transit[back_field][lev,:,:],bins=50)
                cc = 0.5*(edges[:-1]+edges[1:])
                ee = np.cumsum((cc*h/np.sum(cc*h))[::-1])[::-1]         
                nl = [np.argmin(np.abs(ee-x)) for x in [0.9,0.7,0.5,0.3,0.1]]
                CS=ax.contour(xcent,ycent,gaussian_filter(self.transit[back_field][lev,:,:],2),levels=edges[nl])
                strs = ['90%','70%','50%','30%','10%']
                fmt={}
                for l,s in zip(CS.levels,strs):
                    fmt[l] = s
                #plt.clabel(CS,CS.levels[::2],inline=True,fmt=fmt,fontsize=12)
                plt.clabel(CS,inline=True,fmt=fmt,fontsize=12)
            else:
                CS=ax.contour(xcent,ycent,gaussian_filter(self.transit[back_field][lev,:,:],2))                
                plt.clabel(CS)
        # Plot of the main field
        iax=ax.pcolormesh(xedge,yedge,self.transit[field][lev,:,:],
                       clim=[vmin,vmax],cmap=mymap,norm=norm)
        ax.tick_params(labelsize=16)
        plt.title(txt,fontsize=18)
        #plt.xlabel('longitude')
        #plt.ylabel('latitude')
        cax = fig.add_axes([0.91, 0.26, 0.03, 0.5])
        cbar=fig.colorbar(iax,cax=cax)
        cbar.ax.tick_params(labelsize=18)
        if len(fgp)>0:
            plt.savefig('figs/chart-'+fgp+'.png')
        if show:
            plt.show()
    
    def vect(self,lev,thresh=0.00025,type='source',txt="",fgp="",show=True):
        """ Plots a 2d array field with colormap between min and max"""
        try:
            n1=plt.get_fignums()[-1]+1
            plt.figure(plt.get_fignums()[-1]+1,figsize=[13,6])
        except:
            n1=1
            plt.figure(n1+1,figsize=[13,6])
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
        m = Basemap(projection='cyl',llcrnrlat=self.target['range'][1,0],urcrnrlat=self.target['range'][1,1],
            llcrnrlon=self.target['range'][0,0],urcrnrlon=self.target['range'][0,1],resolution='c')
        m.drawcoastlines(color='k'); m.drawcountries(color='k')
        if self.target['range'][1,0] < -50:
            meridians = np.arange(-170.,185.,20.); parallels = np.arange(-80.,85.,20.)
        else:
            meridians = np.arange(-10.,165.,5.); parallels = np.arange(0.,55.,5.)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
        
        #bounds=np.arange(vmin,vmax*(1+0.0001),(vmax-vmin)/mymap.N)
        #norm=colors.BoundaryNorm(bounds,mymap.N)
        #iax=plt.imshow(field.T,interpolation='nearest',extent=source_range.flatten(),
        #               clim=[vmin,vmax],origin='lower',cmap=mymap,norm=norm,aspect=1.)
        # Sparsing the data  
        plt.quiver(xcent_e,ycent_e,dx_e,dy_e,H_e,angles='xy',scale_units='xy',scale=1)    
        plt.title(txt)
        #plt.xlabel('longitude')
        #plt.ylabel('latitude')
        #cax = fig.add_axes([0.91, 0.26, 0.03, 0.5])
        #fig.colorbar(iax,cax=cax)
        if len(fgp)>0:
            plt.savefig('figs/vect-'+fgp+'.png')
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
            plt.savefig('figs/chartv-'+fgp+'.png')
        #fig.suptitle("txt"+" lat="+str(lat))
        if show:
            plt.show()
