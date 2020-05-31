#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script produces masks for the FullAMA domain that fill the domain 
(but perhaps a few pixels depending of the resolution)

See the images in the figs directory to see the generated masks

Series DaylyCycle generates
MaskCartopy2-ERA5.pkl
MaskCartopy2-ERA-I.pkl
MaskCartopy2-JRA-55.pkl
MaskCartopy2-MERRA_2.pkl

Series ST generates 
MaskCartopy2-STCforw.pkl
MaskCartopy2-STCforwfine.pkl
MaskCartopy2-STCforwhyperfine.pkl

Series MonthlyMeans generates
STCmask4-ea.pkl
STCmask4-ei.pkl
STCmask4-m2.pkl
STCmask4-j55.pkl

Created on Tue Jul 10 15:38:07 2018
Adapted many times

@author: sbucci, legras
"""

import numpy as np
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.io.shapereader import natural_earth, Reader
from cartopy import feature
#from cartopy.mpl.patch import geos_to_path
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from matplotlib.path import Path
from shapely.geometry import Point
from shapely.ops import cascaded_union
#import collections
import matplotlib.colors as colors
import pickle,gzip
#import deepdish as dd
from scipy.io import savemat

def main():

    # Set the countries to select. Please notice that each region will be overimposed 
    # on the previous ones.
    # So the order of them in thr following dictionary is IMPORTANT
    
    # Parameter choice
    # if DailyCycle
    # rea among ERA5, ERA5-STC, ERA-I, JRA-55, MERRA-2
    # If MonthlyMeans
    # rea among ea, ei, j55, cfsr, m2
    #rea = 'MERRA-2
    rea = 'm2'
    MonthlyMeans = False
    DailyCycle = False
    STCforw = False
    STC_hyperfine = False
    STC_fine = True
    
    # This must be set to 0
    rollx = 0
    
    #      countries          limits                type         color           label-pos
    regions = {
        1:['CentralAfrica',  [[-11,52],[-1,15]],   'land',      'saddlebrown',  [15,5]],
        2:['NorthAfrica',    [[-11,40],[15,40]],   'land',      'sandybrown',   [5,22]],
        3:['Europe',         [[-9,50],[10,52]],    'country',   'lawngreen',    [10,46]],
        4:['Pen',            [[91,120],[-1,30]],   'land',      'darkviolet',   [100,16]],
        5:['NorthAsia',      None,                 'country',   'mediumpurple', [60,47]],
        6:['IndianSub',      None,                 'country',   'blue',         [72,20]],
        7:['NorthChina',     [[60,150],[35,55]],   'country',   'orange',       [102,38]],
        8:['SouthChina',     [[60,150],[13,35]],   'country',   'red',          [104,30]],
        9:['IndoMalaysia',   [[90,150],[-0.5,20]], 'country',   'firebrick',    [107,2]],
        10:['Bangladesh',    None,                 'country',   'limegreen',    [92,22]],
        11:['Pakistan',      None,                 'country',   'c',            [60,25]],
        12:['JapanKorea',    None,                 'country',   'hotpink',      [132,37]],
        13:['IndianOcean',   [[43,100],[-1,30]],   'ocean',     'seagreen',     [57,10]],
        14:['BoB',           [[78,100],[5,25]],    'ocean',     'yellow',       [85,10]],
        15:['MidPacific',    [[110,161],[22,35]],  'ocean',     'cyan',         [135,28]],
        16:['SCSPhi',        [[99,140],[-1,22]],   'ocean',     'coral',        [125,17]],
        17:['MiddleEast',    None,                 'country',   'bisque',       [37,25]],
        18:['WestAsia',      None,                 'country',   'navy',         [50,32]],
        19:['TibetanPlateau',None,                 'orography', 'lightskyblue', [80,32]],
        20:['GuineaGulf',    [[-11,12],[-1,8]],    'ocean',     'aqua',         [-8,2]],
        21:['Caspian',       [[42,60],[35,48]],    'ocean',     'pink',         [50,40]],
        22:['Mediterranea',  [[-1,42],[30,47]],    'ocean',     'turquoise',    [15,35]],
        23:['Atlantic',      [[-11,-1],[25,51]],   'ocean',     'dodgerblue',   [-10,47]],
        24:['RedSea',        [[30,43],[12,30]],    'ocean',     'yellow',       [37,17]],
        25:['WestPacific',   [[140,161],[-1,22]],  'ocean',     'royalblue',    [140,10]],
        26:['NorthPacific',  [[110,161],[35,51]],  'ocean',     'slateblue',    [143,45]],
        27:['Philippines',   None,                 'country',   'grey',         [115,10]]}
    
    #col=['b','lightskyblue','r','yellow','limegreen','coral','darkviolet','cyan','hotpink','k'] 
    #please see https://i.stack.imgur.com/fMx2j.png
    
    xl = {v[0]: v[4][0] for k,v in regions.items()}
    yl = {v[0]: v[4][1] for k,v in regions.items()}
    
    """" Choose the coordinate and the resolution of the mask
     Resolution in degree in both latitude and longitude
      ACHTUNG ACHTUNG icxy must be such that limits are values of type p*icxy 
      where p is integer (where dores this come?, possible misleading)
      
      the (lon,lat) grid is defined as edges of boxes which centers are the model grid
      the (lonc,latc) grid is the model grid, that is used for plots and output """
    
    """ Version of the masks suited for the analysis of the daily cycle for ERA5 using 
     the intermediate products derived from the MNTH data (by mkLTdailycyc) in 
     the LTdailycyc files and from the reanalysis data in the Heating-diag/STCoutpy 
     directories for ERA-I, JRA-55 and MERRA-2. """
    
    if DailyCycle:
    
        if rea == 'ERA5':
            icy = 0.5
            icx = 0.5
        elif rea == 'ERA5-STC':
            icy = 0.25
            icx = 0.25
        elif rea == 'ERA-I':
            icy = 1
            icx = 1
        elif rea == 'MERRA-2':
            icy = 0.5
            icx = 0.625
        if rea in ['ERA5','ERA5-STC','MERRA-2','ERA-I']:
            lon=np.arange(-10-0.5*icx,160+icx,icx) #degrees from -10 to 160 with step icxy
            lat=np.arange(-0.5*icy,50+icy,icy)
            lonc=np.arange(-10,160+0.5*icx,icx)
            latc=np.arange(0,50+0.5*icy,icy)
        if rea == 'JRA-55':
            icx = 0.5625
            icy = 0.5625
            lonc=np.arange(-9.5625,160,icx)
            lon=np.arange(lonc[0]-0.5*icx,160+icx,icx)
            latc = np.array([0.281, 0.842, 1.404, 1.966, 2.527, 3.089, 3.651, 4.212, 4.774, 5.335, \
                            5.897, 6.459, 7.020, 7.582, 8.144, 8.705, 9.267, 9.828, 10.390, 10.952, \
                           11.513, 12.075, 12.636, 13.198, 13.760, 14.321, 14.883, 15.445, 16.006, 16.568, \
                           17.129, 17.691, 18.253, 18.814, 19.376, 19.938, 20.499, 21.061, 21.622, 22.184, \
                           22.746, 23.307, 23.869, 24.431, 24.992, 25.554, 26.115, 26.677, 27.239, 27.800, \
                           28.362, 28.924, 29.485, 30.047, 30.608, 31.170, 31.732, 32.293, 32.855, 33.416, \
                           33.978, 34.540, 35.101, 35.663, 36.225, 36.786, 37.348, 37.909, 38.471, 39.033, \
                           39.594, 40.156, 40.718, 41.279, 41.841, 42.402, 42.964, 43.526, 44.087, 44.649, \
                           45.211, 45.772, 46.334, 46.895, 47.457, 48.019, 48.580, 49.142, 49.704])
            lat = 0.5*(latc[:-1]+latc[1:])
            lat = np.insert(lat,0,lat[0]-icy)
            lat = np.append(lat,lat[-1]+icy)
    
    elif STCforw:
        icy = 1
        icx = 1
        lon = np.arange(-10,160.6,icx)
        lonc = np.arange(-9.5,160,icx)
        lat = np.arange(0,50.6,icy)
        latc = np.arange(0.5,50,icy)   
    elif MonthlyMeans:
        if rea == 'ea':
            icy = 0.5
            icx = 0.5
            lon  = np.arange(-179.75,180.3,icx)
            lonc = np.arange(-179.5,180.1,icx)
            lat = np.arange(-50.25,50.3,icy)
            latc = np.arange(-50,50.1,icy)    
        elif rea == 'ei':
            icy = 1
            icx = 1
            lon = np.arange(-179.5,180.6,icx)
            lonc = np.arange(-179,180.1,icx)
            lat = np.arange(-50.5,50.6,icy)
            latc = np.arange(-50,50.1,icy)
        elif rea == 'm2':
            icy = 0.5
            icx = 0.625
            lonc = np.arange(-180.,179.5,icx)
            lon = np.arange(-180-0.5*icx,179.8,icx)
            lat = np.arange(-50.25,50.3,icy)
            latc = np.arange(-50,50.1,icy)            
        elif rea == 'cfsr':           
            icy = 2.5
            icx = 2.5
            rollx = int(10/2.5)
            lonc = np.arange(-10,357.6-10,icx)
            lon = np.arange(-10-0.5*icx,359-10,icx)
            latc = np.arange(-50,51,icy)
            lat = np.arange(-50-0.5*icy,51.5,icy)
        elif rea == 'flxhr':
            print ('Use ei grid for flxhr')
            exit()
            icy = 5.
            icx = 5.
            lonc = np.arange(-177.5,180,icx)
            lon = np.arange(-180,181,icx)
            latc = np.arange(-47.5,48.,icy)
            lat = np.arange(-50,51,icy)
        elif rea == 'j55':            
            icx = 0.5625
            icy = 0.5625
            shiftx = 9.5625
            rollx = int(shiftx/icx)
            lonc=np.arange(-shiftx,359.5-shiftx,icx)
            lon=np.arange(lonc[0]-0.5*icx,lonc[-1]+icx,icx)
            latc = np.array([0.281, 0.842, 1.404, 1.966, 2.527, 3.089, 3.651, 4.212, 4.774, 5.335, \
                            5.897, 6.459, 7.020, 7.582, 8.144, 8.705, 9.267, 9.828, 10.390, 10.952, \
                           11.513, 12.075, 12.636, 13.198, 13.760, 14.321, 14.883, 15.445, 16.006, 16.568, \
                           17.129, 17.691, 18.253, 18.814, 19.376, 19.938, 20.499, 21.061, 21.622, 22.184, \
                           22.746, 23.307, 23.869, 24.431, 24.992, 25.554, 26.115, 26.677, 27.239, 27.800, \
                           28.362, 28.924, 29.485, 30.047, 30.608, 31.170, 31.732, 32.293, 32.855, 33.416, \
                           33.978, 34.540, 35.101, 35.663, 36.225, 36.786, 37.348, 37.909, 38.471, 39.033, \
                           39.594, 40.156, 40.718, 41.279, 41.841, 42.402, 42.964, 43.526, 44.087, 44.649, \
                           45.211, 45.772, 46.334, 46.895, 47.457, 48.019, 48.580, 49.142, 49.704, 50.265])
            latc = np.append(-np.flip(latc),latc)
            lat = 0.5*(latc[:-1]+latc[1:])
            lat = np.insert(lat,0,lat[0]-icy)
            lat = np.append(lat,lat[-1]+icy)
    
    elif STC_fine:
        icx = 0.25
        icy = 0.25
        lon = np.arange(-10,160.05,icx)
        lat = np.arange(0,50.05,icy)
        lonc = np.arange(-10+0.5*icx,160,icx)
        latc = np.arange(0.5*icy,50,icy)  
    
    elif STC_hyperfine:
        icx = 0.1
        icy = 0.1
        lon = np.arange(-10,160.05,icx)
        lat = np.arange(0,50.05,icy)
        lonc = np.arange(-9.95,160,icx)
        latc = np.arange(0.05,50,icy)   
    
    # Initialize mask, coordinates grid and dictionaries to interpret the mask
    # The list of coordinates above is here seen as a list of edges 
    print('mask of shape ',len(lonc),len(latc))
    mask=np.zeros([len(lonc),len(latc)])
    # numbering of countries starting from 1
    name_to_number = {v[0]: k for k,v in regions.items()}
    number_to_name = {k: v[0] for k,v in regions.items()}
    
    # Generates the mesh
    x, y = np.meshgrid(lon, lat)
    
    #Creates the mask accordingly to the dictionaries defined in the beginning
    for reg,v in regions.items():
        print(v[2])
        mask=add_mask(x,y,mask,v[2],v[0],reg,v[1])
        
    #Roll longitude dimension if needed
    if rollx != 0:
        lonc = np.roll(lonc,-rollx)
        lon = np.roll(lon,-rollx)
        lonc[-rollx:] += 360
        lon[-rollx:] =+ 360
        mask = np.roll(mask,-rollx,axis=0)
        
    #for i in ['Tibetan-Plateau',]:
    #    print(labs[i])
    #    mask=add_mask(x,y,mask,labs[i],i,name_to_number,limits[i])
    #save the mask    
    mask=mask.T    
    mask0={}
    mask0['mask']= mask
    mask0['lons']= lonc
    mask0['lats']= latc
    mask0['regcode'] = name_to_number
    mask0['regcode_inv'] = number_to_name
    mask0['countries'] = [v[0] for k,v in regions.items()]
    mask0['limits'] = {v[0]: v[1] for k,v in regions.items()}
    mask0['ccreg'] =  {v[0]: v[3] for k,v in regions.items()}  
    mask0['icy']=icy
    mask0['icx']=icx
    mask0['lonmin']=np.min(lonc)
    mask0['lonmax']=np.max(lonc)
    mask0['nlons']=len(lonc)
    mask0['latmin']=np.min(latc)
    mask0['latmax']=np.max(latc)
    mask0['nlats']=len(latc)
    mask0['xl']=xl
    mask0['yl']=yl
    
    #%%
    #for i in ['Tibetan-Plateau',]:
    #    print(regions[19][2])
    #    mask=add_mask(x,y,mask,regions[19][2],i,name_to_number,mask0['limits'][i])
    
    #%%%
    if DailyCycle:
        pickle.dump(mask0,gzip.open('MaskCartopy2-'+rea+'.pkl','wb'),pickle.HIGHEST_PROTOCOL)
        #dd.io.save('MaskCartopy.h5',mask0)
    elif STCforw:
        pickle.dump(mask0,gzip.open('MaskCartopy2-STCforw.pkl','wb'),pickle.HIGHEST_PROTOCOL)
    elif MonthlyMeans:
        pickle.dump(mask0,gzip.open('STCmask4-'+rea+'.pkl','wb'),pickle.HIGHEST_PROTOCOL)
    elif STC_hyperfine:
        pickle.dump(mask0,gzip.open('MaskCartopy2-STChyperfine.pkl','wb'),pickle.HIGHEST_PROTOCOL)
    elif STC_fine:
        pickle.dump(mask0,gzip.open('MaskCartopy2-STCfine.pkl','wb'),pickle.HIGHEST_PROTOCOL)
    
    if DailyCycle:
        mask1 = mask0.copy()
        for reg,ll in mask1['limits'].items():
            if ll is None:
                mask1['limits'][reg] = []
        savemat('MaskCartopy2-'+rea+'.mat',mask1)
    
    
    #%% Plot the mask
    #Read the mask
    # mask0=pickle.load(gzip.open('MaskCartopy2.pkl','rb'))
    mask=mask0['mask'].T
    lonc=mask0['lons']
    latc=mask0['lats']
    name_to_number=mask0['regcode']
    number_to_name=mask0['regcode_inv']
    limits=mask0['limits']    
    ccreg=mask0['ccreg']
    countries=mask0['countries']
    
    plot_mask(lon,lat,mask,number_to_name,limits,ccreg,countries,xl,yl)
    
    plot_mask(lon,lat,mask,number_to_name,limits,ccreg,[],xl,yl)

def add_country(geoms,name):
    shape_records = Reader(natural_earth(resolution='110m',
                                         category='cultural',
                                         name='admin_0_countries')).records()
    for country in shape_records:
        if country.attributes['NAME_LONG'] in name:       
            try:
                geoms += country.geometry
            except TypeError:
                geoms.append(country.geometry)
    return geoms

def get_geometries_countries(country_names):
    print(country_names)
    """
    Get an iterable of Shapely geometries corrresponding to given countries.

    """

    if country_names=='NorthChina':
        country_names2 = ['China']  
    
    elif country_names=='SouthChina':
        country_names2 = ['China','Taiwan']
    
    elif country_names=='JapanKorea':
        country_names2 = ['Japan','Republic of Korea','Dem. Rep. Korea']

    elif country_names=='NorthAsia':
        country_names2 = ['Russian Federation','Mongolia','Kazakhstan']
        
    elif country_names=='IndianSub':
        country_names2 = ['India','Nepal','Bhutan','Sri Lanka']

    elif country_names=='MiddleEast':
        country_names2 = ['Jordan','Israel','Saudi Arabia','Yemen','Oman','Lebanon','Iraq','Syria','United Arab Emirates','Kuwait']

    elif country_names=='WestAsia':
        country_names2 = ['Iran','Afghanistan','Turkmenistan','Uzbekistan','Tajikistan','Kyrgyzstan','Azerbaijan','Georgia','Armenia']

    elif country_names=='IndoMalaysia':
        country_names2 = ['Indonesia','Malaysia','Singapour'] 
    
    elif country_names=='Europe':
        country_names2= ['Italy','France','Greece','Macedonia','Albania','Montenegro','Kosovo','Bulgaria','Serbia','Romania',\
                       'Bosnia and Herzegovina','Slovenia','Hungary','Croatia','Slovakia','Austria','Czech Republic','Cyprus',\
                       'Poland','Ukraine','Belarus','Moldova','Lithuania','Latvia','Estonia','Germany','Spain','Belgium','Switzerland',\
                       'Northern Cyprus','Denmark','Finland','Georgia','Luxembourg','Norway','Sweden','Turkey','Portugal']
    
    else: country_names2 = country_names
    
    # Using the Natural Earth feature interface provided by cartopy.
    # You could use a different source, all you need is the geometries.
    shape_records = Reader(natural_earth(resolution='110m',
                                         category='cultural',
                                         name='admin_0_countries')).records()
    geoms = []
    names = []
    for country in shape_records:
        if country.attributes['NAME_LONG'] in country_names2:
            try:
                print('processing',country.attributes['NAME_LONG'])
                geoms += country.geometry
            except TypeError:
                geoms.append(country.geometry)
            for i in range(len(country.geometry)): names.append(country_names)
    return geoms, ccrs.PlateCarree()._as_mpl_transform, names

def get_geometries_land(name0):
    """
    Get an iterable of Shapely geometries corrresponding to given countries.

    """
    # Using the Natural Earth feature interface provided by cartopy.
    # You could use a different source, all you need is the geometries.
    print(name0)
    shape_records = Reader(shpreader.natural_earth(resolution='110m',
                                      category='Physical',
                                      name='land')).records()
    geoms = []
    names=[]
    for country in shape_records:
        try:
            geoms += country.geometry
        except TypeError:
            geoms.append(country.geometry)
        for i in range(len(country.geometry)): 
            names.append(name0)
    return geoms, ccrs.PlateCarree()._as_mpl_transform,names

def get_geometries_ocean(name0):
    """
    Get an iterable of Shapely geometries corrresponding to given countries.

    """
    # Using the Natural Earth feature interface provided by cartopy.
    # You could use a different source, all you need is the geometries.
    print(name0)
    shape_records = Reader(shpreader.natural_earth(resolution='110m',
                                      category='Physical',
                                      name='ocean')).records()
    geoms = []
    names=[]
    for country in shape_records:
        try:
            geoms += country.geometry
        except TypeError:
            geoms.append(country.geometry)
    if name0=='MPac':
        geoms=add_country(geoms,'Taiwan')

    if name0=='SCSPhi':
        geoms=add_country(geoms,'Philippines')
           
    for i in range(len(geoms)): names.append(name0)
    return geoms, ccrs.PlateCarree()._as_mpl_transform,names

def mask_oro(mask,lon,lat,index,threshold):
        #Load the orography
        ORO=pickle.load(gzip.open('orography_Asia.pkl','rb'))
        oro=ORO['oro']
        oro_lon=ORO['lon']
        oro_lat=ORO['lat']
        oro[oro==-32767]=0
        xoro, yoro = np.meshgrid(oro_lon,oro_lat)
        #reshape the orography accordingly to the lon and lat of the mask
        oro2=np.histogram2d(xoro.flatten(),yoro.flatten(),bins=(lon,lat),weights=oro.flatten())[0]
        oro0=np.histogram2d(xoro.flatten(),yoro.flatten(),bins=(lon,lat))[0]
        oro=oro2/oro0
        oro[np.isnan(oro)]=0.    
        ii=oro>threshold
        mask[ii]=index
        return mask

def add_mask(x,y,mask,lab,name,number,limits,oro=None):
    lon=x[0,0:]
    lat=y[0:,0]
    if np.shape(limits)!=(2,2):    
        lonmin=np.min(x); latmin=np.min(y); lonmax=np.max(x) ; latmax=np.max(y)
        #print('no limit latmin',latmin)
    else:
        lonmin=limits[0][0]; latmin=limits[1][0]; lonmax=limits[0][1] ; latmax=limits[1][1]
        #print('limit latmin',latmin)
    if lab=='country':
        geoms, transform, names = get_geometries_countries(name) 
    elif lab=='ocean':
        geoms, transform, names = get_geometries_ocean(name) 
    elif lab=='land':
        geoms, transform, names= get_geometries_land(name) 
    elif lab=='orography':
        print('orography feeded')
    else:
        print('You are not specifying a correct label for mask: country, ocean or land?')
        
    Points = [Point(xp, yp) for xp, yp in zip(x.ravel(), y.ravel())]
    Points = cascaded_union(Points)
    if lab!='orography': 
        print('lengeoms',len(geoms))
        for g in range(len(geoms)):        
            geom = Points.intersection(geoms[g])
            # There might be some errors in the following step 
            # non determined cause (empty intersections?)
            try:
#                if len(geom.geoms)>0:
#                    for points in geom.geoms:
#                        print(points.x,points.y)    
                A= np.array([(points.x, points.y) for points in geom.geoms])[:]
                xx=A[:,0]
                yy=A[:,1]
                bu=np.histogram2d(xx,yy, bins=(lon,lat),  normed=False)[0]
                ii=np.logical_and(np.logical_and(bu.T>0,np.logical_and(x[:-1,:-1]>=lonmin,x[:-1,:-1]<=lonmax)),np.logical_and(y[:-1,:-1]>=latmin,y[:-1,:-1]<=latmax))        
                if np.sum(ii)>0: print('ii',ii.shape,mask.shape,np.sum(ii))
                mask[ii.T] = number
            except:
                #print('exit of add_mask for',name)
                continue
    else:
        mask=mask_oro(mask,lon,lat,number,3800.)        
        
    plt.imshow(mask.T,origin='lower')
    plt.show()
    return mask

def get_and_plot_mask(file):
    #Read the mask (shoule be the pkl file)
    mask0=pickle.load(gzip.open(file,'rb'))
    mask=mask0['mask'].T
    lonc=mask0['lons']
    latc=mask0['lats']
    number_to_name=mask0['regcode_inv']
    limits=mask0['limits']    
    ccreg=mask0['ccreg']
    countries=mask0['countries']
    xl=mask0['xl']
    yl=mask0['yl']
    plot_mask(lonc,latc,mask,number_to_name,limits,ccreg,countries,xl,yl)

#QUESTO E' DECISAMENTE DA FINIRE E RIFINIRE
# This routine does not plot correctly the rolled part of the mask when this happens
# unclear why 

def plot_mask(x,y,mask,number_to_name,limits,ccreg,countries,xl={},yl={},\
   boundaries=None,savefig=None,title=None,groups=None):
    mask=mask.T
    for l in limits.keys():
        if np.shape(limits[l])!=(2,2):    
            lonmin=np.min(x); latmin=np.min(y); lonmax=np.max(x) ; latmax=np.max(y)
        else:
            lonmin=limits[l][0][0]; latmin=limits[l][1][0]; lonmax=limits[l][0][1] ; latmax=limits[l][1][1]    
        if len(xl)==0:        
            xl[l]=(lonmax-lonmin)/2.+lonmin-(lonmax-lonmin)/2.*0.01
            yl[l]=(latmax-latmin)/2.+latmin-(latmax-latmin)/2.*0.01           

    col=[v for k,v in ccreg.items()]
    ind=np.array(list(number_to_name.keys()))
    cmap = colors.ListedColormap(np.array(col)[np.argsort(ind)])
    bounds=ind[np.argsort(ind)]-0.2
    bounds=list(bounds)
    bounds.append(np.max(ind)+1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    #source_range=np.array([[np.min(x),np.max(x)],[np.min(y),np.max(y)]])
    fs = 15
    # it is unclear how the trick with cm_lon works in imshow but it does
    # the web says that it is tricky to plot data accross dateline with cartopy
    # check https://stackoverflow.com/questions/47335851/issue-w-image-crossing-dateline-in-imshow-cartopy
    cm_lon =0
    # guess that we want to plot accross dateline 
    #if source_range[0,1]> 180: cm_lon = 180
    proj = ccrs.PlateCarree(central_longitude=cm_lon)
    fig = plt.figure(figsize=[15,5])
    fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
    ax = plt.axes(projection = proj)
    listm=number_to_name.keys()
    mask2=np.ma.masked_all(mask.shape)
    for ii in listm: 
        mask2[mask==ii]=mask[mask==ii]
    iax = ax.pcolormesh(x,y,mask2,cmap=cmap,norm=norm)
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
    # to eliminate white borders
    plt.xlim(x[0],x[-1])
    plt.ylim(y[0],y[-1])
    if boundaries is not None:
        plt.xlim(boundaries[0],boundaries[1])
        plt.ylim(boundaries[2],boundaries[3])
    #plt.ylim(0,50)
    if groups==None:
        for reg in countries:
            if xl[reg] is not None:
                plt.text(xl[reg],yl[reg],reg, fontsize =14, color = 'k',bbox={'facecolor':'white', 'alpha':0.5, 'pad':5})
    else:
        for gr in groups:
            if xl[gr] is not None:
                grn = gr
                if gr=='Tibet': grn = 'Tibetan Plateau'
                plt.text(xl[gr],yl[gr],grn, fontsize =14, color = 'k',bbox={'facecolor':'white', 'alpha':0.5, 'pad':5})
    if title is not None:
        try:
            plt.title(title,fontsize=20)
        except:
            print('title argument not suitable')
    else:
        plt.title('Region mask',fontsize=20)
    if savefig is not None:
        try:
            plt.savefig(savefig+'.png',dpi=300,bbox_inches='tight')
            plt.savefig(savefig+'.pdf',dpi=300,bbox_inches='tight')
        except:
            print('savefig argument not suitable')
    plt.show()
    return

if __name__ == '__main__':
    main()
