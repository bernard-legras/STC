#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 19:03:20 2017

Demonstration and test of subgridding

@author: legras
"""
import geosat
#import pickle,gzip
#import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#%%
date=datetime(year=2022,month=1,day=27,hour=15)
# Himawari read and plot
ah=geosat.Himawari(date)
gg=geosat.GeoGrid('HimFull')
ah._get_IR0()
ph=geosat.SatGrid(ah,gg)
ph._sat_togrid('IR0')
ah.show('IR0')
ph.chart('IR0',txt='BT himawari')

# exemple of subgridding and plot with this subgrid
subgg=gg.subgrid([145,195,-30,-15])
print(subgg.box_range)
print(subgg.shapexy)
ph.chart('IR0',txt='subgrid test',subgrid=subgg)

#%%
# tran = ccrs.PlateCarree(central_longitude=1
# proj = tran
# # geostationary plot
# ax = plt.subplot(2,1,2,projection = proj)
# #ax.set_extent(ext,ccrs.PlateCarree())
# ax.set_extent(subbox,ccrs.PlateCarree())
# if  GEOOK:
#     ax = ph.chart('Ash',txt='Ash composit '+utc.strftime('%Y-%m-%d %H:%M'),axf=ax,subgrid=subgg,show=False,cm_lon=cm)
# ax.coastlines('50m')
# gl = ax.gridlines(draw_labels=True,xlocs=xlocs,
#              linewidth=2, color='gray', alpha=0.5, linestyle='--')
# gl.top_labels = False
# gl.right_labels = False
# gl.bottom_labels = False
# gl.left_labels = False
# # to make sure, plot it twice
# ax.plot(lons1+cm,lats1,'k')
# ax.plot(lons1-cm,lats1,'k')


#%%
fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(11,8),
                       subplot_kw={"projection":ccrs.PlateCarree(central_longitude=180)})
fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
geogrid = subgg
field = ph.var['IR0'][geogrid.corner[1]:geogrid.corner[1]+geogrid.box_biny,
                      geogrid.corner[0]:geogrid.corner[0]+geogrid.box_binx]
#x.set_extent()
#ax.set_extent((145,195,-30,-15),crs=ccrs.PlateCarree())
ax0 = ph.chart('IR0',axf=ax[0],show=False,cm_lon=180,subgrid=subgg,bottom=False)
ax1 = ph.chart('IR0',axf=ax[1],show=False,cm_lon=180,subgrid=subgg)
plt.show()

#%%
ah._mk_Ash()
ph._sat_togrid('Ash')
ph.chart('Ash',txt='subgrid test',subgrid=subgg)
