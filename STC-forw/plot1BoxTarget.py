# -*- coding: utf-8 -*-
"""
Created on 25 May 2019

Special version of plot1Box to plot the impacts calculated by stat-forw-target
No accumulation is neeed here

@author: Bernard
"""
import gzip, pickle
import transit as tt
import argparse
import numpy as np
import socket
import os

parser = argparse.ArgumentParser()
#parser.add_argument("-y","--year",type=int,help="year")
#parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
parser.add_argument("-t","--type",choices=["EAD","EAZ","EIZ-FULL","EID-FULL"],help="type")
parser.add_argument("-st","--stream",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')
parser.add_argument("-v","--vert",choices=["theta","baro"],help="vertical discretization")
#parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help="quiet (y) or not (n)")
parser.add_argument("-g","--globa",type=str,choices=["y","n"],help="global (y) or not (n)")
parser.add_argument("-hm","--hmax",type=int,choices=[960,1728],help="max age to be considered (hour)")
parser.add_argument("-ht","--hightype",type=str,choices=['mh','vh','sh'],help="high cloud selection")
#parser.add_argument("-a","--all",type=str,choices=["All","Allx","All6","All7","All8"],help="selection of time period")

# default values for the first start day
#year = 2017
#month = 7
#day = 1
supertype = 'EAD'
hmax = 1728
water_path = True
quiet = False
vert = 'theta'
# 'FullAMA' or 'global'
target = 'FullAMA'
hightype = 'sh'
nerd = True
water_path = True
stream = 'Jul-11'

# Choice of what is plotted (these parameters are not passed as arguments yet)
show=True
plot_hist = True
plot_mage = False
plot_mthet = False
plot_mz = False

args = parser.parse_args()
if args.hmax is not None: hmax = args.hmax
if args.type is not None: supertype = args.type
if args.vert is not None: vert = args.vert
#if args.quiet is not None:
#        if args.quiet=='y': quiet=True
#        else: quiet=False
if args.globa is not None:
        if args.globa=='y':target = 'global'
if args.hightype is not None: hightype = args.hightype
if args.stream is not None: stream = args.stream

if socket.gethostname() == 'gort':
    out_dir = "/dkol/data/STC/STC-FORWBox-meanhigh-OUT"
elif 'ciclad' in socket.gethostname():
    out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"
elif socket.gethostname() == 'Graphium':
    out_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-FORWBox-meanhigh-OUT"
elif socket.gethostname() == 'satie':
    out_dir = "/limbo/data/STC/STC-FORWBox-meanhigh-OUT"
else:
    'This program does not run on this computer'
    exit()
    
if hightype == 'mh':
    suffix = ''
    ht = 'MH'
    suffix2 = '-mh'
elif hightype == 'vh':
    suffix = '_vh'
    ht = 'VH'
    suffix2 = '-vh'
elif hightype == 'sh':
    suffix = '_sh'
    ht = 'SH'
    suffix2 = '-sh'
    
run_type = supertype+'-Box-'+vert+'-'+target

pile_sav_stream  =  os.path.join(out_dir,'pile-save-target-'+run_type+'-'+stream+'-h'+str(hmax)+'.pkl')
with gzip.open(pile_sav_stream,'rb') as f:
    pile = pickle.load(f)
    print(stream,np.sum(pile.transit['hist_t']),np.sum(pile.transit['hist_t_sh']))

## Definition of the archive as a new transit class
#pile = tt.transit(water_path=water_path,vert=vert,target=target)
#
## Accumulating the data from the multiple runs
#for date in dates[All]:
#    pile_save_stream = os.path.join(out_dir,'pile-save-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')
#    with gzip.open(pile_save_stream,'rb') as f:
#        pile_arch = pickle.load(f)
#        print(date,np.sum(pile_arch.transit['hist_t']),np.sum(pile_arch.transit['hist_t_vh']),np.sum(pile_arch.transit['hist_t_sh']))
#    pile.merge(pile_arch)
#    del pile_arch
#print(All,np.sum(pile.transit['hist_t']),np.sum(pile.transit['hist_t_vh']),np.sum(pile.transit['hist_t_sh']))
   
# reading the data
#pile_name = os.path.join(out_dir,'pile-save-stream-'+run_type+'.pkl')
#with gzip.open(pile_name,'rb') as f:  
#    pile = pickle.load(f) 

# making averages
pile.complete()

#print('completion performed')
#pile = tt.transit(water_path = pile2.water_path, vert = pile2.vertype, target = pile2.target['type'])
#pile.transit = pile2.transit

# scaling factor for hist_t and hist_s 
# 
# the scaling for hist_t in the volume of each mesh to get an impact density per volume unit in the x,y,theta space
# this factor is the surface of the mesh used in the target space by which the count must be divided
# it is further divided by the width of the layer (500m for baro plots and 5K for theta plots) 
# and multiplied by the interval between outputs, that is 6 h
# since the mesh is of one degree and that the degree was also the unit used in weighting each pixel in transit
# we get to use a numerical factor 6/(5*cos(lat_target))
# the resulting unit is in hour**2 /K for h
# 
# We use the same factor for hist_s since the count is made for each layer in the target space
# but the surface factor is now the surface of the grid in the source space (100 times the satellite image pixel size 
# but we cancel the cos (lat source factor that was introduced in the scaling of the source))
# Lsat: we transform units to day**2 by dividing by 576
ff_t = 6/(5*np.cos(np.radians(pile.target['ycent'])))/576
ff_s = 6/(5*np.cos(np.radians(pile.source['ycent'])))/576
pile.transit['hist_s'+suffix] *= ff_s[np.newaxis,:,np.newaxis]
pile.transit['hist_t'+suffix] *= ff_t[np.newaxis,:,np.newaxis]

names = {"EAD":"ERA5 diabatic","EAZ":"ERA5 kinematic",
         "EIZ":"ERA-I kinematic","EID":"ERA-I diabatic",
         "EIZ-FULL":"ERA-I kinematic","EID-FULL":"ERA-I diabatic"}

if nerd:
    run_name = run_type+'-'+stream
else:
    run_name = names[supertype]+' '+target
    
#%%
# deactivated plots as the barotropic perspective does not present any real interest
if vert == 'baro':
    pass
#    if plot_hist:
#        pile.chart('hist_t_vh',0,txt=run_name+' target 10 km VH',fgp=run_type+'-target-10km-vh',show=show)      
#        pile.chart('hist_t_vh',3,txt=run_name+' target 12 km VH',fgp=run_type+'-target-12km-vh',show=show)
#        pile.chart('hist_t_vh',7,txt=run_name+' target 14 km VH',fgp=run_type+'-target-14km-vh',show=show)
#        pile.chart('hist_t_vh',11,txt=run_name+' target 16 km VH',vmin=0,vmax=6,fgp=run_type+'-target-16km-vh',show=show)
#        print(np.sum(pile.transit['hist_t_vh'][11,:,:]),np.sum(pile.transit['hist_s_vh'][11,:,:]))
#        pile.chart('hist_t_vh',15,txt=run_name+' target 18 km VH',fgp=run_type+'-target-18km-vh',show=show)
#        pile.chart('hist_s_vh',0,txt=run_name+' source of 10 km VH',fgp=run_type+'-source-10km-vh',show=show)
#        pile.chart('hist_s_vh',3,txt=run_name+' source of 12 km VH',fgp=run_type+'-source-12km-vh',show=show)
#        pile.chart('hist_s_vh',7,txt=run_name+' source of 14 km VH',fgp=run_type+'-source-14km-vh',show=show)
#        pile.chart('hist_s_vh',11,txt=run_name+' source of 16 km VH',vmin=0,vmax=40,fgp=run_type+'-source-16km-vh',show=show)
#        pile.chart('hist_s_vh',15,txt=run_name+' source of 18 km VH',fgp=run_type+'-source-18km-vh',show=show)
#        pile.chartv('hist_t_vh',30,txt=run_name+' target 30N VH',vmin=0,vmax=2000,fgp=run_type+'-target-30N-vh',show=show)
#        pile.chartv('hist_t_vh',35,txt=run_name+' target 35N VH',vmin=0,vmax=2000,fgp=run_type+'-target-35N-vh',show=show)
#        pile.chartv('hist_t_vh',25,txt=run_name+' target 25N VH',vmin=0,vmax=2000,fgp=run_type+'-target-25N-vh',show=show)
#    if plot_mage:        
#        pile.chartv('mage_t',30,txt=run_name+' mage target 30N',fgp=run_type+'-mage-target-30N',show=show)
#        pile.chartv('mage_t',35,txt=run_name+' mage target 35N',fgp=run_type+'-mage-target-35N',show=show)
#        pile.chartv('mage_t',25,txt=run_name+' mage target 25N',fgp=run_type+'-mage-target-25N',show=show)
#        pile.chart('mage_t',11,txt=run_name+' mage target 16 km',fgp=run_type+'-mage-target-16km',show=show)
#        pile.chart('mage_t',7,txt=run_name+' mage target 14 km',fgp=run_type+'-mage-target-14km',show=show)
#        pile.chart('mage_t',15,txt=run_name+' mage target 18 km',fgp=run_type+'-mage-target-18km',show=show)   
#    if plot_mthet:
#        pile.chart('mthet_s',11,vmin=340,vmax=380,txt=run_name+' mthet source 16 km',fgp=run_type+'-mthet-source-16km',show=show)
#        pile.chart('mthet_s',7,vmin=340,vmax=380,txt=run_name+' mthet source 14 km',fgp=run_type+'-mthet-source-14km',show=show)
#        pile.chart('mthet_s',15,vmin=350,vmax=390,txt=run_name+' mthet source 18 km',fgp=run_type+'-mthet-source-18km',show=show)       
#        pile.chart('mdthet_s',11,txt=run_name+' mdthet source 16 km',fgp=run_type+'-mdthet-source-16km',show=show)
#        pile.chart('mdthet_s',7,txt=run_name+' mdthet source 14 km',fgp=run_type+'-mdthet-source-14km',show=show)
#        pile.chart('mdthet_s',15,txt=run_name+' mdthet source 18 km',fgp=run_type+'-mdthet-source-18km',show=show) 
#    if plot_mz:
#        pile.chart('mz_s',11,vmin=13,vmax=16,txt=run_name+' mz source 16 km',fgp=run_type+'-mz-source-16km',show=show)
#        pile.chart('mz_s',7,vmin=13,vmax=15,txt=run_name+' mz source 14 km',fgp=run_type+'-mz-source-14km',show=show)
#        pile.chart('mz_s',15,vmin=14,vmax=17,txt=run_name+' mz source 18 km',fgp=run_type+'-mz-source-18km',show=show)       
#        pile.chart('mdz_s',11,vmin=0,vmax=3,txt=run_name+' mdz source 16 km',fgp=run_type+'-mdz-source-16km',show=show)
#        pile.chart('mdz_s',7,vmin=-1,vmax=1,txt=run_name+' mdz source 14 km',fgp=run_type+'-mdz-source-14km',show=show)
#        pile.chart('mdz_s',15,vmin=0,vmax=3,txt=run_name+' mdz source 18 km',fgp=run_type+'-mdz-source-18km',show=show) 

#%%
elif vert == 'theta':
    if plot_hist:
        pile.chart('hist_t'+suffix,3,txt=run_name+' target 340 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-340K'+suffix2,show=show)
        pile.chart('hist_t'+suffix,5,txt=run_name+' target 350 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-350K'+suffix2,show=show)
        pile.chart('hist_t'+suffix,7,txt=run_name+' target 360 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-360K'+suffix2,show=show)
        pile.chart('hist_t'+suffix,9,txt=run_name+' target 370 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-370K'+suffix2,show=show)
        pile.chart('hist_t'+suffix,11,txt=run_name+' target 380 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-380K'+suffix2,show=show)
        #print('380K',np.sum(pile.transit['hist_t'+suffix][11,:,:]),np.sum(pile.transit['hist_s'+suffix][11,:,:]))
        pile.chart('hist_t'+suffix,13,txt=run_name+' target 390 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-390K'+suffix2,show=show)
        pile.chart('hist_t'+suffix,15,txt=run_name+' target 400 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-400K'+suffix2,show=show)
        pile.chart('hist_t'+suffix,17,txt=run_name+' target 420 K '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-420K'+suffix2,show=show)    
        pile.chart('hist_s'+suffix,3,txt=run_name+' source of 340 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-340K'+suffix2,show=show)
        pile.chart('hist_s'+suffix,5,txt=run_name+' source of 350 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-350K'+suffix2,show=show)
        pile.chart('hist_s'+suffix,7,txt=run_name+' source of 360 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-360K'+suffix2,show=show)
        pile.chart('hist_s'+suffix,9,txt=run_name+' source of 370 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-370K'+suffix2,show=show)
        pile.chart('hist_s'+suffix,11,txt=run_name+' source of 380 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-380K'+suffix2,show=show)
        pile.chart('hist_s'+suffix,13,txt=run_name+' source of 390 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-390K'+suffix2,show=show)
        pile.chart('hist_s'+suffix,15,txt=run_name+' source of 400 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-400K'+suffix2,show=show)
        pile.chart('hist_s'+suffix,17,txt=run_name+' source of 420 K '+ht+' : conv source density (day$^2$ K$^{-1}$)',fgp=run_name+'-source-420K'+suffix2,show=show)
        pile.chartv('hist_t'+suffix,30,txt=run_name+' target 30N '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-30N'+suffix2,show=show)
        pile.chartv('hist_t'+suffix,35,txt=run_name+' target 35N '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-35N'+suffix2,show=show)
        pile.chartv('hist_t'+suffix,25,txt=run_name+' target 25N '+ht+' : conv impact density (day$^2$ K$^{-1}$)',fgp=run_name+'-target-25N'+suffix2,show=show)   
    #%%
    if plot_mage:
        pile.chartv('mage_t'+suffix,30,txt=run_name+' mean age target 30N '+ht+' (day)',fgp=run_name+'-mage-target-30N'+suffix2,show=show)
        pile.chartv('mage_t'+suffix,35,txt=run_name+' mean age target 35N '+ht+' (day)',fgp=run_name+'-mage-target-35N'+suffix2,show=show)
        pile.chartv('mage_t'+suffix,25,txt=run_name+' mean age target 25N '+ht+' (day)',fgp=run_name+'-mage-target-25N'+suffix2,show=show)
        pile.chart('mage_t'+suffix,3,back_field='hist_t'+suffix,cumsum=True,txt=run_name+' mean age target 340 K '+ht+' (day)',fgp=run_name+'-mage-target-340K'+suffix2,show=show)
        pile.chart('mage_t'+suffix,5,back_field='hist_t'+suffix,cumsum=True,txt=run_name+' mean age target 350 K '+ht+' (day)',fgp=run_name+'-mage-target-350K'+suffix2,show=show) 
        pile.chart('mage_t'+suffix,7,back_field='hist_t'+suffix,cumsum=True,txt=run_name+' mean age target 360 K '+ht+' (day)',fgp=run_name+'-mage-target-360K'+suffix2,show=show)
        pile.chart('mage_t'+suffix,9,back_field='hist_t'+suffix,cumsum=True,txt=run_name+' mean age target 370 K '+ht+' (day)',fgp=run_name+'-mage-target-370K'+suffix2,show=show)
        pile.chart('mage_t'+suffix,11,back_field='hist_t'+suffix,cumsum=True,txt=run_name+' mean age target 380 K '+ht+' (day)',fgp=run_name+'-mage-target-380K'+suffix2,show=show)
        pile.chart('mage_t'+suffix,13,back_field='hist_t'+suffix,cumsum=True,txt=run_name+' mean age target 390 K '+ht+' (day)',fgp=run_name+'-mage-target-390K'+suffix2,show=show) 
        pile.chart('mage_t'+suffix,15,back_field='hist_t'+suffix,cumsum=True,txt=run_name+' mean age target 400 K '+ht+' (day)',fgp=run_name+'-mage-target-400K'+suffix2,show=show)    
    #%%
        pile.chart('mage_s'+suffix,3,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean age source 340 K '+ht+' (day)',fgp=run_name+'-mage-source-340K'+suffix2,show=show)
        pile.chart('mage_s'+suffix,5,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean age source 350 K '+ht+' (day)',fgp=run_name+'-mage-source-350K'+suffix2,show=show) 
        pile.chart('mage_s'+suffix,7,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean age source 360 K '+ht+' (day)',fgp=run_name+'-mage-source-360K'+suffix2,show=show)
        pile.chart('mage_s'+suffix,9,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean age source 370 K '+ht+' (day)',fgp=run_name+'-mage-source-370K'+suffix2,show=show)
        pile.chart('mage_s'+suffix,11,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean age source 380 K '+ht+' (day)',fgp=run_name+'-mage-source-380K'+suffix2,show=show)
        pile.chart('mage_s'+suffix,13,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean age source 390 K '+ht+' (day)',fgp=run_name+'-mage-source-390K'+suffix2,show=show) 
        pile.chart('mage_s'+suffix,15,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean age source 400 K '+ht+' (day)',fgp=run_name+'-mage-source-400K'+suffix2,show=show)    
    #%%
    
    #%%
    if plot_mthet:    
        pile.chart('mthet_s'+suffix,7,back_field='hist_s'+suffix,cumsum=True,vmin=350,vmax=370,txt=run_name+' mean θ source of 360 K '+ht+' (K)',fgp=run_name+'-mthet-source-360K'+suffix2,show=show)
        pile.chart('mthet_s'+suffix,11,back_field='hist_s'+suffix,cumsum=True,vmin=350,vmax=380,txt=run_name+' mean θ source of 380 K '+ht+' (K)',fgp=run_name+'-mthet-source-380K'+suffix2,show=show)
        pile.chart('mthet_s'+suffix,15,back_field='hist_s'+suffix,cumsum=True,vmin=360,vmax=380,txt=run_name+' mean θ source of 400 K '+ht+' (K)',fgp=run_name+'-mthet-source-400K'+suffix2,show=show)
        pile.chart('mdthet_s'+suffix,7,back_field='hist_s'+suffix,cumsum=True,vmin=0,vmax=20,txt=run_name+' mean Δθ source of 360 K '+ht+' (K)',fgp=run_name+'-mdthet-source-360K'+suffix2,show=show)
        pile.chart('mdthet_s'+suffix,11,back_field='hist_s'+suffix,cumsum=True,vmin=0,vmax=20,txt=run_name+' mean Δθ source of 380 K '+ht+' (K)',fgp=run_name+'-mdthet-source-380K'+suffix2,show=show)
        pile.chart('mdthet_s'+suffix,15,back_field='hist_s'+suffix,cumsum=True,vmin=5,vmax=35,txt=run_name+' mean Δθ source of 400 K '+ht+' (K)',fgp=run_name+'-mdthet-source-400K'+suffix2,show=show)   
    if plot_mz:
        pile.chart('mz_s'+suffix,7,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean z source 360 K '+ht+' (km)',fgp=run_name+'-mz-source-360K'+suffix2,show=show)
        pile.chart('mz_s'+suffix,11,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean z source 380 K '+ht+' (km)',fgp=run_name+'-mz-source-380K'+suffix2,show=show)
        pile.chart('mz_s'+suffix,15,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean z source 400 K '+ht+' (km)',fgp=run_name+'-mz-source-400K'+suffix2,show=show)
        pile.chart('mdz_s'+suffix,7,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean z source 360 K '+ht+' (km)',fgp=run_name+'-mdz-source-360K'+suffix2,show=show)
        pile.chart('mdz_s'+suffix,11,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean z source 380 K '+ht+' (km)',fgp=run_name+'-mdz-source-380K'+suffix2,show=show)
        pile.chart('mdz_s'+suffix,15,back_field='hist_s'+suffix,cumsum=True,txt=run_name+' mean z source 400 K '+ht+' (km)',fgp=run_name+'-mdz-source-400K'+suffix2,show=show)
