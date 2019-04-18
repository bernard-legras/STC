#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script gathers all the analysis made for a given supertype for the dates from June to August 2017 into a single 
archive that can be used to produce the plots. 
This version is designed for the SAFBox runs with the meanhigh filter in the initialization.
For each supertype, there is an analysis archive for ages truncated at 1 month and for ages truncated at 2 months
Essentially useless since the gathering is now included in plot1Box

Created on Sun Jan 20 19:11:52 2019

@author: Bernard Legras
"""
import pickle,gzip
import socket
import os
import transit as tt
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",choices=["EAD","EAZ","EIZ","EID","EIZ-FULL","EID-FULL"],help="type")
parser.add_argument("-hm","--hmax",type=int,choices=[960,1728],help="max age to be considered (hour)")
parser.add_argument("-v","--vert",choices=["theta","baro"],help="vertical discretization")
parser.add_argument("-g","--globa",type=str,choices=["y","n"],help="global (y) or not (n)")
parser.add_argument("-a","--all",type=str,choices=["All","Allx","All6","All7","All8"],help="selection of time period")

supertype = 'EAD'
hmax = 1728
step = 6
vert = 'theta'
# two choices 'FullAMA' or 'global'
target = 'FullAMA'
water_path = True
All = "-All"

args = parser.parse_args()
if args.type is not None: supertype = args.type
if args.vert is not None: vert = args.vert
if args.hmax is not None: hmax = args.hmax
if args.globa is not None:
        if args.globa =='y': target = 'global'
if 'FULL' not in supertype: target = 'FullAMA'
if args.all is not None: All = "-"+args.all

dates ={"-All":["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],
        "-Allx":["Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],
        "-All6":["Jun-01","Jun-11","Jun-21"],
        "-All7":["Jul-01","Jul-11","Jul-21"],
        "-All8":["Jul-01","Jul-11","Jul-21"]}

run_type = supertype+'-Box-'+vert+'-'+target

# Find the out_dir according to the computer
# ACHTUNG: nosh
if socket.gethostname() == 'gort':
    out_dir = "/dkol/data/STC/STC-FORWBox-meanhigh-OUT"
elif 'ciclad' in socket.gethostname():
    out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"
elif socket.gethostname() == 'satie':
    out_dir = "/limbo/data/STC/STC-FORWBox-meanhigh-OUT"
pile_arch_file = os.path.join(out_dir,'pile-save-stream-'+run_type+All+'-h'+str(hmax)+'.pkl')
    
# Definition of the archive as a new transit class
pile_arch = tt.transit(water_path=water_path,vert=vert,target=target)

for date in dates[All]:
    pile_save_stream = os.path.join(out_dir,'pile-save-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')
    with gzip.open(pile_save_stream,'rb') as f:
        pile = pickle.load(f)
    print(date,np.sum(pile.transit['hist_t']),np.sum(pile.transit['hist_t_vh']),np.sum(pile.transit['hist_t_sh']))
    pile_arch.merge(pile)
    
print(All,np.sum(pile_arch.transit['hist_t']),np.sum(pile_arch.transit['hist_t_vh']),np.sum(pile_arch.transit['hist_t_sh']))
with gzip.open(pile_arch_file,'wb') as f:
    pickle.dump(pile_arch,f,pickle.HIGHEST_PROTOCOL)