# -*- coding: utf-8 -*-
"""
Test of STCdata

Created on Sun Oct 15 11:29:25 2017

@author: Bernard Legras
"""
#from STCdata import STCdata, readAMES1001
import STCdata
from datetime import datetime
#file = 'database\\tdc\\170804_1_tdc.nas'
#tdc = STCdata(file)
#tdc = STCdata.STCdata(file)

# exclude Mainz instruments as files are only templates for now
instruments = ['ucse','tdc','chiwis','flash','fish','amica','copas',\
               'hagar','uhsas','mas']
date = datetime(2017,8,2)

dat={}
for inst in instruments:
    print('reading ',inst)
    if inst in ['pip','cip','ccp-cdp','ccp-cipgs']:
        dat[inst] = STCdata.STCinst(inst,date,fv=True)
    else:
        dat[inst] = STCdata.STCinst(inst,date)

#date = datetime(2017,8,8)
#for inst in ['pip','cip']:
#    dat[inst] = STCdata.STCinst(inst,date,fv=True)