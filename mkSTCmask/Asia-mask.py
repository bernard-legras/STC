# -*- coding: utf-8 -*-
"""

Script to generate the figure of the Asia mask

Created on Tue Nov  5 00:16:02 2019

@author: berna
"""
import pickle,gzip
from group import group
from Nnewmask import plot_mask
from collections import OrderedDict

#We open the existing mask
mask0=pickle.load(gzip.open('MaskCartopy2-STCfine.pkl','rb'))
mask=mask0['mask'].T
lonc=mask0['lons']
latc=mask0['lats']
name_to_number=mask0['regcode']
number_to_name=mask0['regcode_inv']
limits=mask0['limits']    
ccreg=mask0['ccreg']
countries=mask0['countries']
xl=mask0['xl']
yl=mask0['yl']


#We define the colors, the labels and the position for the new defined groups

groups=['Land','Ocean','Tibet']
ccreg_new = OrderedDict([('Land','coral'),('Tibet', 'lightskyblue'),('Ocean', 'steelblue')])
xl_new = OrderedDict([('Land',80),('Tibet', 80),('Ocean', xl['SCSPhi'])])
yl_new = OrderedDict([('Land',22),('Tibet', 32),('Ocean', yl['SCSPhi'])])

for reg in ccreg.keys(): ccreg[reg]='lightgrey'
for reg in xl.keys(): xl[reg]=None
for reg in yl.keys(): yl[reg]=None
for gr in groups:
    for reg in group[gr]:
        #mask[mask==name_to_number[j]]=1
        number_to_name[name_to_number[reg]]=gr
        name_to_number[gr]=name_to_number[reg]
        countries[countries==reg]=gr
        ccreg[reg]=ccreg_new[gr]        
        ccreg[gr]=ccreg_new[gr]        
        xl[reg]=xl_new[gr]
        yl[reg]=yl_new[gr]
        xl[gr]=xl_new[gr]
        yl[gr]=yl_new[gr]

#We plot it
plot_mask(lonc,latc,mask,number_to_name,limits,ccreg,countries,xl,yl,
          boundaries=[40,160,0,50],groups=groups,title='Asia mask',
          savefig='AsiaMask')

