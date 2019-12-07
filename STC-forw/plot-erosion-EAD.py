#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script shows the erosion rate of the impact with altitude 

It uses diag files form the synthetic pile object generated in plot1Box
 
Created on Sun Sep 15 2019

@author: Bernard Legras
"""
import numpy as np
#from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
from os.path import join
from scipy.stats.stats import pearsonr 

#with gzip.open('EID-FULL-Box-theta-global.pkl','rb') as f:
#    [hei_g, tei_g] = pickle.load(f)
#with gzip.open('EID-FULL-Box-theta-FullAMA.pkl','rb') as f:
#    [hei_f, tei_f] = pickle.load(f)
with gzip.open('EAD-Box-theta-FullAMA.pkl','rb') as f:
    [hea_f, tea_f] = pickle.load(f)
    
forw_dir = '.'
dpi = 300

#hei_gr = hei_g[:,90:140,169:339]

#hea_f = np.reshape(hea_f,(20,8500))
#hei_f = np.reshape(hei_f,(20,8500))
#hei_gr = np.reshape(hei_gr,(20,8500))

vcent=np.array([325., 330., 335., 340., 345., 350., 355., 360., 365., 370., 375.,
       380., 385., 390., 395., 400., 405., 410., 415., 420.])

#%%
ss = np.sum(hea_f ,axis=(1,2))
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(1, 1, 1)
iax=ax.plot(ss,vcent)
ax.set_xscale('log')
ax.title('impact')
plt.show()