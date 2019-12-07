#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script generates a comparison of EID FULL, EID, EAD for the sake of justifying to restrain calculations in the FullAMA domain 
to study AMA confinement.
Generates figure 0b of CONF 2019

It uses diag files form the synthetic pile objevt generated in plot1Box
 
Created on Sat Sep 14 18:42:51 2019

@author: Bernard Legras
"""
import numpy as np
#from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
from os.path import join
from scipy.stats.stats import pearsonr 

with gzip.open('EID-FULL-Box-theta-global.pkl','rb') as f:
    [hei_g, tei_g] = pickle.load(f)
with gzip.open('EID-FULL-Box-theta-FullAMA.pkl','rb') as f:
    [hei_f, tei_f] = pickle.load(f)
with gzip.open('EAD-Box-theta-FullAMA.pkl','rb') as f:
    [hea_f, tea_f] = pickle.load(f)
    
forw_dir = '.'
dpi = 300

figsav = False

hei_gr = hei_g[:,90:140,169:339]

hea_f = np.reshape(hea_f,(20,8500))
hei_f = np.reshape(hei_f,(20,8500))
hei_gr = np.reshape(hei_gr,(20,8500))

cei = np.empty(20)
cea = np.empty(20)
rsei = np.empty(20)
rmei = np.empty(20)
rsea  = np.empty(20)
rmea = np.empty(20)

vcent=np.array([325., 330., 335., 340., 345., 350., 355., 360., 365., 370., 375.,
       380., 385., 390., 395., 400., 405., 410., 415., 420.])
#%%
for i in range(20):
    cei[i] = pearsonr(hei_f[i,:],hei_gr[i,:])[0]
    cea[i] = pearsonr(hei_f[i,:],hea_f[i,:])[0]
    rsei[i] = np.sum(hei_f[i,:])/np.sum(hei_gr[i,:]) 
    rsea[i] = np.sum(hea_f[i,:])/np.sum(hei_f[i,:])
    rmei[i] = np.max(hei_f[i,:])/np.max(hei_gr[i,:])
    rmea[i] = np.max(hea_f[i,:])/np.max(hei_f[i,:])  
    print('{:.0f}  {:.1f}  {:.1f}  {:9.1f}  {:9.1f}  {:9.3f}   {:9.3f}   {:9.3f}   {:9.3f}'.format(
            vcent[i], 100*cei[i],100*cea[i], np.sum(hei_f[i,:]), np.sum(hei_gr[i,:]), rsei[i], rsea[i], rmei[i], rmea[i]))
    
#%%
fs=16

fig,ax = plt.subplots(figsize=(6,6))
iax = ax.plot(cei,vcent,'k',cea,vcent,'--k', rsei,vcent,'b',rsea,vcent,'--b', rmei,vcent,'r',rmea,vcent,'--r',linewidth=5)
ax.tick_params(labelsize=fs) 
ax.set_ylim((330,420))
ax.set_xlim((0.2,1.5))
ax.set_ylabel('Potential temperature (K)',fontsize=fs)
#ax.legend([r'$\rho(global,AMA)$',r'toto',
#           r'$\frac{\Sigma \mathrm{AMA}}{\Sigma global}$',r'toto',
#           r'toto',r'toto'],
#           fontsize=14,loc='upper right',)
ax.legend([r'$\rho~(\mathrm{global},\mathrm{AMA})$',r'$\rho~(\mathrm{ERA5},\mathrm{ERAI})$',
           r'${\Sigma~\mathrm{AMA}}~/~{\Sigma~\mathrm{global}}$',r'${\Sigma~\mathrm{ERA5}}~/~{\Sigma~\mathrm{ERAI}}$',
           r'${\max~\mathrm{AMA}}~/~{\max~\mathrm{global}}$',r'${\max~\mathrm{ERA5}}~/~{\max~\mathrm{ERAI}}$'],
           fontsize=16,bbox_to_anchor=(0.7,1))
plt.title('Correlations and ratios',fontsize=fs)
if figsav:
    plt.savefig(join(forw_dir,'figs','correl-ratios-ACP.png'),bbox_inches='tight',dpi=dpi)
    plt.savefig(join(forw_dir,'figs','correl-ratios-ACP.pdf'),bbox_inches='tight',dpi=dpi)
plt.show()

#%%    
plt.plot(hei_f[11,:],hei_gr[11,:],'+')
plt.show()
plt.plot(hei_f[11,:],hea_f[11,:],'+')
plt.show()
#%%
plt.hist2d(hei_f[11,:],hei_gr[11,:])
plt.show()
H2,fedges,gedges = np.histogram2d(hei_f[11,:],hei_gr[11,:])
plt.pcolormesh(fedges,gedges,np.log(H2.T))
plt.show()