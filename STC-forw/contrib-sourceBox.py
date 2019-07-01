#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

From the output of stat-forw-uniq

Calculate the contributions of the defined regions to the source distribution as
defined in transit in the field 'hist-s'.
Plot a comparison of the vertical distribution of the contributions as a fonction of
the altitude of the the target level (in pressure or potential temperature).
In this version we further calculate also accumulation for the age and the altitude

contribs cumulates the hist_s histogram for each run-type 
(that is combination of supertype, target and vert)(although vert is only theta) 
and for each region.

age, thet, z, dthet, dz  do the same for age_s, thet_s, z_s, dthet_s, dz_s

These accumulations are further calculated for groups of regions (which may differ
from those used in diag codes). Should be harmonized.

Print tables of cumulated impact per region and per level

Plots (with target impact altitude as vertical axis) AML, AO, and Tibetan Plateau
an for a given hmax

1) Cumulated convective impact at target isentropic levels (in linear scale)

1bis) Cumulated convective impact at target isentropic levels (in log scale)

1spe) Cumulated impact for All region and EAD between 370K and 420K
and print of the exponential fit coefficients 

2) Mean age with respect to convection

3) Mean potential temperature of the convective source

4) Mean barometric altitude of the convective source

5) Mean potential temperature displacement from the convective source

6) Mean barometric altitude displacement from the convective source

Note: uses MaskCartopy2-STCforw.pkl instead of MaskCartopy2-STCforwfine.pkl 
as in more recent code as this mask fits the resolution of the stat-forw diagnostics. 

Created on Thu Apr  5 18:28:02 2018

@author: Bernard Legras
"""

import os
import pickle, gzip
import numpy as np
import matplotlib.pyplot as plt
import transit as tt
import socket

mask = pickle.load(gzip.open('MaskCartopy2-STCforw.pkl','rb'))
regcode = mask['regcode']
nb_reg = len(regcode)
print(regcode)

AML_list = ['IndianSub','NorthChina','SouthChina','Pen','Bangladesh','Pakistan',
            'JapanKorea']
AO_list  = ['BoB','SCSPhi','IndianOcean','Indonesia','WestPacific','MidPacific']

#%%
# preliminary: estimation of the surface of the region in degree square units

coslat = np.squeeze(np.cos(np.radians(mask['lats'])))
surfreg = {}
for reg,k in regcode.items():
    surf = np.ones(mask['mask'].shape) * coslat[:,np.newaxis]
    surfreg[reg] = np.sum(surf[mask['mask']==k])
    if surfreg[reg] ==0: 
        surfreg[reg] = 1
surfreg['All'] = 0
for reg in regcode:
    surfreg['All'] += surfreg[reg]
surfreg['AO'] = 0
surfreg['AML'] = 0
for reg in AO_list:
    surfreg['AO'] += surfreg[reg]
for reg in AML_list:
    surfreg['AML'] += surfreg[reg]

#%% obsolete
# scaling factor for hist_t and hist_s
# degree square * 10^-9 / layer thickness
# unit 10^9 km h for baro
# unit 10^9 km^2 h K^-1 for theta 
#ff={'baro':1/0.5,'theta':1/5}

#%%
figsave = False

water_path = True
vert = 'theta'
# parameters to vary
hmax = 1728
# hightype = '' for meanhigh (mh)
# hightype = '_sh' for silviahigh (sh)
hightype = '_sh'
#hightype = ''
ss = hightype
All = '-All'
ss2d = {'_sh':['_sh',' SH'],'':['_mh',' MH']} 
ss2 = ss2d[hightype][1]
ss3 = ss2d[hightype][0]

dates ={"-All":["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],
        "-Allx":["Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],
        "-All6":["Jun-01","Jun-11","Jun-21"],
        "-All7":["Jul-01","Jul-11","Jul-21"],
        "-All8":["Aug-01","Aug-11","Aug-21"]}

# Find the main work dir according to the computer
if socket.gethostname() == 'gort':
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif 'ciclad' in socket.gethostname():
    out_dir = '/data/legras/STC/STC-FORWBox-meanhigh-OUT'
elif 'satie' in socket.gethostname():
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif 'Graphium' in socket.gethostname():
    out_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-FORWBox-meanhigh-OUT" 

contribs = {}
#impact = {}
age = {}
thet = {}
z ={}
dthet = {}
dz = {}
#%%
for supertype in ["EAD","EAZ","EIZ-FULL","EID-FULL"]:
    for target in ['FullAMA','global']:
        if (target == 'global') & ~('FULL' in supertype): continue
        for vert in ['theta',]:
            run_type = supertype+'-Box-'+vert+'-'+target
            pile_name = os.path.join(out_dir,'pile-save-stream-'+run_type+All+'-h'+str(hmax)+'.pkl')
            try:
                with gzip.open(pile_name,'rb') as f:
                    pile = pickle.load(f)
            except:
                print('cannot open ',pile_name)
                print('building')
                pile = tt.transit(water_path=water_path,vert=vert,target=target)
                for date in dates[All]:
                    pile_save_stream = os.path.join(out_dir,'pile-save-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')
                    with gzip.open(pile_save_stream,'rb') as f:
                        pile_arch = pickle.load(f)
                        print(date,np.sum(pile_arch.transit['hist_t']),np.sum(pile_arch.transit['hist_t_vh']),np.sum(pile_arch.transit['hist_t_sh']))
                        pile.merge(pile_arch)
                    del pile_arch
            print('processing ',run_type)
            
            # it is appropriate to use H_s here since we want to link the impact to
            # the source region
            # Notice that the sum of H_t and H_s on a given level are equal
            
            H_s = pile.transit['hist_s'+ss]
            #H_t = pile.transit['hist_t'+ss]
            age_s = pile.transit['totage_s'+ss]
            thet_s = pile.transit['totthet_s'+ss] 
            z_s = pile.transit['totz_s'+ss]
            dthet_s = pile.transit['totdthet_s'+ss] 
            dz_s = pile.transit['totdz_s'+ss]

            dims = H_s.shape
            nlev = dims[0]
            contribs[run_type] = {}
            #impact[run_type] = {}
            age[run_type] = {}
            thet[run_type] = {}
            z[run_type] = {}
            dthet[run_type] = {}
            dz[run_type] = {}
            for reg,k in regcode.items():
                contribs[run_type][reg] = np.empty(nlev)
                age[run_type][reg] = np.empty(nlev)
                thet[run_type][reg] = np.empty(nlev)
                z[run_type][reg] = np.empty(nlev)
                dthet[run_type][reg] = np.empty(nlev)
                dz[run_type][reg] = np.empty(nlev)
                for lev in np.arange(nlev):
                    contribs[run_type][reg][lev] = np.sum(H_s[lev,mask['mask'] == k])
                    age[run_type][reg][lev] = np.sum(age_s[lev,mask['mask'] == k])
                    thet[run_type][reg][lev] = np.sum(thet_s[lev,mask['mask'] == k])
                    z[run_type][reg][lev] = np.sum(z_s[lev,mask['mask'] == k])
                    dthet[run_type][reg][lev] = np.sum(dthet_s[lev,mask['mask'] == k])
                    dz[run_type][reg][lev] = np.sum(dz_s[lev,mask['mask'] == k])    
                #impact[run_type][reg] = contribs[run_type][reg] * ff[vert] /surfreg[reg]
            for supereg in ['All','AML','AO']:
                contribs[run_type][supereg] = np.zeros(nlev)
                age[run_type][supereg] = np.zeros(nlev)
                thet[run_type][supereg] = np.zeros(nlev)
                z[run_type][supereg] = np.zeros(nlev)
                dthet[run_type][supereg] = np.zeros(nlev)
                dz[run_type][supereg] = np.zeros(nlev)
            for reg in regcode:
                contribs[run_type]['All'] += contribs[run_type][reg]
                age[run_type]['All'] += age[run_type][reg]
                thet[run_type]['All'] += thet[run_type][reg]
                z[run_type]['All'] += z[run_type][reg]
                dthet[run_type]['All'] += dthet[run_type][reg]
                dz[run_type]['All'] += dz[run_type][reg]
            for reg in AML_list:
                contribs[run_type]['AML'] += contribs[run_type][reg]
                age[run_type]['AML'] += age[run_type][reg]
                thet[run_type]['AML'] += thet[run_type][reg]
                z[run_type]['AML'] += z[run_type][reg]
                dthet[run_type]['AML'] += dthet[run_type][reg]
                dz[run_type]['AML'] += dz[run_type][reg]
            for reg in AO_list:
                contribs[run_type]['AO'] += contribs[run_type][reg]
                age[run_type]['AO'] += age[run_type][reg]
                thet[run_type]['AO'] += thet[run_type][reg]
                z[run_type]['AO'] += z[run_type][reg]
                dthet[run_type]['AO'] += dthet[run_type][reg]
                dz[run_type]['AO'] += dz[run_type][reg]
            contribs[run_type]['Monsoon'] = contribs[run_type]['AML']+contribs[run_type]['AO']+contribs[run_type]['TibetanPlateau']
            age[run_type]['Monsoon'] = age[run_type]['AML']+age[run_type]['AO']+age[run_type]['TibetanPlateau']
            thet[run_type]['Monsoon'] = thet[run_type]['AML']+thet[run_type]['AO']+thet[run_type]['TibetanPlateau']
            z[run_type]['Monsoon'] = z[run_type]['AML']+z[run_type]['AO']+z[run_type]['TibetanPlateau']
            dthet[run_type]['Monsoon'] = dthet[run_type]['AML']+dthet[run_type]['AO']+dthet[run_type]['TibetanPlateau']
            dz[run_type]['Monsoon'] = dz[run_type]['AML']+dz[run_type]['AO']+dz[run_type]['TibetanPlateau']
            #impact[run_type]['AML']  =  contribs[run_type]['AML'] * ff[vert] /surfreg['AML']
            #impact[run_type]['AO']  =  contribs[run_type]['AO'] * ff[vert] /surfreg['AO']           
            for reg in contribs[run_type].keys():
                age[run_type][reg] /=  contribs[run_type][reg]
                thet[run_type][reg] /=  contribs[run_type][reg]
                z[run_type][reg] /=  contribs[run_type][reg]                          
                dthet[run_type][reg] /=  contribs[run_type][reg]
                dz[run_type][reg] /=  contribs[run_type][reg]
#%%                      
with open('contribsBox'+All+'-h'+str(hmax)+ss3+'.pkl','wb') as f:                        
    pickle.dump([contribs,age,thet,z,dthet,dz],f)
    
#%% Set TeX to write the labels and annotations
plt.rc('text', usetex=True)
plt.rc('font', family='serif') 

baro = np.array([ 10.25,  10.75,  11.25,  11.75,  12.25,  12.75,  13.25,  13.75,
        14.25,  14.75,  15.25,  15.75,  16.25,  16.75,  17.25,  17.75,
        18.25,  18.75,  19.25,  19.75])
theta= np.array([327.5, 332.5, 337.5, 342.5, 347.5, 352.5, 357.5, 362.5, 367.5,
       372.5, 377.5, 382.5, 387.5, 392.5, 397.5, 402.5, 407.5, 412.5,417.5, 422.5])
       
# parameters for plots
# resolution of the png
dpi = 300
                    
#%% Plots for potential temperature target (fig 1)
    
#rts = ['EAD-Box-theta-FullAMA-All-h'+str(hmax),
#       'EID-FULL-Box-theta-FullAMA-All-h'+str(hmax),
#       'EID-FULL-Box-theta-global-All-h'+str(hmax),
#       'EAZ-Box-theta-FullAMA-All-h'+str(hmax),
#       'EIZ-FULL-Box-theta-FullAMA-All-h'+str(hmax),
#       'EIZ-FULL-Box-theta-global-All-h'+str(hmax)]
#
#rts = ['EAD-Box-theta-FullAMA',
#       'EID-FULL-Box-theta-FullAMA',
#       'EID-FULL-Box-theta-global',
#       'EAZ-Box-theta-FullAMA',
#       'EIZ-FULL-Box-theta-FullAMA',
#       'EIZ-FULL-Box-theta-global']
#
#fig=plt.figure(figsize=(11,6))
#fig.suptitle('Convective impact density at target isentropic levels [h'+str(hmax)+ss2+']',fontsize=18)
#ax=plt.subplot(1,3,1)
#ax.plot(impact[rts[0]]['AML'],theta,'k',impact[rts[1]]['AML'],theta,'b',impact[rts[2]]['AML'],theta,'r',
#         impact[rts[3]]['AML'],theta,'g',impact[rts[4]]['AML'],theta,'c',impact[rts[5]]['AML'],theta,'m',
#         linewidth=4)
#plt.ylim((330,400))
#ax.tick_params(labelsize=16)
##plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
##            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
##            fontsize=16,loc='upper left')
#plt.title(r'Land Asia no Tib. Pl.',fontsize=16)
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
#plt.xlabel(r'impact (hour K$^{-1}$)',fontsize=16)
#ax=plt.subplot(1,3,2)
#ax.plot(impact[rts[0]]['AO'],theta,'k',impact[rts[1]]['AO'],theta,'b',impact[rts[2]]['AO'],theta,'r',
#         impact[rts[3]]['AO'],theta,'g',impact[rts[4]]['AO'],theta,'c',impact[rts[5]]['AO'],theta,'m',
#         linewidth=4)
#ax.tick_params(labelsize=16)
#plt.ylim((330,400))
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=16,loc='upper center')
##plt.ylabel(r'target potential temperature (km)',fontsize=16)
#plt.xlabel(r'impact (hour K$^{-1}$)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#ax=plt.subplot(1,3,3)
#ax.plot(impact[rts[0]]['TibetanPlateau'],theta,'k',impact[rts[1]]['TibetanPlateau'],theta,'b',impact[rts[2]]['TibetanPlateau'],theta,'r',
#         impact[rts[3]]['TibetanPlateau'],theta,'g',impact[rts[4]]['TibetanPlateau'],theta,'c',impact[rts[5]]['TibetanPlateau'],theta,'m',
#         linewidth=4)
#ax.tick_params(labelsize=16)
#plt.ylim((330,400))
#plt.title(r'TibetanPlateauan plateau ',fontsize=16)
##plt.ylabel(r'target potential temperature (km)',fontsize=16)
#plt.xlabel(r'impact (hour K$^{-1}$)',fontsize=16)
#if figsave:
#    plt.savefig(os.path.join('figs','fig1-impact-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi)
#plt.show()

#%% Table of the cumulated impact per region and per layer range 
    
rts = ['EAD-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-global-All-h'+str(hmax),
       'EAZ-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-global-All-h'+str(hmax)]

rts = ['EAD-Box-theta-FullAMA',
       'EID-FULL-Box-theta-FullAMA',
       'EID-FULL-Box-theta-global',
       'EAZ-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-global']

# factor Delta t / Delta theta with Delta t = 6 h and Delta theta = 5 K
# units is then hour**2 degree**2 K**-1
# to change the time units to days we use an additional factor 1/24**2
# then divide per 1000 to scale properly
fc =  0.001*(6/5) * 1/576
# Converting degree to km with a  further 0.001 factor
fc *= 0.001*(6371 * np.pi / 180)**2
# integration in K over 5K steps
fc2 = fc*5
# Print cumulated sum over the levels
print()
print("cumulated impact unit km**2 day**2")
print("total cumulated impact")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL-AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'])/np.sum(contribs[rts[i]]['Monsoon'])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'])/np.sum(contribs[rts[i]]['Monsoon'])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'])/np.sum(contribs[rts[i]]['Monsoon'])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'])/np.sum(contribs[rts[i]]['Monsoon'])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'])/np.sum(contribs[rts[i]]['Monsoon'])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'])/np.sum(contribs[rts[i]]['Monsoon'])) for i in [3,4,5]])
print()
print("cumulated impact for theta < 350 K")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL-AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][0:4])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][0:4])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][0:4])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][0:4])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][0:4])/np.sum(contribs[rts[i]]['Monsoon'][0:4])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][0:4])/np.sum(contribs[rts[i]]['Monsoon'][0:4])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][0:4])/np.sum(contribs[rts[i]]['Monsoon'][0:4])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][0:4])/np.sum(contribs[rts[i]]['Monsoon'][0:4])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][0:4])/np.sum(contribs[rts[i]]['Monsoon'][0:4])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][0:4])/np.sum(contribs[rts[i]]['Monsoon'][0:4]) )for i in [3,4,5]])
print()
print("cumulated impact for 350 K < theta < 370 K")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][5:8])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][5:8])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][5:8])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][5:8])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][5:8])/np.sum(contribs[rts[i]]['Monsoon'][5:8])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][5:8])/np.sum(contribs[rts[i]]['Monsoon'][5:8])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][5:8])/np.sum(contribs[rts[i]]['Monsoon'][5:8])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][5:8])/np.sum(contribs[rts[i]]['Monsoon'][5:8])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][5:8])/np.sum(contribs[rts[i]]['Monsoon'][5:8])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][5:8])/np.sum(contribs[rts[i]]['Monsoon'][5:8])) for i in [3,4,5]])
print()
print("cumulated impact for 370 K < theta")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][9:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][9:20])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][9:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][9:20])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][9:20])/np.sum(contribs[rts[i]]['Monsoon'][9:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][9:20])/np.sum(contribs[rts[i]]['Monsoon'][9:20])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][9:20])/np.sum(contribs[rts[i]]['Monsoon'][9:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][9:20])/np.sum(contribs[rts[i]]['Monsoon'][9:20])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][9:20])/np.sum(contribs[rts[i]]['Monsoon'][9:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][9:20])/np.sum(contribs[rts[i]]['Monsoon'][9:20])) for i in [3,4,5]])

print()
print("cumulated impact for 380 K < theta")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][11:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][11:20])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][11:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][11:20])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][11:20])/np.sum(contribs[rts[i]]['Monsoon'][11:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][11:20])/np.sum(contribs[rts[i]]['Monsoon'][11:20])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][11:20])/np.sum(contribs[rts[i]]['Monsoon'][11:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][11:20])/np.sum(contribs[rts[i]]['Monsoon'][11:20])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][11:20])/np.sum(contribs[rts[i]]['Monsoon'][11:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][11:20])/np.sum(contribs[rts[i]]['Monsoon'][11:20])) for i in [3,4,5]])

print()
print("cumulated impact for 390 K < theta")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][13:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][13:20])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][13:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][13:20])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][13:20])/np.sum(contribs[rts[i]]['Monsoon'][13:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][13:20])/np.sum(contribs[rts[i]]['Monsoon'][13:20])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][13:20])/np.sum(contribs[rts[i]]['Monsoon'][13:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][13:20])/np.sum(contribs[rts[i]]['Monsoon'][13:20])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][13:20])/np.sum(contribs[rts[i]]['Monsoon'][13:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][13:20])/np.sum(contribs[rts[i]]['Monsoon'][13:20])) for i in [3,4,5]])

print()
print("cumulated impact for 400 K < theta")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][15:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][15:20])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][15:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][15:20])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][15:20])/np.sum(contribs[rts[i]]['Monsoon'][15:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][15:20])/np.sum(contribs[rts[i]]['Monsoon'][15:20])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][15:20])/np.sum(contribs[rts[i]]['Monsoon'][15:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][15:20])/np.sum(contribs[rts[i]]['Monsoon'][15:20])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][15:20])/np.sum(contribs[rts[i]]['Monsoon'][15:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][15:20])/np.sum(contribs[rts[i]]['Monsoon'][15:20])) for i in [3,4,5]])

print()
print("cumulated impact for 410 K < theta")
print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
print("FULL AMA domain                   ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][17:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][17:20])) for i in [3,4,5]])
print("Monsoon domain                    ",
      'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][17:20])) for i in [0,1,2]],\
      '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Monsoon'][17:20])) for i in [3,4,5]])
print("Land Asia except Tibetan plateau  ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][17:20])/np.sum(contribs[rts[i]]['Monsoon'][17:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AML'][17:20])/np.sum(contribs[rts[i]]['Monsoon'][17:20])) for i in [3,4,5]])
print("Seas surrounding Asia             ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][17:20])/np.sum(contribs[rts[i]]['Monsoon'][17:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['AO'][17:20])/np.sum(contribs[rts[i]]['Monsoon'][17:20])) for i in [3,4,5]])
print("Tibetan plateau                   ",
      'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][17:20])/np.sum(contribs[rts[i]]['Monsoon'][17:20])) for i in [0,1,2]],\
      '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['TibetanPlateau'][17:20])/np.sum(contribs[rts[i]]['Monsoon'][17:20])) for i in [3,4,5]])


#%%  Plots for potential temperature target (fig 1bis)

plt.figure(figsize=(11,6))
plt.suptitle('Cumulated convective impact at target isentropic levels [h'+str(hmax)+ss2+All+']',fontsize=18)
ax=plt.subplot(1,3,1)
ax.plot(fc*contribs[rts[0]]['AML'],theta,'k',fc*contribs[rts[1]]['AML'],theta,'b',fc*contribs[rts[2]]['AML'],theta,'r',
         fc*contribs[rts[3]]['AML'],theta,'g',fc*contribs[rts[4]]['AML'],theta,'c',fc*contribs[rts[5]]['AML'],theta,'m',
         linewidth=4)
plt.ylim((330,400))
#plt.xlim((0,100))
ax.tick_params(labelsize=16)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=16,loc='upper left')
plt.title(r'Land Asia no Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'impact (10$^{6}$ km$^2$ day$^2$ K$^{-1}$)',fontsize=16)
ax=plt.subplot(1,3,2)
ax.plot(fc*contribs[rts[0]]['AO'],theta,'k',fc*contribs[rts[1]]['AO'],theta,'b',fc*contribs[rts[2]]['AO'],theta,'r',
         fc*contribs[rts[3]]['AO'],theta,'g',fc*contribs[rts[4]]['AO'],theta,'c',fc*contribs[rts[5]]['AO'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
#plt.xlim((0,100))
plt.ylim((330,400))
plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
            fontsize=16,loc='upper center')
#plt.ylabel(r'target potential temperature (km)',fontsize=16)
plt.xlabel(r'impact (10$^{6}$ km$^2$ day$^2$ K$^{-1}$)',fontsize=16)
plt.title(r'Seas surrounding Asia',fontsize=16)
ax=plt.subplot(1,3,3)
ax.plot(fc*contribs[rts[0]]['TibetanPlateau'],theta,'k',fc*contribs[rts[1]]['TibetanPlateau'],theta,'b',fc*contribs[rts[2]]['TibetanPlateau'],theta,'r',
         fc*contribs[rts[3]]['TibetanPlateau'],theta,'g',fc*contribs[rts[4]]['TibetanPlateau'],theta,'c',fc*contribs[rts[5]]['TibetanPlateau'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
plt.ylim((330,400))
plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target potential temperature (km)',fontsize=16)
plt.xlabel(r'impact (10$^{6}$ km$^2$ day$^K$ K$^{-1}$)',fontsize=16)
if figsave:
    plt.savefig(os.path.join('figs','fig1bis-impact-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi,bbox_inches='tight')
plt.show()

#%%  Plots for potential temperature target (fig 1ter)
# same as fig 1bis but with log axis for the impact

plt.figure(figsize=(11,6))
plt.suptitle('Cumulated convective impact at target isentropic levels [h'+str(hmax)+ss2+All+']',fontsize=18)
ax=plt.subplot(1,3,1)
ax.semilogx(fc*contribs[rts[0]]['AML'],theta,'k',fc*contribs[rts[1]]['AML'],theta,'b',fc*contribs[rts[2]]['AML'],theta,'r',
         fc*contribs[rts[3]]['AML'],theta,'g',fc*contribs[rts[4]]['AML'],theta,'c',fc*contribs[rts[5]]['AML'],theta,'m',
         linewidth=4)
plt.ylim((330,400))
plt.xlim((1,600))
#plt.xlim((0,100))
ax.tick_params(labelsize=16)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=16,loc='upper left')
plt.title(r'Land Asia no Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'impact (10$^{6}$ deg$^2$ day$^2$ K$^{-1}$)',fontsize=16)
ax=plt.subplot(1,3,2)
ax.semilogx(fc*contribs[rts[0]]['AO'],theta,'k',fc*contribs[rts[1]]['AO'],theta,'b',fc*contribs[rts[2]]['AO'],theta,'r',
         fc*contribs[rts[3]]['AO'],theta,'g',fc*contribs[rts[4]]['AO'],theta,'c',fc*contribs[rts[5]]['AO'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
#plt.xlim((0,100))
plt.ylim((330,400))
plt.xlim((1,600))
plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
            fontsize=16,loc='upper center')
#plt.ylabel(r'target potential temperature (km)',fontsize=16)
plt.xlabel(r'impact (10$^{6}$ deg$^2$ day$^2$ K$^{-1}$)',fontsize=16)
plt.title(r'Seas surrounding Asia',fontsize=16)
ax=plt.subplot(1,3,3)
ax.semilogx(fc*contribs[rts[0]]['TibetanPlateau'],theta,'k',fc*contribs[rts[1]]['TibetanPlateau'],theta,'b',fc*contribs[rts[2]]['TibetanPlateau'],theta,'r',
         fc*contribs[rts[3]]['TibetanPlateau'],theta,'g',fc*contribs[rts[4]]['TibetanPlateau'],theta,'c',fc*contribs[rts[5]]['TibetanPlateau'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
plt.ylim((330,400))
plt.xlim((1,00))
plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target potential temperature (km)',fontsize=16)
plt.xlabel(r'impact (10$^{6}$ deg$^2$ day$^K$ K$^{-1}$)',fontsize=16)
if figsave:
    plt.savefig(os.path.join('figs','fig1ter-impact-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi,bbox_inches='tight')
plt.show()

#%% special plot for 'all' in FullAMA for EAD
plt.figure(figsize = (2,5))
plt.semilogx(fc*contribs[rts[0]]['All'],theta,'k')
plt.ylim((370,420))
plt.xlim((1,100))
ax.tick_params(labelsize=18)
plt.show()
# Find the slope
print(np.polyfit(theta[9:],np.log(fc*contribs[rts[0]]['All'][9:]),1))


#%% Plot the age (fig 2)

rts = ['EAD-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-global-All-h'+str(hmax),
       'EAZ-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-global-All-h'+str(hmax)]

rts = ['EAD-Box-theta-FullAMA',
       'EID-FULL-Box-theta-FullAMA',
       'EID-FULL-Box-theta-global',
       'EAZ-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-global']

  
plt.figure(figsize=(11,6))
plt.suptitle('Mean age with respect to convection [h'+str(hmax)+ss2+All+']',fontsize=18)
ax=plt.subplot(1,3,1)
ax.plot(age[rts[0]]['AML'],theta,'k',age[rts[1]]['AML'],theta,'b',age[rts[2]]['AML'],theta,'r',
         age[rts[3]]['AML'],theta,'g',age[rts[4]]['AML'],theta,'c',age[rts[5]]['AML'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'mean age (day)',fontsize=16)
plt.xlim((0,40))
plt.ylim((330,400))
ax=plt.subplot(1,3,2)
ax.plot(age[rts[0]]['AO'],theta,'k',age[rts[1]]['AO'],theta,'b',age[rts[2]]['AO'],theta,'r',
         age[rts[3]]['AO'],theta,'g',age[rts[4]]['AO'],theta,'c',age[rts[5]]['AO'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
            fontsize=16,loc='upper center')
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'mean age (day)',fontsize=16)
plt.title(r'Seas surrounding Asia',fontsize=16)
plt.xlim((0,40))
plt.ylim((330,400))
ax=plt.subplot(1,3,3)
ax.plot(age[rts[0]]['TibetanPlateau'],theta,'k',age[rts[1]]['TibetanPlateau'],theta,'b',age[rts[2]]['TibetanPlateau'],theta,'r',
         age[rts[3]]['TibetanPlateau'],theta,'g',age[rts[4]]['TibetanPlateau'],theta,'c',age[rts[5]]['TibetanPlateau'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'mean age (day)',fontsize=16)
plt.xlim((0,40))
plt.ylim((330,400))
if figsave:
    plt.savefig(os.path.join('figs','fig2-age-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi,bbox_inches='tight')
plt.show()

#%% Plot thet (fig 3)

rts = ['EAD-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-global-All-h'+str(hmax),
       'EAZ-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-global-All-h'+str(hmax)]

rts = ['EAD-Box-theta-FullAMA',
       'EID-FULL-Box-theta-FullAMA',
       'EID-FULL-Box-theta-global',
       'EAZ-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-global']


plt.figure(figsize=(11,6))
plt.suptitle('Mean potential temperature of the convective source [h'+str(hmax)+ss2+All+']',fontsize=18)
plt.subplot(1,3,1)
plt.plot(thet[rts[0]]['AML'],theta,'k',thet[rts[1]]['AML'],theta,'b',thet[rts[2]]['AML'],theta,'r',
         thet[rts[3]]['AML'],theta,'g',thet[rts[4]]['AML'],theta,'c',thet[rts[5]]['AML'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\theta$ (K)',fontsize=16)
plt.xlim((350,380))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(thet[rts[0]]['AO'],theta,'k',thet[rts[1]]['AO'],theta,'b',thet[rts[2]]['AO'],theta,'r',
         thet[rts[3]]['AO'],theta,'g',thet[rts[4]]['AO'],theta,'c',thet[rts[5]]['AO'],theta,'m',
         linewidth=4)
plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
            fontsize=16,loc='lower right')
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.tick_params(labelsize=16)
plt.xlabel(r'source $\theta$ (K)',fontsize=16)
plt.title(r'Seas surrounding Asia',fontsize=16)
plt.xlim((350,380))
plt.ylim((330,400))
plt.subplot(1,3,3)
plt.plot(thet[rts[0]]['TibetanPlateau'],theta,'k',thet[rts[1]]['TibetanPlateau'],theta,'b',thet[rts[2]]['TibetanPlateau'],theta,'r',
         thet[rts[3]]['TibetanPlateau'],theta,'g',thet[rts[4]]['TibetanPlateau'],theta,'c',thet[rts[5]]['TibetanPlateau'],theta,'m',
         linewidth=4)
plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.tick_params(labelsize=16)
plt.xlabel(r'source $\theta$ (K)',fontsize=16)
plt.xlim((350,380))
plt.ylim((330,400))
if figsave:
    plt.savefig(os.path.join('figs','fig3-theta_source-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi,bbox_inches='tight')
plt.show()

#%% Plot z (fig 4)

rts = ['EAD-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-global-All-h'+str(hmax),
       'EAZ-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-global-All-h'+str(hmax)]

rts = ['EAD-Box-theta-FullAMA',
       'EID-FULL-Box-theta-FullAMA',
       'EID-FULL-Box-theta-global',
       'EAZ-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-global']

    
plt.figure(figsize=(11,6))
plt.suptitle('Mean barometric altitude of the convective source [h'+str(hmax)+ss2+All+']',fontsize=18)
plt.subplot(1,3,1)
plt.plot(z[rts[0]]['AML'],theta,'k',z[rts[1]]['AML'],theta,'b',z[rts[2]]['AML'],theta,'r',
         z[rts[3]]['AML'],theta,'g',z[rts[4]]['AML'],theta,'c',z[rts[5]]['AML'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl;',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source baro. alt. (km)',fontsize=16)
plt.xlim((11,17))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(z[rts[0]]['AO'],theta,'k',z[rts[1]]['AO'],theta,'b',z[rts[2]]['AO'],theta,'r',
         z[rts[3]]['AO'],theta,'g',z[rts[4]]['AO'],theta,'c',z[rts[5]]['AO'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
            fontsize=16,loc='lower right')
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source baro. alt. (km)',fontsize=16)
plt.title(r'Seas surrounding Asia',fontsize=16)
plt.xlim((11,17))
plt.ylim((330,400))
plt.subplot(1,3,3)
plt.plot(z[rts[0]]['TibetanPlateau'],theta,'k',z[rts[1]]['TibetanPlateau'],theta,'b',z[rts[2]]['TibetanPlateau'],theta,'r',
         z[rts[3]]['TibetanPlateau'],theta,'g',z[rts[4]]['TibetanPlateau'],theta,'c',z[rts[5]]['TibetanPlateau'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source baro. alt. (km)',fontsize=16)
plt.xlim((11,17))
plt.ylim((330,400))
if figsave:
    plt.savefig(os.path.join('figs','fig4-z_source-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi,bbox_inches='tight')
plt.show()

#%% Plot dthet (fig 5)

rts = ['EAD-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-global-All-h'+str(hmax),
       'EAZ-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-global-All-h'+str(hmax)]

rts = ['EAD-Box-theta-FullAMA',
       'EID-FULL-Box-theta-FullAMA',
       'EID-FULL-Box-theta-global',
       'EAZ-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-global']


plt.figure(figsize=(11,6))
plt.suptitle('Mean potential temperature displacement from the convective source [h'+str(hmax)+ss2+All+']',fontsize=18)
plt.subplot(1,3,1)
plt.plot(dthet[rts[0]]['AML'],theta,'k',dthet[rts[1]]['AML'],theta,'b',dthet[rts[2]]['AML'],theta,'r',
         dthet[rts[3]]['AML'],theta,'g',dthet[rts[4]]['AML'],theta,'c',dthet[rts[5]]['AML'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
plt.xlim((-30,50))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(dthet[rts[0]]['AO'],theta,'k',dthet[rts[1]]['AO'],theta,'b',dthet[rts[2]]['AO'],theta,'r',
         dthet[rts[3]]['AO'],theta,'g',dthet[rts[4]]['AO'],theta,'c',dthet[rts[5]]['AO'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
            fontsize=16,loc='lower right')
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
plt.title(r'Seas surrounding Asia',fontsize=16)
plt.xlim((-30,50))
plt.ylim((330,400))
plt.subplot(1,3,3)
plt.plot(dthet[rts[0]]['TibetanPlateau'],theta,'k',dthet[rts[1]]['TibetanPlateau'],theta,'b',dthet[rts[2]]['TibetanPlateau'],theta,'r',
         dthet[rts[3]]['TibetanPlateau'],theta,'g',dthet[rts[4]]['TibetanPlateau'],theta,'c',dthet[rts[5]]['TibetanPlateau'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
plt.xlim((-30,50))
plt.ylim((330,400))
if figsave:
    plt.savefig(os.path.join('figs','fig5-dtheta-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi,bbox_inches='tight')
plt.show()

#%% Plot dz (fig 6)
    
rts = ['EAD-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EID-FULL-Box-theta-global-All-h'+str(hmax),
       'EAZ-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-FullAMA-All-h'+str(hmax),
       'EIZ-FULL-Box-theta-global-All-h'+str(hmax)]

rts = ['EAD-Box-theta-FullAMA',
       'EID-FULL-Box-theta-FullAMA',
       'EID-FULL-Box-theta-global',
       'EAZ-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-FullAMA',
       'EIZ-FULL-Box-theta-global']

plt.figure(figsize=(11,6))
plt.suptitle('Mean barometric altitude displacement from the convective source [h'+str(hmax)+ss2+All+']',fontsize=18)
plt.subplot(1,3,1)
plt.plot(dz[rts[0]]['AML'],theta,'k',dz[rts[1]]['AML'],theta,'b',dz[rts[2]]['AML'],theta,'r',
         dz[rts[3]]['AML'],theta,'g',dz[rts[4]]['AML'],theta,'c',dz[rts[5]]['AML'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
plt.xlim((-4,4))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(dz[rts[0]]['AO'],theta,'k',dz[rts[1]]['AO'],theta,'b',dz[rts[2]]['AO'],theta,'r',
         dz[rts[3]]['AO'],theta,'g',dz[rts[4]]['AO'],theta,'c',dz[rts[5]]['AO'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
            fontsize=16,loc='lower right')
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
plt.title(r'Seas surrounding Asia',fontsize=16)
plt.xlim((-4,4))
plt.ylim((330,400))
plt.subplot(1,3,3)
plt.plot(dz[rts[0]]['TibetanPlateau'],theta,'k',dz[rts[1]]['TibetanPlateau'],theta,'b',dz[rts[2]]['TibetanPlateau'],theta,'r',
         dz[rts[3]]['TibetanPlateau'],theta,'g',dz[rts[4]]['TibetanPlateau'],theta,'c',dz[rts[5]]['TibetanPlateau'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
plt.xlim((-4,4))
plt.ylim((330,400))
if figsave:
    plt.savefig(os.path.join('figs','fig6-dz-theta-h'+str(hmax)+All+ss3+'.png'),dpi=dpi,bbox_inches='tight')
plt.show()

#%% Plots for barometric altitude target
##%% Plot the contribution
#   
#rts = ['EAD-N-baro-source','EID-FULL-N-baro-sourceB','EID-FULL-N-baro-target',
#       'EAZ-N-baro-source','EIZ-FULL-N-baro-sourceB','EIZ-FULL-N-baro-target']
#plt.figure(figsize=(12,8))
#ax=plt.subplot(1,3,1)
#ax.plot(contribs[rts[0]]['AML'],baro,'k',contribs[rts[1]]['AML'],baro,'b',contribs[rts[2]]['AML'],baro,'r',
#         contribs[rts[3]]['AML'],baro,'g',contribs[rts[4]]['AML'],baro,'c',contribs[rts[5]]['AML'],baro,'m',
#         linewidth=4)
#ax.tick_params(labelsize=16)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'convective impact ($10^9$km h)',fontsize=16)
#ax=plt.subplot(1,3,2)
#ax.plot(contribs[rts[0]]['AO'],baro,'k',contribs[rts[1]]['AO'],baro,'b',contribs[rts[2]]['AO'],baro,'r',
#         contribs[rts[3]]['AO'],baro,'g',contribs[rts[4]]['AO'],baro,'c',contribs[rts[5]]['AO'],baro,'m',
#         linewidth=4)
#ax.tick_params(labelsize=16)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],fontsize=14)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'convective impact ($10^9$km h)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#ax=plt.subplot(1,3,3)
#ax.plot(contribs[rts[0]]['TibetanPlateau'],baro,'k',contribs[rts[1]]['TibetanPlateau'],baro,'b',contribs[rts[2]]['TibetanPlateau'],baro,'r',
#         contribs[rts[3]]['TibetanPlateau'],baro,'g',contribs[rts[4]]['TibetanPlateau'],baro,'c',contribs[rts[5]]['TibetanPlateau'],baro,'m',
#         linewidth=4)
#ax.tick_params(labelsize=16)
#plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'convective impact ($10^9$km h)',fontsize=16)
#plt.show()
#    
##%% Plot the age
#   
#rts = ['EAD-N-baro-source','EID-FULL-N-baro-sourceB','EID-FULL-N-baro-target',
#       'EAZ-N-baro-source','EIZ-FULL-N-baro-sourceB','EIZ-FULL-N-baro-target']
#plt.figure(figsize=(12,8))
#plt.subplot(1,3,1)
#plt.plot(age[rts[0]]['AML'],baro,'k',age[rts[1]]['AML'],baro,'b',age[rts[2]]['AML'],baro,'r',
#         age[rts[3]]['AML'],baro,'g',age[rts[4]]['AML'],baro,'c',age[rts[5]]['AML'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'mean age (day)',fontsize=16)
#plt.xlim((0,30))
#plt.subplot(1,3,2)
#plt.plot(age[rts[0]]['AO'],baro,'k',age[rts[1]]['AO'],baro,'b',age[rts[2]]['AO'],baro,'r',
#         age[rts[3]]['AO'],baro,'g',age[rts[4]]['AO'],baro,'c',age[rts[5]]['AO'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'mean age (day)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((0,30))
#plt.subplot(1,3,3)
#plt.plot(age[rts[0]]['TibetanPlateau'],baro,'k',age[rts[1]]['TibetanPlateau'],baro,'b',age[rts[2]]['TibetanPlateau'],baro,'r',
#         age[rts[3]]['TibetanPlateau'],baro,'g',age[rts[4]]['TibetanPlateau'],baro,'c',age[rts[5]]['TibetanPlateau'],baro,'m',
#         linewidth=4)
#plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'mean age (day)',fontsize=16)
#plt.xlim((0,30))
#plt.show()
#
##%% Plot thet 
#    
#rts = ['EAD-N-baro-source','EID-FULL-N-baro-sourceB','EID-FULL-N-baro-target',
#       'EAZ-N-baro-source','EIZ-FULL-N-baro-sourceB','EIZ-FULL-N-baro-target']
#plt.figure(figsize=(12,8))
#plt.subplot(1,3,1)
#plt.plot(thet[rts[0]]['AML'],baro,'k',thet[rts[1]]['AML'],baro,'b',thet[rts[2]]['AML'],baro,'r',
#         thet[rts[3]]['AML'],baro,'g',thet[rts[4]]['AML'],baro,'c',thet[rts[5]]['AML'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\theta$ (K)',fontsize=16)
#plt.xlim((350,430))
#plt.subplot(1,3,2)
#plt.plot(thet[rts[0]]['AO'],baro,'k',thet[rts[1]]['AO'],baro,'b',thet[rts[2]]['AO'],baro,'r',
#         thet[rts[3]]['AO'],baro,'g',thet[rts[4]]['AO'],baro,'c',thet[rts[5]]['AO'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\theta$ (K)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((350,430))
#plt.subplot(1,3,3)
#plt.plot(thet[rts[0]]['TibetanPlateau'],baro,'k',thet[rts[1]]['TibetanPlateau'],baro,'b',thet[rts[2]]['TibetanPlateau'],baro,'r',
#         thet[rts[3]]['TibetanPlateau'],baro,'g',thet[rts[4]]['TibetanPlateau'],baro,'c',thet[rts[5]]['TibetanPlateau'],baro,'m',
#         linewidth=4)
#plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\theta$ (K)',fontsize=16)
#plt.xlim((350,430))
#plt.show()
#
##%% Plot z
#
#    
#rts = ['EAD-N-baro-source','EID-FULL-N-baro-sourceB','EID-FULL-N-baro-target',
#       'EAZ-N-baro-source','EIZ-FULL-N-baro-sourceB','EIZ-FULL-N-baro-target']
#plt.figure(figsize=(12,8))
#plt.subplot(1,3,1)
#plt.plot(z[rts[0]]['AML'],baro,'k',z[rts[1]]['AML'],baro,'b',z[rts[2]]['AML'],baro,'r',
#         z[rts[3]]['AML'],baro,'g',z[rts[4]]['AML'],baro,'c',z[rts[5]]['AML'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source baro. alt. (km)',fontsize=16)
#plt.xlim((13,19))
#plt.subplot(1,3,2)
#plt.plot(z[rts[0]]['AO'],baro,'k',z[rts[1]]['AO'],baro,'b',z[rts[2]]['AO'],baro,'r',
#         z[rts[3]]['AO'],baro,'g',z[rts[4]]['AO'],baro,'c',z[rts[5]]['AO'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source baro. alt. (km)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((13,19))
#plt.subplot(1,3,3)
#plt.plot(z[rts[0]]['TibetanPlateau'],baro,'k',z[rts[1]]['TibetanPlateau'],baro,'b',z[rts[2]]['TibetanPlateau'],baro,'r',
#         z[rts[3]]['TibetanPlateau'],baro,'g',z[rts[4]]['TibetanPlateau'],baro,'c',z[rts[5]]['TibetanPlateau'],baro,'m',
#         linewidth=4)
#plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source baro. alt. (km)',fontsize=16)
#plt.xlim((13,19))
#plt.show()
#
##%% Plot dthet 
#    
#rts = ['EAD-N-baro-source','EID-FULL-N-baro-sourceB','EID-FULL-N-baro-target',
#       'EAZ-N-baro-source','EIZ-FULL-N-baro-sourceB','EIZ-FULL-N-baro-target']
#plt.figure(figsize=(12,8))
#plt.subplot(1,3,1)
#plt.plot(dthet[rts[0]]['AML'],baro,'k',dthet[rts[1]]['AML'],baro,'b',dthet[rts[2]]['AML'],baro,'r',
#         dthet[rts[3]]['AML'],baro,'g',dthet[rts[4]]['AML'],baro,'c',dthet[rts[5]]['AML'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
#plt.xlim((-10,70))
#plt.subplot(1,3,2)
#plt.plot(dthet[rts[0]]['AO'],baro,'k',dthet[rts[1]]['AO'],baro,'b',dthet[rts[2]]['AO'],baro,'r',
#         dthet[rts[3]]['AO'],baro,'g',dthet[rts[4]]['AO'],baro,'c',dthet[rts[5]]['AO'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((-10,70))
#plt.subplot(1,3,3)
#plt.plot(dthet[rts[0]]['TibetanPlateau'],baro,'k',dthet[rts[1]]['TibetanPlateau'],baro,'b',dthet[rts[2]]['TibetanPlateau'],baro,'r',
#         dthet[rts[3]]['TibetanPlateau'],baro,'g',dthet[rts[4]]['TibetanPlateau'],baro,'c',dthet[rts[5]]['TibetanPlateau'],baro,'m',
#         linewidth=4)
#plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
#plt.xlim((-10,70))
#plt.show()
#
##%% Plot dz
#    
#rts = ['EAD-N-baro-source','EID-FULL-N-baro-sourceB','EID-FULL-N-baro-target',
#       'EAZ-N-baro-source','EIZ-FULL-N-baro-sourceB','EIZ-FULL-N-baro-target']
#plt.figure(figsize=(12,8))
#plt.subplot(1,3,1)
#plt.plot(dz[rts[0]]['AML'],baro,'k',dz[rts[1]]['AML'],baro,'b',dz[rts[2]]['AML'],baro,'r',
#         dz[rts[3]]['AML'],baro,'g',dz[rts[4]]['AML'],baro,'c',dz[rts[5]]['AML'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
#plt.xlim((-3,4))
#plt.subplot(1,3,2)
#plt.plot(dz[rts[0]]['AO'],baro,'k',dz[rts[1]]['AO'],baro,'b',dz[rts[2]]['AO'],baro,'r',
#         dz[rts[3]]['AO'],baro,'g',dz[rts[4]]['AO'],baro,'c',dz[rts[5]]['AO'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((-3,4))
#plt.subplot(1,3,3)
#plt.plot(dz[rts[0]]['TibetanPlateau'],baro,'k',dz[rts[1]]['TibetanPlateau'],baro,'b',dz[rts[2]]['TibetanPlateau'],baro,'r',
#         dz[rts[3]]['TibetanPlateau'],baro,'g',dz[rts[4]]['TibetanPlateau'],baro,'c',dz[rts[5]]['TibetanPlateau'],baro,'m',
#         linewidth=4)
#plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
#plt.xlim((-3,4))
#plt.show()