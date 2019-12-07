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

Plots (with target impact altitude as vertical axis) AsiaLand, AO, and Tibetan Plateau
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

Modified to be addapted to the new distribution of regions

@author: Bernard Legras
"""

import os
import pickle, gzip
import numpy as np
import matplotlib.pyplot as plt
import transit as tt
import socket
from group import group

with gzip.open(os.path.join('..','mkSTCmask','MaskCartopy2-STCforw.pkl'),'rb') as f: 
    mask =pickle.load(f)
regcode = mask['regcode']
nb_reg = len(regcode)
print(regcode)

#%%
# preliminary: estimation of the surface of the region in degree square units

coslat = np.squeeze(np.cos(np.radians(mask['lats'])))
surfreg = {}
for reg,k in regcode.items():
    surf = np.ones(mask['mask'].shape) * coslat[:,np.newaxis]
    surfreg[reg] = np.sum(surf[mask['mask']==k])
    #if surfreg[reg] ==0: 
    #    surfreg[reg] = 1
        
for gr in group.keys():
    surfreg[gr] = 0
    for reg in group[gr]:
        surfreg[gr] += surfreg[reg]

#%% obsolete
# scaling factor for hist_t and hist_s
# degree square * 10^-9 / layer thickness
# unit 10^9 km h for baro
# unit 10^9 km^2 h K^-1 for theta 
#ff={'baro':1/0.5,'theta':1/5}

#%%
figsave = True

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
            for gr in ['Asia','Land','Ocean','Tibet','All']:
                contribs[run_type][gr] = np.zeros(nlev)
                age[run_type][gr] = np.zeros(nlev)
                thet[run_type][gr] = np.zeros(nlev)
                z[run_type][gr] = np.zeros(nlev)
                dthet[run_type][gr] = np.zeros(nlev)
                dz[run_type][gr] = np.zeros(nlev)
                for reg in group[gr]:
                    contribs[run_type][gr] += contribs[run_type][reg]
                    age[run_type][gr] += age[run_type][reg]
                    thet[run_type][gr] += thet[run_type][reg]
                    z[run_type][gr] += z[run_type][reg]
                    dthet[run_type][gr] += dthet[run_type][reg]
                    dz[run_type][gr] += dz[run_type][reg]
       
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

baro = np.arange(10,20.5,0.5)
theta= np.arange(325,425,5)
       
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
#ax.plot(impact[rts[0]]['Land'],theta,'k',impact[rts[1]]['Land'],theta,'b',impact[rts[2]]['Land'],theta,'r',
#         impact[rts[3]]['Land'],theta,'g',impact[rts[4]]['Land'],theta,'c',impact[rts[5]]['Land'],theta,'m',
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
#ax.plot(impact[rts[0]]['Ocean'],theta,'k',impact[rts[1]]['Ocean'],theta,'b',impact[rts[2]]['Ocean'],theta,'r',
#         impact[rts[3]]['Ocean'],theta,'g',impact[rts[4]]['Ocean'],theta,'c',impact[rts[5]]['Ocean'],theta,'m',
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
#ax.plot(impact[rts[0]]['Tibet'],theta,'k',impact[rts[1]]['Tibet'],theta,'b',impact[rts[2]]['Tibet'],theta,'r',
#         impact[rts[3]]['Tibet'],theta,'g',impact[rts[4]]['Tibet'],theta,'c',impact[rts[5]]['Tibet'],theta,'m',
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

dd = {0:["total cumulated",range(0,20)],
      1:["cumulated impact for theta target < 350 K",range(0,5)],
      2:["cumulated impact for 350K <= theta target < 370 K",range(5,9)],
      3:["cumulated for 370K <= theta target",range(9,20)],
      4:["cumulated impact for 380K <= theta target",range(11,20)],
      5:["cumulated impact for 390K <= theta target",range(13,20)],
      6:["cumulated impact for 400K <= theta target",range(15,20)],
      7:["cumulated for 410K <= theta target",range(17,20)],}
print()
print("cumulated impact unit km**2 day**2")
for i in dd:
    print()
    print(dd[i][0])
    print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
    sel = dd[i][1]
    print("FULL AMA domain                   ",
          'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][sel])) for i in [0,1,2]],\
          '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['All'][sel])) for i in [3,4,5]])
    print("Monsoon domain                    ",
          'D',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Asia'][sel])) for i in [0,1,2]],\
          '  Z',*['{:6.1f}'.format(np.sum(fc2*contribs[rts[i]]['Asia'][sel])) for i in [3,4,5]])
    print("Land Asia except Tibetan plateau  ",
          'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['Land'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [0,1,2]],\
          '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['Land'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [3,4,5]])
    print("Seas surrounding Asia             ",
          'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['Ocean'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [0,1,2]],\
          '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['Ocean'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [3,4,5]])
    print("Tibetan plateau                   ",
          'D',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['Tibet'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [0,1,2]],\
          '  Z',*['{:5.1f}%'.format(np.sum(100*contribs[rts[i]]['Tibet'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [3,4,5]])


#%% Mean potential temperature of the source for impact levels 

dd = {0:["mean theta of the source",range(0,20)],
      1:["mean theta source for theta target < 350 K",range(0,5)],
      2:["mean theta source for 350K <= theta target < 370 K",range(5,9)],
      3:["mean theta source for 370K <= theta target",range(9,20)],
      4:["mean theta source for 380K <= theta target",range(11,20)],
      5:["mean theta source for 390K <= theta target",range(13,20)],
      6:["mean theta source for 400K <= theta target",range(15,20)],
      7:["mean theta source for 410K <= theta target",range(17,20)],}
for i in dd:
    print()
    print(dd[i][0])
    print('order: EAD, EID, EID-FULL, EAZ, EIZ, EIZ-FULL' )
    sel = dd[i][1]
    print("FULL-AMA domain                   ",
          'D',*['{:6.1f}'.format(np.sum(thet[rts[i]]['All'][sel]*contribs[rts[i]]['All'][sel])/np.sum(contribs[rts[i]]['All'][sel])) for i in [0,1,2]],\
          '  Z',*['{:6.1f}'.format(np.sum(thet[rts[i]]['All'][sel]*contribs[rts[i]]['All'][sel])/np.sum(contribs[rts[i]]['All'][sel])) for i in [3,4,5]])
    print("Monsoon domain                    ",
          'D',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Asia'][sel]*contribs[rts[i]]['Asia'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [0,1,2]],\
          '  Z',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Asia'][sel]*contribs[rts[i]]['Asia'][sel])/np.sum(contribs[rts[i]]['Asia'][sel])) for i in [3,4,5]])
    print("Land Asia except Tibetan plateau  ",
          'D',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Land'][sel]*contribs[rts[i]]['Land'][sel])/np.sum(contribs[rts[i]]['Land'][sel])) for i in [0,1,2]],\
          '  Z',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Land'][sel]*contribs[rts[i]]['Land'][sel])/np.sum(contribs[rts[i]]['Land'][sel])) for i in [3,4,5]])
    print("Seas surrounding Asia             ",
          'D',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Ocean'][sel]*contribs[rts[i]]['Ocean'][sel])/np.sum(contribs[rts[i]]['Ocean'][sel])) for i in [0,1,2]],\
          '  Z',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Ocean'][sel]*contribs[rts[i]]['Ocean'][sel])/np.sum(contribs[rts[i]]['Ocean'][sel])) for i in [3,4,5]])
    print("Tibetan plateau                   ",
          'D',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Tibet'][sel]*contribs[rts[i]]['Tibet'][sel])/np.sum(contribs[rts[i]]['Tibet'][sel])) for i in [0,1,2]],\
          '  Z',*['{:6.1f}'.format(np.sum(thet[rts[i]]['Tibet'][sel]*contribs[rts[i]]['Tibet'][sel])/np.sum(contribs[rts[i]]['Tibet'][sel])) for i in [3,4,5]])


#%%  Plots for potential temperature target (fig 1bis)

plt.figure(figsize=(11,6))
plt.suptitle('Cumulated convective impact at target isentropic levels [h'+str(hmax)+ss2+All+']',fontsize=18)
ax=plt.subplot(1,3,1)
ax.plot(fc*contribs[rts[0]]['Land'],theta,'k',fc*contribs[rts[1]]['Land'],theta,'b',fc*contribs[rts[2]]['Land'],theta,'r',
         fc*contribs[rts[3]]['Land'],theta,'g',fc*contribs[rts[4]]['Land'],theta,'c',fc*contribs[rts[5]]['Land'],theta,'m',
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
ax.plot(fc*contribs[rts[0]]['Ocean'],theta,'k',fc*contribs[rts[1]]['Ocean'],theta,'b',fc*contribs[rts[2]]['Ocean'],theta,'r',
         fc*contribs[rts[3]]['Ocean'],theta,'g',fc*contribs[rts[4]]['Ocean'],theta,'c',fc*contribs[rts[5]]['Ocean'],theta,'m',
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
ax.plot(fc*contribs[rts[0]]['Tibet'],theta,'k',fc*contribs[rts[1]]['Tibet'],theta,'b',fc*contribs[rts[2]]['Tibet'],theta,'r',
         fc*contribs[rts[3]]['Tibet'],theta,'g',fc*contribs[rts[4]]['Tibet'],theta,'c',fc*contribs[rts[5]]['Tibet'],theta,'m',
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
ax.semilogx(fc*contribs[rts[0]]['Land'],theta,'k',fc*contribs[rts[1]]['Land'],theta,'b',fc*contribs[rts[2]]['Land'],theta,'r',
         fc*contribs[rts[3]]['Land'],theta,'g',fc*contribs[rts[4]]['Land'],theta,'c',fc*contribs[rts[5]]['Land'],theta,'m',
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
ax.semilogx(fc*contribs[rts[0]]['Ocean'],theta,'k',fc*contribs[rts[1]]['Ocean'],theta,'b',fc*contribs[rts[2]]['Ocean'],theta,'r',
         fc*contribs[rts[3]]['Ocean'],theta,'g',fc*contribs[rts[4]]['Ocean'],theta,'c',fc*contribs[rts[5]]['Ocean'],theta,'m',
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
ax.semilogx(fc*contribs[rts[0]]['Tibet'],theta,'k',fc*contribs[rts[1]]['Tibet'],theta,'b',fc*contribs[rts[2]]['Tibet'],theta,'r',
         fc*contribs[rts[3]]['Tibet'],theta,'g',fc*contribs[rts[4]]['Tibet'],theta,'c',fc*contribs[rts[5]]['Tibet'],theta,'m',
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
ax.plot(age[rts[0]]['Land'],theta,'k',age[rts[1]]['Land'],theta,'b',age[rts[2]]['Land'],theta,'r',
         age[rts[3]]['Land'],theta,'g',age[rts[4]]['Land'],theta,'c',age[rts[5]]['Land'],theta,'m',
         linewidth=4)
ax.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'mean age (day)',fontsize=16)
plt.xlim((0,40))
plt.ylim((330,400))
ax=plt.subplot(1,3,2)
ax.plot(age[rts[0]]['Ocean'],theta,'k',age[rts[1]]['Ocean'],theta,'b',age[rts[2]]['Ocean'],theta,'r',
         age[rts[3]]['Ocean'],theta,'g',age[rts[4]]['Ocean'],theta,'c',age[rts[5]]['Ocean'],theta,'m',
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
ax.plot(age[rts[0]]['Tibet'],theta,'k',age[rts[1]]['Tibet'],theta,'b',age[rts[2]]['Tibet'],theta,'r',
         age[rts[3]]['Tibet'],theta,'g',age[rts[4]]['Tibet'],theta,'c',age[rts[5]]['Tibet'],theta,'m',
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
plt.plot(thet[rts[0]]['Land'],theta,'k',thet[rts[1]]['Land'],theta,'b',thet[rts[2]]['Land'],theta,'r',
         thet[rts[3]]['Land'],theta,'g',thet[rts[4]]['Land'],theta,'c',thet[rts[5]]['Land'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\theta$ (K)',fontsize=16)
plt.xlim((350,380))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(thet[rts[0]]['Ocean'],theta,'k',thet[rts[1]]['Ocean'],theta,'b',thet[rts[2]]['Ocean'],theta,'r',
         thet[rts[3]]['Ocean'],theta,'g',thet[rts[4]]['Ocean'],theta,'c',thet[rts[5]]['Ocean'],theta,'m',
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
plt.plot(thet[rts[0]]['Tibet'],theta,'k',thet[rts[1]]['Tibet'],theta,'b',thet[rts[2]]['Tibet'],theta,'r',
         thet[rts[3]]['Tibet'],theta,'g',thet[rts[4]]['Tibet'],theta,'c',thet[rts[5]]['Tibet'],theta,'m',
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
plt.plot(z[rts[0]]['Land'],theta,'k',z[rts[1]]['Land'],theta,'b',z[rts[2]]['Land'],theta,'r',
         z[rts[3]]['Land'],theta,'g',z[rts[4]]['Land'],theta,'c',z[rts[5]]['Land'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl;',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source baro. alt. (km)',fontsize=16)
plt.xlim((11,17))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(z[rts[0]]['Ocean'],theta,'k',z[rts[1]]['Ocean'],theta,'b',z[rts[2]]['Ocean'],theta,'r',
         z[rts[3]]['Ocean'],theta,'g',z[rts[4]]['Ocean'],theta,'c',z[rts[5]]['Ocean'],theta,'m',
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
plt.plot(z[rts[0]]['Tibet'],theta,'k',z[rts[1]]['Tibet'],theta,'b',z[rts[2]]['Tibet'],theta,'r',
         z[rts[3]]['Tibet'],theta,'g',z[rts[4]]['Tibet'],theta,'c',z[rts[5]]['Tibet'],theta,'m',
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
plt.plot(dthet[rts[0]]['Land'],theta,'k',dthet[rts[1]]['Land'],theta,'b',dthet[rts[2]]['Land'],theta,'r',
         dthet[rts[3]]['Land'],theta,'g',dthet[rts[4]]['Land'],theta,'c',dthet[rts[5]]['Land'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
plt.xlim((-30,50))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(dthet[rts[0]]['Ocean'],theta,'k',dthet[rts[1]]['Ocean'],theta,'b',dthet[rts[2]]['Ocean'],theta,'r',
         dthet[rts[3]]['Ocean'],theta,'g',dthet[rts[4]]['Ocean'],theta,'c',dthet[rts[5]]['Ocean'],theta,'m',
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
plt.plot(dthet[rts[0]]['Tibet'],theta,'k',dthet[rts[1]]['Tibet'],theta,'b',dthet[rts[2]]['Tibet'],theta,'r',
         dthet[rts[3]]['Tibet'],theta,'g',dthet[rts[4]]['Tibet'],theta,'c',dthet[rts[5]]['Tibet'],theta,'m',
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
plt.plot(dz[rts[0]]['Land'],theta,'k',dz[rts[1]]['Land'],theta,'b',dz[rts[2]]['Land'],theta,'r',
         dz[rts[3]]['Land'],theta,'g',dz[rts[4]]['Land'],theta,'c',dz[rts[5]]['Land'],theta,'m',
         linewidth=4)
plt.tick_params(labelsize=16)
plt.title(r'Land Asia excluding Tib. Pl.',fontsize=16)
plt.ylabel(r'target potential temperature (K)',fontsize=16)
plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
plt.xlim((-4,4))
plt.ylim((330,400))
plt.subplot(1,3,2)
plt.plot(dz[rts[0]]['Ocean'],theta,'k',dz[rts[1]]['Ocean'],theta,'b',dz[rts[2]]['Ocean'],theta,'r',
         dz[rts[3]]['Ocean'],theta,'g',dz[rts[4]]['Ocean'],theta,'c',dz[rts[5]]['Ocean'],theta,'m',
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
plt.plot(dz[rts[0]]['Tibet'],theta,'k',dz[rts[1]]['Tibet'],theta,'b',dz[rts[2]]['Tibet'],theta,'r',
         dz[rts[3]]['Tibet'],theta,'g',dz[rts[4]]['Tibet'],theta,'c',dz[rts[5]]['Tibet'],theta,'m',
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
#ax.plot(contribs[rts[0]]['Land'],baro,'k',contribs[rts[1]]['Land'],baro,'b',contribs[rts[2]]['Land'],baro,'r',
#         contribs[rts[3]]['Land'],baro,'g',contribs[rts[4]]['Land'],baro,'c',contribs[rts[5]]['Land'],baro,'m',
#         linewidth=4)
#ax.tick_params(labelsize=16)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'convective impact ($10^9$km h)',fontsize=16)
#ax=plt.subplot(1,3,2)
#ax.plot(contribs[rts[0]]['Ocean'],baro,'k',contribs[rts[1]]['Ocean'],baro,'b',contribs[rts[2]]['Ocean'],baro,'r',
#         contribs[rts[3]]['Ocean'],baro,'g',contribs[rts[4]]['Ocean'],baro,'c',contribs[rts[5]]['Ocean'],baro,'m',
#         linewidth=4)
#ax.tick_params(labelsize=16)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],fontsize=14)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'convective impact ($10^9$km h)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#ax=plt.subplot(1,3,3)
#ax.plot(contribs[rts[0]]['Tibet'],baro,'k',contribs[rts[1]]['Tibet'],baro,'b',contribs[rts[2]]['Tibet'],baro,'r',
#         contribs[rts[3]]['Tibet'],baro,'g',contribs[rts[4]]['Tibet'],baro,'c',contribs[rts[5]]['Tibet'],baro,'m',
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
#plt.plot(age[rts[0]]['Land'],baro,'k',age[rts[1]]['Land'],baro,'b',age[rts[2]]['Land'],baro,'r',
#         age[rts[3]]['Land'],baro,'g',age[rts[4]]['Land'],baro,'c',age[rts[5]]['Land'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'mean age (day)',fontsize=16)
#plt.xlim((0,30))
#plt.subplot(1,3,2)
#plt.plot(age[rts[0]]['Ocean'],baro,'k',age[rts[1]]['Ocean'],baro,'b',age[rts[2]]['Ocean'],baro,'r',
#         age[rts[3]]['Ocean'],baro,'g',age[rts[4]]['Ocean'],baro,'c',age[rts[5]]['Ocean'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'mean age (day)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((0,30))
#plt.subplot(1,3,3)
#plt.plot(age[rts[0]]['Tibet'],baro,'k',age[rts[1]]['Tibet'],baro,'b',age[rts[2]]['Tibet'],baro,'r',
#         age[rts[3]]['Tibet'],baro,'g',age[rts[4]]['Tibet'],baro,'c',age[rts[5]]['Tibet'],baro,'m',
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
#plt.plot(thet[rts[0]]['Land'],baro,'k',thet[rts[1]]['Land'],baro,'b',thet[rts[2]]['Land'],baro,'r',
#         thet[rts[3]]['Land'],baro,'g',thet[rts[4]]['Land'],baro,'c',thet[rts[5]]['Land'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\theta$ (K)',fontsize=16)
#plt.xlim((350,430))
#plt.subplot(1,3,2)
#plt.plot(thet[rts[0]]['Ocean'],baro,'k',thet[rts[1]]['Ocean'],baro,'b',thet[rts[2]]['Ocean'],baro,'r',
#         thet[rts[3]]['Ocean'],baro,'g',thet[rts[4]]['Ocean'],baro,'c',thet[rts[5]]['Ocean'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\theta$ (K)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((350,430))
#plt.subplot(1,3,3)
#plt.plot(thet[rts[0]]['Tibet'],baro,'k',thet[rts[1]]['Tibet'],baro,'b',thet[rts[2]]['Tibet'],baro,'r',
#         thet[rts[3]]['Tibet'],baro,'g',thet[rts[4]]['Tibet'],baro,'c',thet[rts[5]]['Tibet'],baro,'m',
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
#plt.plot(z[rts[0]]['Land'],baro,'k',z[rts[1]]['Land'],baro,'b',z[rts[2]]['Land'],baro,'r',
#         z[rts[3]]['Land'],baro,'g',z[rts[4]]['Land'],baro,'c',z[rts[5]]['Land'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source baro. alt. (km)',fontsize=16)
#plt.xlim((13,19))
#plt.subplot(1,3,2)
#plt.plot(z[rts[0]]['Ocean'],baro,'k',z[rts[1]]['Ocean'],baro,'b',z[rts[2]]['Ocean'],baro,'r',
#         z[rts[3]]['Ocean'],baro,'g',z[rts[4]]['Ocean'],baro,'c',z[rts[5]]['Ocean'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source baro. alt. (km)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((13,19))
#plt.subplot(1,3,3)
#plt.plot(z[rts[0]]['Tibet'],baro,'k',z[rts[1]]['Tibet'],baro,'b',z[rts[2]]['Tibet'],baro,'r',
#         z[rts[3]]['Tibet'],baro,'g',z[rts[4]]['Tibet'],baro,'c',z[rts[5]]['Tibet'],baro,'m',
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
#plt.plot(dthet[rts[0]]['Land'],baro,'k',dthet[rts[1]]['Land'],baro,'b',dthet[rts[2]]['Land'],baro,'r',
#         dthet[rts[3]]['Land'],baro,'g',dthet[rts[4]]['Land'],baro,'c',dthet[rts[5]]['Land'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
#plt.xlim((-10,70))
#plt.subplot(1,3,2)
#plt.plot(dthet[rts[0]]['Ocean'],baro,'k',dthet[rts[1]]['Ocean'],baro,'b',dthet[rts[2]]['Ocean'],baro,'r',
#         dthet[rts[3]]['Ocean'],baro,'g',dthet[rts[4]]['Ocean'],baro,'c',dthet[rts[5]]['Ocean'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta \theta$ (K)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((-10,70))
#plt.subplot(1,3,3)
#plt.plot(dthet[rts[0]]['Tibet'],baro,'k',dthet[rts[1]]['Tibet'],baro,'b',dthet[rts[2]]['Tibet'],baro,'r',
#         dthet[rts[3]]['Tibet'],baro,'g',dthet[rts[4]]['Tibet'],baro,'c',dthet[rts[5]]['Tibet'],baro,'m',
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
#plt.plot(dz[rts[0]]['Land'],baro,'k',dz[rts[1]]['Land'],baro,'b',dz[rts[2]]['Land'],baro,'r',
#         dz[rts[3]]['Land'],baro,'g',dz[rts[4]]['Land'],baro,'c',dz[rts[5]]['Land'],baro,'m',
#         linewidth=4)
#plt.title(r'Land Asia excluding TibetanPlateau',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
#plt.xlim((-3,4))
#plt.subplot(1,3,2)
#plt.plot(dz[rts[0]]['Ocean'],baro,'k',dz[rts[1]]['Ocean'],baro,'b',dz[rts[2]]['Ocean'],baro,'r',
#         dz[rts[3]]['Ocean'],baro,'g',dz[rts[4]]['Ocean'],baro,'c',dz[rts[5]]['Ocean'],baro,'m',
#         linewidth=4)
#plt.legend([r'\textbf{ERA5 diabatic AMA}',r'ERA-I diabatic AMA',r'\textit{ERA-I diabatic Global}',
#            r'\textbf{ERA5 kinematic AMA}',r'ERA-I kinematic AMA',r'\textit{ERA-I kinematic Global}'],
#            fontsize=14,loc='lower right')
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
#plt.title(r'Seas surrounding Asia',fontsize=16)
#plt.xlim((-3,4))
#plt.subplot(1,3,3)
#plt.plot(dz[rts[0]]['Tibet'],baro,'k',dz[rts[1]]['Tibet'],baro,'b',dz[rts[2]]['Tibet'],baro,'r',
#         dz[rts[3]]['Tibet'],baro,'g',dz[rts[4]]['Tibet'],baro,'c',dz[rts[5]]['Tibet'],baro,'m',
#         linewidth=4)
#plt.title(r'Tibetan plateau ',fontsize=16)
#plt.ylabel(r'target barometric altitude (km)',fontsize=16)
#plt.xlabel(r'source $\Delta$ baro. alt. (km)',fontsize=16)
#plt.xlim((-3,4))
#plt.show()