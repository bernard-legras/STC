# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 23:26:52 2017

@author:  Bernard Legras
"""
import numpy as np
import pickle, gzip
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
import deepdish as dd
import socket
from os.path import join
# flags
I_DEAD = 0x200000
I_HIT = 0x400000
I_OLD = 0x800000
I_CROSSED = 0x2000000
I_DBORNE =  0x1000000

# color list with 20 colors

listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']
mymap=colors.ListedColormap(listcolors)

# %%
""" Defines the domain
"""
domain=np.array([[-179.,181.],[-30.,50.]])
binx = 360; biny = 80

deltay = (domain[1,1]-domain[1,0])/biny
deltax = (domain[0,1]-domain[0,0])/binx
ycent = np.arange(domain[1,0] + 0.5*deltay,domain[1,1],deltay)
xcent = np.arange(domain[0,0] + 0.5*deltax,domain[0,1],deltax)
yedge = np.arange(domain[1,0],domain[1,1]+0.1*deltay,deltay)
xedge = np.arange(domain[0,0],domain[0,1]+0.1*deltax,deltax)
# Generate the grid of points
xg = np.tile(xcent,(biny,1))
yg = np.tile(ycent,(binx,1)).T
# Size of the grid
bloc_size = binx * biny
xg = np.reshape(xg,bloc_size)
yg = np.reshape(yg,bloc_size)

if socket.gethostname() == 'Graphium':
    BACK_DIR = 'C:\\cygwin64\\home\\berna\\data\\STC\\STC-BACK-OUT-SAF-OPAQ'
elif socket.gethostname() == 'grapelli':
    BACK_DIR = '/limbo/data/STC/STC-BACK-OUT-SAF-OPAQ'
elif socket.gethostname() in ['couperin','zappa','coltrane','puccini']:
     BACK_DIR = '/net/grapelli/limbo/data/STC/STC-BACK-OUT-SAF-OPAQ'
elif socket.gethostname() == 'gort':
    BACK_DIR = '/dkol/data/STC/STC-BACK-OUT-SAF-OPAQ'

def main():
        # load the big file containing ended       
    theta=400
    fname1 = join(BACK_DIR,'BACK-EIZ-FULL-Jul-2017-'+str(theta)+'K.hdf5z')
    fname2 = join(BACK_DIR,'BACK-EIZ-FULL-Aug-2017-'+str(theta)+'K.hdf5z')

    h_hits_m = {}
    h_dborne_m = {}
    h_old_m = {}
    h_dead_m = {}
    H_m = {}
    
    i=0
    for fname in [fname1,fname2]:   
        ended = dd.io.load(fname)
    
        # number of parcels launched per bin
        bin_size = int(len(ended['src']['x'])/bloc_size)
        # percentage of parcels hiting a cloud
        j_hits = np.where(ended['flag_source'] & I_HIT == I_HIT)[0] % bloc_size
        h_hits_m[i] = 100*np.histogram(j_hits,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        del(j_hits)
        # percentage of deadborne parcels
    #    j_dborne = np.where(ended['flag_source'] & I_DBORNE == I_DBORNE)[0] % bloc_size
    #    h_dborne = 100*np.histogram(j_dborne,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
    #    del(j_dborne)
        # percentage of parcels ending as too old
        j_old = np.where(ended['flag_source'] & I_OLD == I_OLD)[0] % bloc_size
        h_old_m[i] = 100*np.histogram(j_old,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        del(j_old)
        # percentage of parcels ending by crossing the edges
        j_dead = np.where(ended['flag_source'] & I_CROSSED == I_CROSSED)[0] % bloc_size
        h_dead_m[i] = 100*np.histogram(j_dead,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        del(j_dead)
        x_hits = ended['src']['x'][ended['flag_source'] & I_HIT == I_HIT]
        y_hits = ended['src']['y'][ended['flag_source'] & I_HIT == I_HIT]
        H_m[i]=np.histogram2d(y_hits,x_hits,bins=[biny,binx],range=domain[::-1])[0]
        i += 1
    
    h_hits = 0.5*(h_hits_m[0] + h_hits_m[1])
    #h_dborne = 0.5*(h_dborne_m[0] + h_dborne_m[1])
    h_old = 0.5*(h_old_m[0] + h_old_m[1])
    h_dead = 0.5*(h_dead_m[0] + h_dead_m[1])
    # percentage of still alive parcels
    #j_alive = np.where(ended['flag_source'] == 0)[0] % bloc_size
    #h_alive = 100*np.histogram(j_alive,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
    #del(j_alive)
    # %%
    # plot of the statistics
    chart(np.reshape(h_hits,[biny,binx]),vmin=0,vmax=100,txt="Percentage of convective hits from EIZ FULL "+str(theta)+" K Jul-Aug 2017",
          fgp='EIZ-FULL-percentage-hits-'+str(theta)+'K')
    #chart(np.reshape(h_dborne,[biny,binx]),txt="percentage of deadborne")
    chart(np.reshape(h_old,[biny,binx]),vmin=0,vmax=100,txt="percentage of ending by age")
    chart(np.reshape(h_dead,[biny,binx]),vmin=0,vmax=100,txt="percentage of escape")
    #chart(np.reshape(h_alive,[biny,binx]),txt="percentage of still alive")

    # %%
    # statistics of the convectives sources
    H = 0.5*(H_m[0] + H_m[1])
    chart(H,txt="distribution of convective sources from EIZ FULL "+str(theta)+" K Jul-Aug 2017",fgp='EIZ-FULL-distrib-sources-'+str(theta)+'K')


#    # %% Moisture
#    # load MLS profile data
#    prof_file = '../MLS/MeanMLSProf-H2O-2016-07.pkl'
#    prof = pickle.load(gzip.open(prof_file,'r'))
#    # CORRECTIVE STEPS
#    # copy and correction of water vapour (as satratio was in kg/kg)
#    rvs=29/18*ended['rvs']
#    """ take care of the remaining parcels (temporary fix before removal from
#     convsrc3 processing)
#    About 3500 parcels are in this case with no info kept in ended.
#    Put them all at a single location.
#    This is not going to have a role as it represents 0,01%"""
#    jpnull=np.where(ended['p']==0)[0]
#    ended['p'][jpnull] = 10000.
#    ended['y'][jpnull] = 25.
#    # set rvs for the DEADBORNE parcels
#    j_db = np.where(ended['flag'] & I_DBORNE == I_DBORNE)[0]
#    y_end = ended['y'][j_db]
#    logp_end = np.log(ended['p'][j_db])
#    idy = np.digitize(y_end,prof['LatEdges'])-1
#    idp = np.digitize(logp_end,prof['LogPressureEdges'])-1
#    rvs[j_db] =np.minimum(rvs[j_db],prof['H2O'][idp,idy])
#    # 1) Plot of the raw distribution of moisture
#    num_blocs = int(len(ended['x'])/bloc_size)
#    rvs_r = np.mean(np.reshape(rvs,[num_blocs,bloc_size]),0)
#    chart(np.reshape(1.e6*rvs_r,[biny,binx]),txt='raw moisture')
#    # 2) contribution of hits to water vapour
#    j_hits = np.where(ended['flag'] & I_HIT)[0] % bloc_size
#    hh = np.maximum(np.histogram(j_hits,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0],1)
#    rvs_hits = rvs[ended['flag'] & I_HIT == I_HIT]
#    rvs_c = np.zeros(bloc_size)
#    for j in range(len(rvs_hits)):
#        rvs_c[j_hits[j]] += rvs_hits[j]
#    rvs_c = rvs_c/hh
#    chart(np.reshape(1.e6*rvs_c,[biny,binx]),txt="hit moisture")

#    #%%
#    # 3) correction of non hits rvs with MLS values on the edge
#    prof_file = '../MLS/MeanMLSProf-H2O-2016-07.pkl'
#    prof = pickle.load(gzip.open(prof_file,'r'))
#    j_nohits = np.where(ended['flag'] & I_HIT ==0)[0]
#    y_end = ended['y'][j_nohits]
#    logp_end = np.log(ended['p'][j_nohits])
#    idy = np.digitize(y_end,prof['LatEdges'])-1
#    idp = np.digitize(logp_end,prof['LogPressureEdges'])-1
#
#    rvs[j_nohits] =np.minimum(rvs[j_nohits],prof['H2O'][idp,idy])
#    rvs_t = np.mean(np.reshape(rvs,[num_blocs,bloc_size]),0)
#    chart(np.reshape(1.e6*rvs_t,[biny,binx]),txt='moisture v1')
#
#    # 4) now set the rvs on hit at the mls value same alt and lat
#    j_hits=np.where(ended['flag'] & I_HIT == I_HIT)[0]
#    y_end = ended['y'][j_hits]
#    logp_end = np.log(ended['p'][j_hits])
#    idy = np.digitize(y_end,prof['LatEdges'])-1
#    idp = np.digitize(logp_end,prof['LogPressureEdges'])-1
#
#    rvs[j_hits] =np.minimum(rvs[j_hits],prof['H2O'][idp,idy])
#    rvs_tt = np.sum(np.reshape(rvs,[num_blocs,bloc_size]),0)/num_blocs
#    chart(np.reshape(1.e6*rvs_tt,[biny,binx]),txt='moisture v2')

# %%
def chart(field,txt="",vmin=0,vmax=0,fgp=""):
    """ Plots a 2d array field with colormap between min and max"""
    if(len(field.shape)>2):
        print("The field should be 2d")
        return -1
    try:
        n1=plt.get_fignums()[-1]+1
        fig=plt.figure(plt.get_fignums()[-1]+1,figsize=[13,6])
    except:
        n1=1
        fig=plt.figure(n1+1,figsize=[13,6])
    m = Basemap(projection='cyl',llcrnrlat=domain[1,0],urcrnrlat=domain[1,1],
        llcrnrlon=domain[0,0],urcrnrlon=domain[0,1],resolution='c')
    m.drawcoastlines(color='w'); m.drawcountries(color='k')
    meridians = np.arange(-180.,180.,30.); parallels = np.arange(-50.,50.,20.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
    if vmin==0:
        vmin=np.min(field)
    if vmax==0:
        vmax=np.max(field)
    bounds=np.arange(vmin,vmax*(1+0.0001),(vmax-vmin)/mymap.N)
    norm=colors.BoundaryNorm(bounds,mymap.N)
    iax=plt.imshow(field,interpolation='nearest',extent=domain.flatten(),
                   clim=[vmin,vmax],origin='lower',cmap=mymap,norm=norm,aspect=1.)
    plt.title(txt,fontsize=18)
    #plt.xlabel('longitude')
    #plt.ylabel('latitude')
    cax = fig.add_axes([0.91, 0.26, 0.03, 0.5])
    cbar=fig.colorbar(iax,cax=cax)
    cbar.ax.tick_params(labelsize=18)
    #if len(fgp)>0:
        #plt.savefig('figs/chart-'+fgp+'.png')   
    plt.show()

# %%
if __name__ == '__main__':
    main()