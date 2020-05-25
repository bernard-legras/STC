# -*- coding: utf-8 -*-
"""

Analysis of the output of the convsrcSAF processing of the backward trajectories.

This code needs to be commented and a few problems to be fixed

TO DO : reactivate and update the water vapor part

Created on Wed Feb  8 23:26:52 2017

@author:  Bernard Legras
"""
import numpy as np
import pickle, gzip
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as clb
#from cartopy import feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from scipy.ndimage import gaussian_filter
import deepdish as dd
import socket
from os.path import join
import constants as cst

global h_hits

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
domain = np.array([[-10.,160.],[0.,50.],])
# bins of 1° used in the source domain
# and corresponding to the discretization in the target domain (initial time)
# for FullAMA initialization
binx = 170; biny = 50
deltay = (domain[1,1]-domain[1,0])/biny
deltax = (domain[0,1]-domain[0,0])/binx
ycent = np.arange(domain[1,0] + 0.5*deltay,domain[1,1],deltay)
xcent = np.arange(domain[0,0] + 0.5*deltax,domain[0,1],deltax)
yedge = np.arange(domain[1,0],domain[1,1]+0.1*deltay,deltay)
xedge = np.arange(domain[0,0],domain[0,1]+0.1*deltax,deltax)

# for global initialization
domain_g = np.array([[-179,181.],[-40,40.],])
binx_g = 360; biny_g = 80
deltay_g = (domain_g[1,1]-domain_g[1,0])/biny_g
deltax_g = (domain_g[0,1]-domain_g[0,0])/binx_g
ycent_g = np.arange(domain_g[1,0] + 0.5*deltay_g,domain_g[1,1],deltay_g)
xcent_g = np.arange(domain_g[0,0] + 0.5*deltax_g,domain_g[0,1],deltax_g)
yedge_g = np.arange(domain_g[1,0],domain_g[1,1]+0.1*deltay_g,deltay_g)
xedge_g = np.arange(domain_g[0,0],domain_g[0,1]+0.1*deltax_g,deltax_g)

#
# Generate the grid of points
#xg = np.tile(xcent,(biny,1))
#yg = np.tile(ycent,(binx,1)).T
# Size of the grid
#bloc_size = binx * biny
#xg = np.reshape(xg,bloc_size)
#yg = np.reshape(yg,bloc_size)

if socket.gethostname() == 'Graphium':
    BACK_DIR = 'C:\\cygwin64\\home\\berna\\data\\STC\\STC-BACK-OUT-SAF-silviahigh'
elif socket.gethostname() == 'satie':
    BACK_DIR = '/limbo/data/STC/STC-BACK-OUT-SAF-silviahigh'
elif socket.gethostname() in ['couperin','zappa','coltrane','puccini']:
     BACK_DIR = '/net/grapelli/limbo/data/STC/STC-BACK-OUT-SAF-OPAQ'
elif socket.gethostname() == 'gort':
    BACK_DIR = '/dkol/data/STC/STC-BACK-OUT-SAF-silviahigh'

def main():
    # load the big file containing ended
    theta=380
    # possible choices for supertype: EAD, EAZ, EID, EIZ, EID-FULL, EIZ-FULL
    supertype='EAD'
    sept = True
    check_high = True
    fname7 = join(BACK_DIR,'BACK-'+supertype+'-Jul-2017-'+str(theta)+'K.hdf5z')
    fname8 = join(BACK_DIR,'BACK-'+supertype+'-Aug-2017-'+str(theta)+'K.hdf5z')
    fname9 = join(BACK_DIR,'BACK-'+supertype+'-Sep-2017-'+str(theta)+'K.hdf5z')

    if 'FULL' in supertype:
        binx_t = binx_g
        biny_t = biny_g
        glob = True
        shrink = 0.48
    else:
        binx_t = binx
        biny_t = biny
        glob = False
        shrink = 0.63
    bloc_size = binx_t * biny_t

    h_hits_m = {}
    h_badhits_m = {}
    h_dborne_m = {}
    h_old_m = {}
    h_dead_m = {}
    age_m = {}
    H_m = {}
    H_m_high = {}
    Age_s_m = {}
    Theta_s_m = {}
    h_theta_m = {}
    nb_bottom_exit_m = {}
    i=0
    if sept:
        fnames = (fname7,fname8,fname9)
        months = 'Jul-Aug-Sep'
    else:
        fnames = (fname7,fname8)
        months = 'Jul-Aug'

    for fname in fnames:
        ended = dd.io.load(fname)

        # %% The new calculation has additional steps / analys.py due to the
        # introduction of the filtering and the theta distribution of the souurces.
        # This filtering is not optimal: it should be done in convscrSAF
        # as there is no second chance for parcels to hit a cloud after this
        # spurious hit. The estimate of the impact of this filter should be done
        # by turning it off.
        # statistics of the convectives sources
        hits = ended['flag_source'] & I_HIT == I_HIT
        excode = (ended['flag_source'] >> 13) & 0xFF
        nb_bottom_exit_m[i] =  np.sum((excode ==7) | (excode == 1))
        x_hits = ended['src']['x'][hits]
        y_hits = ended['src']['y'][hits]
        t_hits = ended['src']['t'][hits]
        p_hits = ended['src']['p'][hits]
        theta_hits = t_hits*(p_hits/cst.p0)**(-cst.kappa)
        del t_hits,p_hits
        age_hits = ended['src']['age'][hits]
        #filter the spurious sources at high latitudes (consistent with forward analysis)
        # adding however suppression of spurious European sources
        good_sources = ((y_hits<40) | (theta_hits<360)) & ((y_hits<36) | (x_hits>40))
        x_hits = x_hits[good_sources]
        y_hits = y_hits[good_sources]
        theta_hits = theta_hits[good_sources]
        age_hits = age_hits[good_sources]
        # make a selection of very high clouds in the remaining
        highthet = theta_hits>390

        # histograms of hits and of ages in the source domain (weighting of age by H_m done below)
        H_m[i] = np.histogram2d(y_hits,x_hits,bins=[biny,binx],range=domain[::-1])[0]
        H_m_high[i] = np.histogram2d(y_hits[highthet],x_hits[highthet],bins=[biny,binx],range=domain[::-1])[0]

        Age_s_m[i] = np.histogram2d(y_hits,x_hits,bins=[biny,binx],weights=age_hits,range=domain[::-1])[0]
        # 1D histogram of the potential temperature of the hits
        h_theta_m[i] = np.histogram(theta_hits,bins=100,range=(320,420))[0]
        # histogram of theta hits in the source domain
        Theta_s_m[i] = np.histogram2d(y_hits,x_hits,bins=[biny,binx],weights=theta_hits,range=domain[::-1])[0]

        # histograms in the target domain
        # number of parcels launched per bin in the target domain
        bin_size = int(len(ended['src']['x'])/bloc_size)
        # proportion of parcels hiting a good cloud (this is correct because the
        # histogram is made on the index that lives in the [0,numpart-1] interval)
        j_hits = (np.where(hits)[0])[good_sources]
        h_hits_m[i] = np.histogram(j_hits % bloc_size,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        # cumulated age / bin_size in the target space (weighting by h_hits done below)
        age_m[i] = np.histogram(j_hits % bloc_size,bins=bloc_size,range=(-0.5,bloc_size-0.5),
             weights=age_hits)[0]/bin_size
        del j_hits
        # proportion of deadborne parcels
        j_dborne = np.where(ended['flag_source'] & I_DBORNE == I_DBORNE)[0] % bloc_size
        h_dborne_m[i] = np.histogram(j_dborne,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        del(j_dborne)
        # proportion of parcels ending as too old
        j_old = np.where(ended['flag_source'] & I_OLD == I_OLD)[0] % bloc_size
        h_old_m[i] = np.histogram(j_old,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        del(j_old)
        # proportion of parcels ending by crossing the edges
        j_dead = np.where(ended['flag_source'] & I_CROSSED == I_CROSSED)[0] % bloc_size
        h_dead_m[i] = np.histogram(j_dead,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        del(j_dead)
        # proportion of bad sources
        j_badhits = (np.where(hits)[0])[~good_sources]
        h_badhits_m[i] = np.histogram(j_badhits % bloc_size,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
        del j_badhits
        i += 1

        if check_high:
            # Let us show what are the points with high h_theta
            # useful to understand the behaviour at theta = 400K
            high_hits = (theta_hits>390)
            print('high hits number',np.sum(high_hits))
            x_high = x_hits[high_hits]
            y_high = y_hits[high_hits]
            print('mean high',np.mean(x_high),np.mean(y_high))
            plt.plot(x_high,y_high,'+')
            plt.show()

    if sept:
        h_hits = (h_hits_m[0] + h_hits_m[1] + h_hits_m[2])/3
        h_dborne = (h_dborne_m[0] + h_dborne_m[1] + h_dborne_m[2])/3
        h_old = (h_old_m[0] + h_old_m[1] +  h_old_m[2])/3
        h_dead = (h_dead_m[0] + h_dead_m[1] + h_dead_m[2])/3
        npart = h_hits.copy()
        npart[npart==0]=1
        age = (age_m[0] + age_m[1] + age_m[2])/3/npart/86400
        nb_bottom_exit = (nb_bottom_exit_m[0]+nb_bottom_exit_m[1]+nb_bottom_exit_m[2])/3
        h_badhits = (h_badhits_m[0] + h_badhits_m[1] + h_badhits_m[2])/3
        del npart
        h_theta = h_theta_m[0]+h_theta_m[1]+h_theta_m[2]
        H = (H_m[0] + H_m[1] + H_m[2])/3
        H_high = (H_m_high[0] + H_m_high[1] + H_m_high[2])/3
        npart = H.copy()
        npart[npart==0] = 1
        Age_s = (Age_s_m[0]+Age_s_m[1]+Age_s_m[2])/3/npart/86400
        Theta_s = (Theta_s_m[0]+Theta_s_m[1]+Theta_s_m[2])/3/npart
    else:
        h_hits = 0.5*(h_hits_m[0] + h_hits_m[1])
        h_dborne = 0.5*(h_dborne_m[0] + h_dborne_m[1])
        h_old = 0.5*(h_old_m[0] + h_old_m[1])
        h_dead = 0.5*(h_dead_m[0] + h_dead_m[1])
        npart = h_hits.copy()
        npart[npart==0]=1
        age = 0.5*(age_m[0] + age_m[1])/npart/86400
        nb_bottom_exit = 0.5*(nb_bottom_exit_m[0]+nb_bottom_exit_m[1])
        h_badhits = 0.5*(h_badhits_m[0] + h_badhits_m[1])
        del npart
        h_theta = h_theta_m[0]+h_theta_m[1]
        H = 0.5*(H_m[0] + H_m[1])
        H_high = 0.5*(H_m_high[0] + H_m_high[1])
        npart = H.copy()
        npart[npart==0] = 1
        Age_s = 0.5*(Age_s_m[0]+Age_s_m[1])/npart/86400
        Theta_s = 0.5*(Theta_s_m[0]+Theta_s_m[1])/npart

    # make a save for further processing
    pickle.dump([H,Age_s,Theta_s],gzip.open('TMP.pkl','wb'))
    # print stats
    print('stats for '+supertype+' at '+str(theta)+'K')
    print('mean hit percentage               ',100*np.mean(h_hits))
    print('mean weigthed age                 ',np.mean(age*h_hits)/np.mean(h_hits))
    print('mean escape percentage            ',100*np.mean(h_dead))
    print('bottom ratio                      ',nb_bottom_exit/np.sum(h_dead)/bin_size)
    print('mean deadborne percentage         ',100*np.mean(h_dborne))
    print('mean old                          ',100*np.mean(h_old))
    print('bad hits                          ',100*np.mean(h_badhits))
    # percentage of still alive parcels
    #j_alive = np.where(ended['flag_source'] == 0)[0] % bloc_size
    #h_alive = 100*np.histogram(j_alive,bins=bloc_size,range=(-0.5,bloc_size-0.5))[0]/bin_size
    #del(j_alive)
    # %%
    # plot of the statistics
    chart(100*np.reshape(h_hits,[biny_t,binx_t]),vmin=0,vmax=100,glob=glob,shrink=shrink,
          txt="Percentage of convective hits from "+supertype+" "+str(theta)+" K "+months+" 2017",
          fgp=supertype+'-percentage-hits-'+months+'-'+str(theta)+'K')
    chart(100*np.reshape(h_old,[biny_t,binx_t]),glob=glob,shrink=shrink,
          txt="Percentage ending by age from "+supertype+" "+str(theta)+" K "+months+" 2017",
           fgp=supertype+'-percentage-age-end-'+months+'-'+str(theta)+'K')
    chart(100*np.reshape(h_badhits,[biny_t,binx_t]),vmin=0,vmax=100,glob=glob,shrink=shrink,
          txt="Percentage of bad convective hits from "+supertype+" "+str(theta)+" K "+months+" 2017",
          fgp=supertype+'-percentage-badhits-'+months+'-'+str(theta)+'K',coast_color='w')
    if supertype == 'EAD':
        vmin = {340:4,350:0,360:0,370:9,380:12,390:23,400:25}
        vmax = {340:14,350:10,360:10,370:19,380:22,390:33,400:35}
    else:
        vmin = {340:None,350:None,360:None,370:None,380:None,390:None,400:None}
        vmax = vmin
    chart(np.reshape(age,[biny_t,binx_t]),glob=glob,shrink=shrink,vmin=vmin[theta],vmax=vmax[theta],
          txt="mean age (day) from "+supertype+" "+str(theta)+" K "+months+" 2017",
          fgp=supertype+'-age-target-'+months+'-'+str(theta)+'K',
          back_field = np.reshape(h_hits,[biny_t,binx_t]), cumsum=True, truncate=True)
    if not glob:
        chart(100*np.reshape(h_dborne,[biny_t,binx_t]),glob=glob,shrink=shrink,
          txt="Percentage of deadbornes from "+supertype+" "+str(theta)+" K "+months+" 2017",
          fgp=supertype+'-percentage-deadborne-'+months+'-'+str(theta)+'K',coast_color='w')
        chart(100*np.reshape(h_dead,[biny_t,binx_t]),vmin=0,vmax=100,glob=glob,shrink=shrink,
          txt="Percentage of escape from "+supertype+" "+str(theta)+" K "+months+" 2017",
           fgp=supertype+'-percentage-escape-'+months+'-'+str(theta)+'K')


    #chart(np.reshape(h_alive,[biny,binx]),txt="percentage of still alive")
    # Rescale the distribution of sources before plotting
    #   Calculate the area and the scale factor
    area = np.sum(np.cos(np.deg2rad(ycent)))*len(xcent)
    d0 = np.sum(H)/area
    #   Rescale H with d0 and the cosine factor
    H = H / (np.cos(np.deg2rad(ycent))[:,np.newaxis]*d0)
    chart(H,txt="distribution of convective sources from "+supertype+" "+str(theta)+" K "+months+" 2017",
          fgp=supertype+'-distrib-sources-'+months+'-'+str(theta)+'K',coast_color='w',TP=True)
    H_high = H_high / (np.cos(np.deg2rad(ycent))[:,np.newaxis]*d0)
    chart(H_high,txt="distribution of high convective sources from "+supertype+" "+str(theta)+" K "+months+" 2017",
          fgp=supertype+'-distrib-high-sources-'+months+'-'+str(theta)+'K',coast_color='w',TP=True)
    if supertype == 'EAD':
        vmin = {340:5,350:0,360:0,370:7,380:7,390:20,400:25}
        vmax = {340:25,350:10,360:10,370:17,380:27,390:40,400:45}
    else:
        vmin = {340:None,350:None,360:None,370:None,380:None,390:None,400:None}
        vmax = vmin
    chart(Age_s,txt="mean age (day) from "+supertype+" "+str(theta)+" K "+months+" 2017",vmin=vmin[theta],vmax=vmax[theta],
          fgp=supertype+'-age-source-'+months+'-'+str(theta)+'K',back_field=H,cumsum=True,truncate=True,TP=True)
    chart(Theta_s,txt="mean theta (day) from "+supertype+" "+str(theta)+" K "+months+" 2017",
          fgp=supertype+'-theta-source-'+months+'-'+str(theta)+'K',back_field=H,cumsum=True,
          vmin=340,vmax=390,truncate=True,TP=True)
    # plot of the pdf of theta source
    plt.figure(figsize=(3.9,4.0))
    thetav = np.arange(320.5,420,1)
    plt.plot(h_theta/np.sum(h_theta),thetav,linewidth=6)
    plt.xlabel(r'pdf $\theta$ source (K$^{-1}$)',fontsize=22)
    plt.ylabel(r'$\theta$ (K)',fontsize=22)
    plt.ylim((340,410))
    plt.tick_params(labelsize=22)
    plt.savefig('figs/plot-source-theta-'+supertype+'-'+months+'-'+str(theta)+'K.png',dpi=300,bbox_inches='tight')
    plt.show()

    # backup of hits for further processing
    if supertype == 'EAD':
        pickle.dump([100*np.reshape(h_hits,[biny_t,binx_t]),xcent,ycent],open('hits_'+str(theta)+'.pkl','wb'))


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
def chart(field,txt="",vmin=None,vmax=None,fgp="",back_field=None,cumsum=False,
          shrink=0.63,glob=False,truncate=False,coast_color='k',TP=False):
    """ Plots a 2d array field with colormap between min and max"""
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fs = 18
    if(len(field.shape)>2):
        print("The field should be 2d")
        return -1
    try:
        n1=plt.get_fignums()[-1]+1
        fig=plt.figure(plt.get_fignums()[-1]+1,figsize=[13,6])
    except:
        n1=1
        fig=plt.figure(n1+1,figsize=[13,6])
    ax = fig.add_subplot(111)
    cm_lon =0
    # ACHTUNG !!: if we are to plot accross dateline, set cm to 180
    proj = ccrs.PlateCarree(central_longitude=cm_lon)
    fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
    ax = plt.axes(projection = proj)
    if vmin==None:
        vmin=np.min(field)
    if vmax==None:
        vmax=np.max(field)
        vmax=np.floor(vmax)+1
    print(vmin,vmax,(vmax-vmin)/mymap.N)
    bounds=np.arange(vmin,vmax*(1+0.0001),(vmax-vmin)/mymap.N)
    norm=colors.BoundaryNorm(bounds,mymap.N)
    if back_field is not None:
            if cumsum:
                # Here we plot levels for areas that contain a percentage of the sum
                h,edges = np.histogram(back_field,bins=200)
                cc = 0.5*(edges[:-1]+edges[1:])
                ee = np.cumsum((cc*h/np.sum(cc*h))[::-1])[::-1]
                nl = [np.argmin(np.abs(ee-x)) for x in [0.9,0.7,0.5,0.3,0.1]]
                nl95 = np.argmin(np.abs(ee-0.95))
                if glob:
                    CS=ax.contour(xcent_g,ycent_g,gaussian_filter(back_field,2),
                              transform=proj,levels=edges[nl],linewidths=3)
                else:
                    CS=ax.contour(xcent,ycent,gaussian_filter(back_field,2),
                              transform=proj,levels=edges[nl],linewidths=3)
                strs = ['90%','70%','50%','30%','10%']
                fmt={}
                for l,s in zip(CS.levels,strs):
                    fmt[l] = s
                #plt.clabel(CS,CS.levels[::2],inline=True,fmt=fmt,fontsize=12)
                plt.clabel(CS,inline=True,fmt=fmt,fontsize=fs)
                if truncate:
                    field = np.ma.array(field)
                    field[gaussian_filter(back_field,2) < edges[nl95]] = np.ma.masked
            else:
                if glob:
                    CS=ax.contour(xcent_g,ycent_g,gaussian_filter(back_field,2),
                              transform=proj,linewidths=3)
                else:
                    CS=ax.contour(xcent,ycent,gaussian_filter(back_field,2),
                              transform=proj,linewidths=3)
                plt.clabel(CS)
    if glob:
        iax=ax.pcolormesh(xedge_g,yedge_g,field,clim=[vmin,vmax],transform=proj,cmap=mymap,norm=norm)
    else:
        iax=ax.pcolormesh(xedge,yedge,field,clim=[vmin,vmax],transform=proj,cmap=mymap,norm=norm)
    # plot coastline
    ax.coastlines('50m',color=coast_color)
    if TP:
        tibetc = pickle.load(open(join('..','mkSTCmask','TibetContour.pkl'),'rb'))
        ax.plot(tibetc[:,0],tibetc[:,1],'lime')
    # Beautify layout
    xlocs = None
    if cm_lon == 180: xlocs = [0,30,60,90,120,150,180,-150,-120,-90,-60,-30]
    gl = ax.gridlines(draw_labels=True, xlocs=xlocs,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': fs}
    gl.ylabel_style = {'size': fs}
    #gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    # Eliminate white borders
    if glob:
        plt.xlim(xedge_g[0],xedge_g[-1])
        plt.ylim(yedge_g[0],yedge_g[-1])
    else:
        plt.xlim(xedge[0],xedge[-1])
        plt.ylim(yedge[0],yedge[-1])
    plt.title(txt,fontsize=fs)

    # plot adjusted colorbar
    cax,kw = clb.make_axes(ax,location='right',pad=0.02,shrink=shrink,fraction=0.10,aspect=12)
    cbar = fig.colorbar(iax,cax=cax,**kw)
    cbar.ax.tick_params(labelsize=fs)

    if len(fgp)>0:
        pass
        #plt.savefig('figs/chart-'+fgp+'.png',dpi=300,bbox_inches='tight')
    plt.show()

# %%
if __name__ == '__main__':
    main()