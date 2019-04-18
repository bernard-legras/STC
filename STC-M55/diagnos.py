# -*- coding: utf-8 -*-
"""
Diagnostic of the backward M55 runs

Created on Sun Oct 15 19:21:15 2017

@author: Bernard Legras
"""

from datetime import datetime
from os.path import join
import numpy as np
import socket
import pickle,gzip

#import matplotlib
#matplotlib.use('Agg') # to avoid requiring X output
#matplotlib.use('Qt4Agg') # anaconda spyder on Windows 10
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import collections
from STCdata import STCinst

# defines the source_range (used only in the last section)
source_range=np.array([[6.,160.],[0.,50.]])
source_binx=77;source_biny=25
source_binx=154;source_biny=50
xedges=np.arange(source_range[0,0],source_range[0,1]+0.001,(source_range[0,1]-source_range[0,0])/source_binx)
yedges=np.arange(source_range[1,0],source_range[1,1]+0.001,(source_range[1,1]-source_range[1,0])/source_biny)

# color list with 20 colors
listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']
mymap=colors.ListedColormap(listcolors)

#ccreg={'India':'blue','Tibet':'lightskyblue','China':'red','Bob':'yellow',
#       'IO':'limegreen','SCSPhi':'coral','Pen':'darkviolet','MPac':'cyan',
#       'NEAsia':'hotpink','Tot':'k','AMA':'gold'}
ccreg=collections.OrderedDict([('India', 'blue'),('Tibet','lightskyblue'),
                               ('China', 'red'),('Bob', 'yellow'),('IO','limegreen'),
                               ('SCSPhi', 'coral'),('Pen','darkviolet'),('MPac','cyan'),
                               ('NEAsia', 'hotpink'),('Tot', 'k'),('AMA', 'gold')])
regs=list(ccreg.keys())[0:-2]

# excluded regions in this study
exreg = ('NAus','SAus','CSAf','SEP','ITA','CAm','SAm')

p0 = 100000.
I_DEAD = 0x200000
I_HIT = 0x400000
I_CROSSED = 0x2000000
I_DBORNE =  0x10000000

# number of parcels launched at each time
segl = 1000

savefig=True

def main():

    if socket.gethostname() == 'Graphium':
        STC_out = 'C:\\cygwin64\\home\\berna\\data\\STC\\STC-M55-OUT'
    elif 'ciclad' in socket.gethostname():
        STC_out ='tbd'
    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        STC_out = '/homedata/legras/STC/STC-M55'
    elif socket.gethostname() == 'satie':
        STC_out = '/limbo/data/STC/STC-M55-OUT'
    elif socket.gethostname() == 'zappa':
        STC_out ='tbd'
    elif socket.gethostname() == 'gort':
        STC_out = '/dkol/data/STC/STC-M55-out'
    else:
         print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')

    all_dates = [datetime(2017,7,27),datetime(2017,7,29),datetime(2017,7,31),
                 datetime(2017,8,2), datetime(2017,8,4), datetime(2017,8,6),
                 datetime(2017,8,8), datetime(2017,8,10)]

    all_advect = ['OPZ','EAZ','EAD','EIZ','EID']

    platform = 'M55'

    suffix = 'D01'

    date = datetime(2017,8,10)

    # read the regional mask
    #** mod1 17/10/2017
    with gzip.open('STCmask2-m55.pkl') as f:
        mask = pickle.load(f)

    #for date in all_dates:
    for advect in ['OPZ','EAZ','EAD']:

        #set txt
        txt1 = platform+'  '+advect+'  '+suffix+'  '+date.strftime('%d/%m/%Y  ')
        txt2 = date.strftime('%Y-%m-%d/')+platform+'-'+advect+'-'+suffix+'-'

        # read the data
        filename = join(STC_out,platform+date.strftime('-%Y%m%d-')+advect+'-'+suffix+'.pkl')

        # read tdc file
        tdc = STCinst('tdc',date)
        # read prod0 file
        with gzip.open(filename,'rb') as f:
            prod0 = pickle.load(f)
        numpart = int(len(prod0['flag_source'])/segl)
        print('numpart ',numpart)

        # isolate convective parcels and plot a density plot
        convect = (prod0['flag_source'] & I_HIT) != 0
        crossed = (prod0['flag_source'] & I_CROSSED) != 0
        prop0 = np.sum(convect)/len(prod0['flag_source'])
        txt = txt1+'convective prop {:5.1f} %'.format(100*prop0)

        # fig 1: figure showing the density plot (derived from histbas)
        #histbas(prod0['src']['x'][convect],prod0['src']['y'][convect],date,txt=txt)
        fig1, ax =plt.subplots(figsize=[15,6])
        m = Basemap(projection='cyl',llcrnrlat=source_range[1,0],
                    urcrnrlat=source_range[1,1],llcrnrlon=source_range[0,0],
                    urcrnrlon=source_range[0,1],resolution='c')
        m.drawcoastlines(color='k'); m.drawcountries(color='k')
        meridians = np.arange(5.,160.,5.); parallels = np.arange(0.,50.,5.)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
        H,xedges,yedges=np.histogram2d(prod0['src']['x'][convect],
                                       prod0['src']['y'][convect],
                                       bins=[source_binx,source_biny],range=source_range)
        # in order to align the colormap with the figure
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right',size='3%',pad=0.05)
        aa = np.ma.masked_equal(H.T,0)
        im=ax.imshow(np.ma.log10(aa),interpolation='nearest',extent=source_range.flatten(),
                     origin='lower',cmap='jet',aspect=1.)
        ax.set_title(txt+" (log10 #)")
        fig1.colorbar(im,cax=cax,orientation='vertical')
        if savefig: plt.savefig('figs/'+txt2+date.strftime('source-density-%Y%m%d.png'))
        #if savefig: plt.savefig(date.strftime('figs/source-density-%Y%m%d.pdf'))
        plt.show()

        # plot the proportion of convective parcels as a function of time
        # reshape the convective and crossed flags
        convect2 = np.reshape(convect,(numpart,segl))
        crossed2 = np.reshape(crossed,(numpart,segl))
        # proportion of convective parcels
        prop1 = np.sum(convect2,axis=1)/segl
        # proportion of parcels which have left the domain
        prop2 = np.sum(crossed2,axis=1)/segl
        # proportion of parcels which have stayed in the domain (AMA)
        prop3 = 1 - prop1 - prop2
        # get start and end times of the TRACZILLA run
        start = tdc.STCsegmt[tdc.flight][1]
        end = tdc.STCsegmt[tdc.flight][2]
        # make a time vector
        xx = np.arange(start,end+1)
        # plot prop1 and prop2 (now obsolete)
#        fig2 = plt.figure(figsize=[13,5])
#        plt.plot(xx,100*prop1,xx,100*(prop1+prop2))
#        plt.xlabel('time (UTC s)')
#        plt.ylabel('percentage')
#        plt.title(txt)
#        plt.savefig(date.strftime('prop1-%Y%m%d.png'))
#        plt.savefig(date.strftime('prop1-%Y%m%d.pdf'))
#        plt.show()

        # analysis of convection from regions
        # determine the region of origin for all parcels
        idx = np.clip(np.floor((prod0['src']['x']-mask['lonmin'])/mask['icxy']).astype(int),0,mask['nlons']-1)
        idy = np.clip(np.floor((prod0['src']['y']-mask['latmin'])/mask['icxy']).astype(int),0,mask['nlats']-1)
        region = mask['mask'][idy,idx]

        print(len(region))
        #** mod1 17/10/2017
        # make a report on the sources of convective parcels
        ppreg = {}
        for reg in mask['regcode'].keys():
            if reg not in exreg:
                selec = (region==mask['regcode'][reg]) & convect
                ppreg[reg] = 100*np.sum(selec)/(numpart*segl)
                print('{:>10}: {:6.4f}'.format(reg,ppreg[reg]))
        #
        # dictionaries for mean age, age spectrum, regional counter, and max of
        # the regional counter
        agemean = {}
        agespec = {}
        totreg = {}
        maxtotreg = {}
        for reg in regs:
            # make a mask in 2D for the convective parcels originating from the selected regions
            # trick: the mask is itself masked in order to get only convective parcels
            region2 = np.ma.array(np.reshape(region == mask['regcode'][reg],(numpart,segl)),mask=~convect2)
            # sum the masked mask to get regional count
            totreg[reg] = np.ma.sum(region2,axis=1)/segl
            maxtotreg[reg] = np.max(totreg[reg])
            # use the masked mask to calculate mean regional age
            # first build a reshaped age which is masked following region2
            age2 = np.ma.array(np.reshape(prod0['src']['age']/86400,(numpart,segl)),mask=~convect2)
            age2[region2==False] = np.ma.masked
            # then build a replicate of the track time
            times = np.ma.array(np.tile(xx,segl).reshape(segl,numpart).T,mask=age2.mask)
            # calculate mean age as the mean over the ages of parcel coming from the region
            agemean[reg] = np.ma.mean(age2,axis=1)
            # calculate the age spectrum over the ages of parcels as a 2d histogram
            # where the x dimension is given by replicated xx (much faster than calculating
            # numpart 1d histograms)
            agespec[reg],_,_ = np.histogram2d(times.compressed(),age2.compressed(),
                   bins =(len(xx),np.arange(-0.25,30.25,0.5)))

        # plot of the proportions by filling colors between cumulated curves
        fig3,(ax1,ax2) = plt.subplots(nrows=2,figsize=[13,10],sharex=True)
        nregs = len(regs)
        # building the incremented curves
        pl=np.zeros([len(regs)+2,numpart])
        pl[0,:] = totreg[regs[0]]
        for ii in range(1,nregs):
            pl[ii,:]=pl[ii-1,:]+totreg[regs[ii]]
        # add the last ones "Tot conv" and "AMA"
        pl[nregs,:]=prop1
        pl[nregs+1,:]=prop1+prop3
        i = nregs+1
        for reg in reversed(list(ccreg.keys())):
           ax1.fill_between(xx,pl[i,:]*100.,0,label=reg,color=ccreg[reg])
           i -= 1
        #ax.plot(xx,1regs00*prop1,color='k')
        #ax1.set_xlabel('time (UTC s)')
        ax1.set_ylabel('percentage')
        ax1.set_title(txt)
        #ax.set_xticklabels([])
        #ax.legend(loc='upper right')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),
            fancybox=True, shadow=True)
        ax1.set_xlim([np.min(xx),np.max(xx)])
        ax1.axis('on')

        ax2.set_ylim(-2,31)
        ax2.set_xlim([np.min(xx),np.max(xx)])
        for reg in regs :
            aa = np.ma.masked_where(totreg[reg]<0.05,agemean[reg])
            ax2.plot(xx,aa,ccreg[reg])

        # now add the cloud proxy from mas if available
        # this is shown at the bottom of the bottom panel
        try:
            mas = STCinst('mas',date)
            br = np.ma.masked_equal(mas.var[8],999999.0)
            brs = np.sign(br-1.2)
            brs[brs==0] = 1
            brs[brs<0] = np.ma.masked
            brs = -brs
            ax2.plot(mas.x,brs,'|')
        except:
            pass

        #ax2.plot(xx,agemean['India'],'b',xx,agemean['China'],'r',
        #         xx,agemean['Bob'],'y',xx,agemean['SCSPhi'],'coral',
        #         xx,agemean['Pen'],'darkviolet',xx,agemean['Tibet'],'lightskyblue')
        ax2.set_xlabel('time (UTC s)')
        ax2.set_ylabel('age (day)')
        ax2.set_title(txt1+'mean age from source with contribution > 5%')

        if savefig: plt.savefig('figs/'+txt2+date.strftime('prop-n-age-%Y%m%d.png'))
        #if savefig: plt.savefig(date.strftime('figs/prop-n-age-%Y%m%d.pdf'))
        plt.show()

        # Now show the age spectra for regions contibuting more than 20% at some tim
        nplot=0
        # need to know how much rows in the plot
        for reg in totreg.keys():
            if maxtotreg[reg] > 0.2: nplot += 1
        fig4, axes = plt.subplots(nrows=nplot,figsize=[13,nplot*4],sharex=True)
        ii = 0
        for reg in totreg.keys():
            if maxtotreg[reg] > 0.2:
                aa = np.ma.masked_where(totreg[reg]<0.05,agemean[reg])
                bb = np.ma.masked_equal(agespec[reg].T,0)
                axes[ii].set_ylim(-1,31)
                divider = make_axes_locatable(axes[ii])
                cax = divider.append_axes('right',size='3%',pad=0.05)
                pl=axes[ii].plot(xx,aa,'k',alpha=0.5,linewidth=2)
                im=axes[ii].imshow(np.ma.log10(bb),interpolation='nearest',aspect='auto',
                       origin='lower',cmap='jet',extent=[xx[0],xx[-1],0,30])
                fig4.colorbar(im,cax=cax,orientation='vertical')
                #axes[ii].colorbar(iax)
                axes[ii].axis('on')
                axes[ii].set_title(txt1+'age spectrum for region '+reg+' (log 10 #)')
                ii += 1
        if savefig: plt.savefig('figs/'+txt2+date.strftime('age-spectrum-%Y%m%d.png'))
        #if savefig: plt.savefig(date.strftime('figs/age-spectrum-%Y%m%d.pdf'))
        plt.show()

if __name__ == '__main__':
    main()
