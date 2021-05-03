#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preparation of shifted track runs SK-2 and SK-3 for Sergey

To be improved if used more extensively

Created on Sun Jun 10 03:59:50 2018

@author: Bernard Legras
"""

from datetime import datetime, timedelta
import numpy as np
import math
from scipy.interpolate import RectBivariateSpline, interp1d
from STCdata import STCinst
from ECMWF_N import ECMWF
import constants as cst
import io107

def interp3d_thet(data,theta):
    """ Produces an interpolation of the temperature field to a 2d regular grid 
    on a given potential temperature level. Same horizontal grid as in the ECMWF field.
    Input arguments:
        data : as read from read_ECMWF
        theta   : target potential temperature (K)
    """
    # Interpolate the ECMWF grid onto the target potential temperature
    Tg = np.empty(shape=[data.nlat,data.nlon])
    Pg = np.empty(shape=[data.nlat,data.nlon])
    for j in range(data.nlat):
        if j%10==0:
            print(j)
        for i in range(data.nlon):
			# - sign because interp wants growing abscissa
            #Tg[j,i] = np.interp(-theta,-data.var['PT'][:,j,i],data.var['T'][:,j,i])
            #Pg[j,i] = np.exp(np.interp(-theta,-data.var['PT'][:,j,i],np.log(data.var['P'][:,j,i])))
            # much longer but much better too
            f1T = interp1d(-data.var['PT'][:,j,i],data.var['T'][:,j,i],kind='cubic',assume_sorted=False)
            f1P = interp1d(-data.var['PT'][:,j,i],np.log(data.var['P'][:,j,i]),kind='cubic',assume_sorted=False)
            Tg[j,i] = f1T(-theta)
            Pg[j,i] = math.exp(f1P(-theta))
            #Tg[j,i] = krogh_interpolate(-data.var['PT'][:,j,i],data.var['T'][:,j,i],-theta)
            #Pg[j,i] = krogh_interpolate(-data.var['PT'][:,j,i],data.var['P'][:,j,i],-theta)      
    # Then do a 2d interpolation on this surface to the target grid        
    fT =  RectBivariateSpline(data.attr['lats'],data.attr['lons'],Tg)
    fP =  RectBivariateSpline(data.attr['lats'],data.attr['lons'],Pg)
    return fT, fP

#main
nsel = 1000    
# Case 1: SK-2
begSeq = 15000 
endSeq = 15300
date = datetime(2017,7,29)
theta = 410
dTheta = 2
name = 'SK-2'
# Case 2: SK-3
begSeq = 15001
endSeq = 15700
date = datetime(2017,7,29)
theta = 400
dTheta = 3
name = 'SK-3'

outputDir = '.'

#READ TRACK from file
#file = 'ucse/'+date.strftime('%y%m%d_1_ucse.nas')
data = STCinst('ucse',date)
utc = data.x
lats = data.var['Lat']
lons = data.var['Long']

# not needed but we keep it for future usage
#split the sequences into intervals with hourly boundaries
#nInt = 1+ int(endSeq/3600) - int(begSeq/3600)
#begIn = []
#endIn = []
#begIn.append(begSeq)
#for n in range(nInt):
#    if 3600*(int(begIn[n]/3600)+1) > endSeq:
#        endIn.append(endSeq)
#        break
#    else:
#        endIn.append(3600*(int(begIn[n]/3600)+1)-1)
#        begIn.append(endIn[-1]+1)
#
#print(nInt)
#print(begIn)
#print(endIn)

# Prepare the output file
stpd = date + timedelta(days=1)
part0 = {}
part0['lhead'] = 3
part0['outnfmt'] = 107
part0['mode'] =  3
part0['stamp_date'] = stpd.year*10**10 + stpd.month*10**8 + stpd.day*10**6 + stpd.hour*10**4 + stpd.minute*100
print('stamp_date',part0['stamp_date'])
part0['itime'] = 0
part0['step'] = 450
part0['idx_orgn'] = 1
part0['itime'] = 0
part0['nact_lastO'] = 0
part0['nact_lastNM'] = 0
part0['nact_lastNH'] = 0
part0['flag'] = np.empty(0,dtype=int)
part0['ir_start'] = np.empty(0,dtype=int)
part0['x'] = np.empty(0,dtype=float)
part0['y'] = np.empty(0,dtype=float)
part0['t'] = np.empty(0,dtype=float)
part0['p'] = np.empty(0,dtype=float)
part0['idx_back'] = np.empty(0,dtype=int)
numpart = 0

# Loop on the values of theta
for targetTheta in [theta-dTheta, theta, theta+dTheta]:
    # get the first time
    date1 = date + timedelta(hours=int(begSeq/3600))
    predat = ECMWF('STC',date1)
    predat._get_T()
    predat._mkp()
    predat._mkthet()
    (fT1,fP1) = interp3d_thet(predat,targetTheta)
    predat.close()
    
    #for n in range(nInt):
    # finds the boundary
    id1 = np.where(utc == begSeq)[0][0]
    id2 = np.where(utc == endSeq)[0][0]
    # get lat, lon and release time
    yy = data.var['Lat'][id1:id2+1]
    xx = data.var['Long'][id1:id2+1] 
    ir_start = utc[id1:id2+1] - 86400
    ir_start  = ir_start.astype(np.int)
    #weigthing factor
    w2 = (utc[id1:id2+1] % 3600) / 3600
    # Bloc_size
    block_size = len(xx)*nsel
    
    # get the next time
    date2 = date1 + timedelta(hours=1)
    posdat = ECMWF('STC',date2)
    posdat._get_T()
    posdat._mkp()
    posdat._mkthet()
    (fT2,fP2) = interp3d_thet(posdat,targetTheta)
    posdat.close()
        
    #%% get the interpolated value of TT and PP
    TT = (1-w2) * fT1.ev(yy,xx) + w2 * fT2.ev(yy,xx)
    PP = (1-w2) * fP1.ev(yy,xx) + w2 * fP2.ev(yy,xx)    
    
    # Check theta 
    THETA = TT *(cst.p0/PP)**cst.kappa 
    print(np.min(THETA-targetTheta),np.max(THETA-targetTheta))
    
    #replicate nsel times each field   
    part0['x'] = np.append(part0['x'],np.repeat(xx,nsel))
    part0['y'] = np.append(part0['y'],np.repeat(yy,nsel))
    part0['t'] = np.append(part0['t'],np.repeat(TT,nsel))
    part0['p'] = np.append(part0['p'],np.repeat(PP,nsel))
    part0['ir_start'] = np.append(part0['ir_start'],np.repeat(ir_start,nsel))
    # flag value is 14 + 0x10 + 0x20 = 62
    part0['flag'] = np.append(part0['flag'],np.full(block_size,62,dtype=int))
    idx1 = numpart+1
    numpart += block_size
    part0['idx_back'] = np.append(part0['idx_back'],np.arange(idx1,numpart+1,dtype=int))

# Terminate
part0['numpart'] = numpart
part0['nact'] = numpart

# write the result as part_000 file in the workin directory 
newpart0 = outputDir+'/'+name+'-part_000'
io107.writeidx107(newpart0,part0)


