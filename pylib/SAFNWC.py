#!/usr/bin/env python
# encoding: utf-8

'''

Module for reading SAFNWC file hdf5, both himawari and msg3
Supported products : Cloud Mask (CMa), Cloud Type (CT), Cloud Top Temperature and Height (CTTH)

S. Bucci - Created on Tue Apr 11 17:21:26 CET 2017

'''
import tables
import numpy as np
#from datetime import datetime
#path='/net/grapelli/limbo/data/sats/'
#path='/home/sbucci/Data'
#from geosat import *
import geosat
import os
#import matplotlib.colors as mcl
#import pprint
#import pylab
import matplotlib.pyplot as plt
import re
class SAFNWC(geosat.GeoSat):
    
   def __init__(self,date,sat,typ):
     
        '''
        filename : name of the SAFNWC file
        '''
        if sat=='himawari':
            nam='hima08'
            let='J'      
            geosat.read_mask_himawari()
            masksat=geosat.mask_sat['himawari']
            
        elif sat=='msg3':
            nam='MSG3'
            let='M'
            geosat.read_mask_MSG()        
            masksat=geosat.mask_sat['msg']
   
        elif sat=='msg1':
            nam='MSG1'
            let='I'
            geosat.read_mask_MSG()        
            masksat=geosat.mask_sat['msg']
   
        filename='SAFNWC_'+nam+'_'+typ+date.strftime("%Y%m%d%H%M")+'_globe'+let+'______.h5'
        #print(filename)
        try:
            self.time = re.search(typ+'(.+?)'+'_globe',filename).group(1)
        except AttributeError:
                # AAA, ZZZ not found in the original string
            self.time = '' # apply your error handling
        fullname = os.path.join(geosat.root_dir,sat,'safnwc',date.strftime("%Y"),
                            date.strftime("%Y_%m_%d"),filename)
        print (fullname)
        try: 
            self.h5 = tables.open_file(fullname, mode='r')
        except:
            print('No file ',filename)
            return
        self.swathnode = '/'
        self.mask=masksat
                
   def close(self):
        '''
        close an opened object.
        '''
        self.h5.close()
        
           
   def CMa(self):
        if 'CMa' not in self.var.keys():
            self._CMa()
            del self.mask
        return self
          
   def CMa_QUALITY(self):
        if 'CMa_QUALITY' not in self.var.keys():
            self._CMa_QUALITY()
        return self
   
   def CMa_TEST(self):
        if 'CMa_TEST' not in self.var.keys():
            self._CMa_TEST()
        return self
    
   def CMa_DUST(self):
        if 'CMa_DUST' not in self.var.keys():
            self._CMa_DUST()
        return self
    
   def CMa_VOLCANIC(self):
        if 'CMa_VOLCANIC' not in self.var.keys():
            self._CMa_VOLCANIC()
        return self

   def CT(self):
        if 'CT' not in self.var.keys():
            self._CT()
        return self
        
   def CT_QUALITY(self):
        if 'CT_QUALITY' not in self.var.keys():
            self._CT_QUALITY()
        return self

   def CT_PHASE(self):
        if 'CT_PHASE' not in self.var.keys():
            self._CT_PHASE()
        return self

   def CTTH_PRESS(self):
        if 'CTTH_PRESS' not in self.var.keys():
            self._CTTH_PRESS(self)
        return self
        
   def CTTH_QUALITY(self):
        if 'CTTH_QUALITY' not in self.var.keys():
            self._CTTH_QUALITY()
        return self

   def CTTH_TEMPER(self):
        if 'CTTH_TEMPER' not in self.var.keys():
            self._CTTH_TEMPER(self)
        return self

   def CTTH_HEIGHT(self):
        if 'CTTH_HEIGHT' not in self.var.keys():
            self._CTTH_HEIGHT(self)
        return self

   def CTTH_EFFECT(self):
        if 'CTTH_EFFECT' not in self.var.keys():
            self._CTTH_EFFECT(self)
        return self
    
   def _merge(self,other):
       for var in other.var.keys():
           self.var[var] = other.var[var]
       for atr in other.attr.keys():
           self.attr[atr] = other.attr[atr]
        
class SAFNWC_CMa(SAFNWC):
    '''
    Class to read SAF NWC Cloud mask products
    '''
    
    def __init__(self,date,sat):
        typ='CMa__'  
        SAFNWC.__init__(self,date,sat,typ)  
        self.var={}

        self.attr={}   
        self.sat=sat
        

    def _CMa(self):
        '''
        returns numpy array containing the product from SAF NWC file
        '''
        self.attr['CMa']={}
        data = self.h5.get_node(self.swathnode, 'CMa')
        data = data.read()
        self.var['CMa'] = np.ma.array(data)
        self.var['CMa'].__setmask__(self.mask)
        self.var['CMa']._sharedmask=False
        self.attr['CMa']['units']=('0: Non processed','1: Cloud free','2: Cloud contaminated','3: Cloud filled','4: Snow/Ice contaminated','5: Undefined')
        data = self.h5.get_node(self.swathnode, '01-PALETTE')
        self.attr['CMa']['PALETTE']=data.read()
        
        return #self.var['CMa']
        
    
    def _CMa_QUALITY(self):
        '''
        returns numpy array containing the quality information on the product from SAF NWC file
        '''
        self.attr['CMa_QUALITY']={}
        data = self.h5.get_node(self.swathnode, 'CMa_QUALITY')
        var=data.read()
        q=np.zeros((6,len(var),len(var)))
        q[0,:,:]=var&0x7
        q[1,:,:]=var&0x18 >>3 
        q[2,:,:]=var&0x60 >>5
        q[3,:,:]=var&0x180>>7
        q[4,:,:]=var&0x200>>9
        q[5,:,:]=var&0x400>>10
        m=[self.mask,self.mask,self.mask,self.mask,self.mask,self.mask]
        self.var['CMa_QUALITY']=np.ma.array(q,mask=m)
        self.var['CMa_QUALITY']._sharedmask=False
        self.attr['CMa_QUALITY']['units']=('0: illumination and viewing conditions','1: NWP input data status','2: SEVIRI input data status',\
        '3: Quality of the processing','4: Temporal processing indicator','5: HRV processing indicator')
        self.attr['CMa_QUALITY']['flags']=(('0: Undefined (Space)','1: Night','2: Twilight','3: Day','4: Sunglint'),\
                   ('0: Undefined (Space)','1: All NWP parameters available (no low level inversion)',\
                    '2: All NWP parameters available (low level inversion)','3: At least one NWP parameter missing'),\
                    ('0: Undefined (Space)','1: All useful SEVIRI channels available','2: At least one useful SEVIRI channel missing',\
                    '3: At least one mandatory SAVIRI channel missing'),\
                     ('0: Non processed','1: Good quality','2: Poor quality','3: Reclassified after spatial smoothing'),\
                     ('0: Not performed','1: Performed'),\
                     ('0: Not performed','1: Performed'))
        return
        
    def _CMa_TEST(self):
        self.attr['CMa_TEST']={}
        data = self.h5.get_node(self.swathnode, 'CMa_TEST')
        var = data.read()
        t=np.zeros((16,len(var),len(var)))
        t[0,:,:]=var&0x1
        t[1,:,:]=var&(0x2)>>1
        t[2,:,:]=var&(0x4)>>2
        t[3,:,:]=var&(0x8)>>3
        t[4,:,:]=var&(0x10)>>4
        t[5,:,:]=var&(0x20)>>5
        t[6,:,:]=var&(0x40)>>6
        t[7,:,:]=var&(0x80)>>7
        t[8,:,:]=var&(0x100)>>8
        t[9,:,:]=var&(0x200)>>9
        t[10,:,:]=var&(0x400)>>10
        t[11,:,:]=var&(0x800)>>11
        t[12,:,:]=var&(0x1000)>>12
        t[13,:,:]=var&(0x2000)>>13
        t[14,:,:]=var&(0x4000)>>14
        t[15,:,:]=var&(0x8000)>>15
        m=[self.mask,self.mask,self.mask,self.mask,self.mask,\
           self.mask,self.mask,self.mask,self.mask,self.mask,\
           self.mask,self.mask,self.mask,self.mask,self.mask]
        self.var['CMa_TEST']=np.ma.array(t,mask=m)
        self.var['CMa_TEST']._sharedmask=False
        self.attr['CMa_TEST']['units']=('0: T10.8um or SST','1: R0.6um (land) or R0.8um (sea)','2: Sunglint test using 3.9 um',\
        '3: Local Spatial Texture','4: T10.8um-T12.0um','5: T10.8um-T3.9um or T12.0um-T3.9 um','6: T3.9um-T10.8um',\
        '7: Spatial smoothing (reclassify isolated cloud-free pixels)','8: T8.7um-T3.9um','9: R1.6um (sea)',\
        '10: T8.7um-T10.8um or T10.8um-T8.7um','11: Snow using R1.6um or T3.9um','12: HRV-based test','13: Stationary cloud in twilight',\
        '14: Spatial expansion of stationary cloud in twilight','15: Temporal-differencing')
        return
                
    def _CMa_DUST(self):
        self.attr['CMa_DUST']={}
        data = self.h5.get_node(self.swathnode, 'CMa_DUST')
        self.var['CMa_DUST']=np.ma.array(data.read())
        self.var['CMa_DUST'].__setmask__(self.mask)
        self.var['CMa_DUST']._sharedmask=False
        self.attr['CMa_DUST']['units']=('0: Non processed','1: Dust','2: Non dust','3: Undefined (separability problems)')
        data = self.h5.get_node(self.swathnode, '02-PALETTE')
        self.attr['CMa_DUST']['PALETTE']=data.read()
        
        return 
    def _CMa_VOLCANIC(self):
        self.attr['CMa_VOLCANIC']={}
        data = self.h5.get_node(self.swathnode, 'CMa_VOLCANIC')
        self.var['CMa_VOLCANIC']=np.ma.array(data.read())
        self.var['CMa_VOLCANIC'].__setmask__(self.mask)
        self.var['CMa_VOLCANIC']._sharedmask=False
        self.attr['CMa_VOLCANIC']['units']=('0: Non processed','1: Volcanic plume','2: Non volcanic plume','3: Undefined (separability problems)')
        data = self.h5.get_node(self.swathnode, '03-PALETTE')
        self.attr['CMa_VOLCANIC']['PALETTE']=data.read()
        return 
        
    
class SAFNWC_CT(SAFNWC):
    '''
    Class to read SAF NWC Cloud Type products
    '''

    def __init__(self,date,sat):
        typ='CT___'
        SAFNWC.__init__(self,date,sat,typ)  
        self.var={}
        self.attr={}   
        self.sat=sat
        
        
    def _CT(self):
        '''
        returns numpy array containing the Cloud Type product from SAF NWC file
        '''
        self.attr['CT']={}
        data = self.h5.get_node(self.swathnode, 'CT')
        self.var['CT']=np.ma.array(data.read())
        self.var['CT'].__setmask__(self.mask)
        #self.var['CT']._sharedmask=False
#        self.var['CT']=np.ma.array(ct,mask=m)
#        self.var['CT']._sharedmask=False
        self.attr['CT']['units']=('0: Non processed','1: Cloud free land','2: Cloud free sea','3: land+snow','4: sea+snow','5: very low cumuliform',\
        '6: very low stratiform','7: low cumulifrom','8: low stratiform','9: medium cumuliform','10: medium stratiform',\
        '11: high opaque cumuliform','12: high opaque stratiform','13: very high opaque cumuliform','14: very high opaque stratiform',\
        '15: high semitransparent thin clouds','16: high semitransparent meanly thick clouds','17: high semitransparent thick clouds',\
        '18: high semitransparent above low or medium clouds','19: fractional clouds','20: undefined')
        data = self.h5.get_node(self.swathnode, '01-PALETTE')
        self.attr['CT']['PALETTE']=data.read()
        return

    def _CT_QUALITY(self):
        self.attr['CT_QUALITY']={}
        data = self.h5.get_node(self.swathnode, 'CT_QUALITY')
        var = data.read()
        q=np.zeros((5,len(var),len(var)))
        q[0,:,:]=var&0x7
        q[1,:,:]=(var&0x18)>>3
        q[2,:,:]=(var&0x60)>>5
        q[3,:,:]=(var&0x180)>>7
        q[4,:,:]=(var&0x200)>>9
        m=[self.mask,self.mask,self.mask,self.mask,self.mask]
        self.var['CT_QUALITY']=np.ma.array(q,mask=m)
        self.var['CT_QUALITY']._sharedmask=False
        self.attr['CT_QUALITY']['units']=('0: Illumination and viewing conditions','1: NWP input data status','2: SEVIRI input data status',\
        '3: Quality of the processing','4: Existence of separation between stratiform and cumuliform')
        self.attr['CT_QUALITY']['flags']=(('0: Undefined (Space)','1: Night','2: Twilight','3: Day','4: Sunglint'),\
                   ('0: Undefined (Space)','1: All NWP parameters available (no low level inversion)',\
                    '2: All NWP parameters available (low level inversion)','3: At least one NWP parameter missing'),\
                    ('0: Undefined (Space)','1: All useful SEVIRI channels available','2: At least one useful SEVIRI channel missing',\
                    '3: At least one mandatory SAVIRI channel missing'),\
                     ('0: Non processed','1: Good quality','2: Poor quality','3: Reclassified after spatial smoothing'),\
                     ('0: Not performed','1: Performed'))
        return
     
   

    def _CT_PHASE(self):
        self.attr['CT_PHASE']={}
        data = self.h5.get_node(self.swathnode, 'CT_PHASE')
        self.var['CT_PHASE']=np.ma.array(data.read())
        self.var['CT_PHASE'].__setmask__(self.mask)
        self.var['CT_PHASE']._sharedmask=False
        self.attr['CT_PHASE']['units']=('0: Non processed or no cloud','1: water cloud','2: ice cloud','3: undefined')
        data = self.h5.get_node(self.swathnode, '02-PALETTE')
        self.attr['CT_PHASE']['PALETTE']=data.read()

        return 
        
class SAFNWC_CTTH(SAFNWC):
    '''
    Class to read SAF NWC Cloud Type products
    '''

    def __init__(self,date,sat):
        typ='CTTH_'
        SAFNWC.__init__(self,date,sat,typ)  
        self.var={}
        self.attr={}   
        self.sat=sat
                
    def _CTTH_PRESS(self):
        self.attr['CTTH_PRESS']={}
        data = self.h5.get_node(self.swathnode, 'CTTH_PRESS')
        var= data.read()
        var=(var&0x3FF).astype(int)
        gain=25 #hPa/count
        intercept=-250 #hPa
        self.var['CTTH_PRESS']=np.ma.array(gain*var+intercept)
        self.var['CTTH_PRESS'].__setmask__(self.mask)
        self.var['CTTH_PRESS']._sharedmask=False
        self.attr['CTTH_PRESS']['units']='hPa'
        data = self.h5.get_node(self.swathnode, '01-PALETTE')
        self.attr['CTTH_PRESS']['PALETTE']=data.read()
        self.attr['CTTH_PRESS']['gain'] = gain
        self.attr['CTTH_PRESS']['intercept'] = intercept

        return

    def _CTTH_QUALITY(self):
        self.attr['CTTH_QUALITY']={}
        data = self.h5.get_node(self.swathnode, 'CTTH_QUALITY')
        var=data.read()
        q=np.empty([6,len(var),len(var)],dtype=np.uint8)
        q[0,:,:]=var&0x3
        q[1,:,:]=(var&0x4)>>2 
        q[2,:,:]=(var&0x38)>>3
        q[3,:,:]=(var&0xc0)>>6
        q[4,:,:]=(var&0xf00)>>8
        q[5,:,:]=(var&0x3000)>>12
        m=[self.mask,self.mask,self.mask,self.mask,self.mask,self.mask]
        self.var['CTTH_QUALITY']=np.ma.array(q,mask=m)
        self.var['CTTH_QUALITY']._sharedmask=False
        self.attr['CTTH_QUALITY']['units']=('0: processing status','1: RTTOV IR availability','2: NWP input data status','3: SEVIRI input data status','4: method','5: Quality of processing')
        self.attr['CTTH_QUALITY']['flags']=(('0: Non processed','1: Non processed beacause cloud free','2: Processed, cloudy but no results','3: Processed with results'),\
                   ('0: Non available','1: Available'),\
                    ('0: Undefined (Space)','1: All NWP parameters available (no thermal inversion)','2: All NWP parameters available (thermal inversion present)',\
                    '3: Some NWP pressure levels missing (no thermal inversion)','4: Some NWP pressure levels missing (thermal inversion present)',\
                    '5: At least one mandatory NWP information missing'),\
                    ('0: Undefined (Space)','1: All useful SEVIRI channels available','2: At least one useful SEVIRI channel missing',\
                    '3: At least one mandatory SAVIRI channel missing'),\
                     ('0: Non processed','1: Opaqe cloud, using RTTOV','2: Opaqe cloud, not using RTTOV','3: Intercept method 10.8um/13.4um','4: Intercept method 10.8um/6.2um',\
                      '5: Intercept method 10.8um/7.3um','6: Radiance Ratioing method 10.8um/13.4um','7: Radiance Ratioing method 10.8um/6.2um','8: Radiance Ratioing method 10.8um/7.3um',\
                      '9: Spare','10: Spare','11: Spare','12: Spare','13: Opaqe cloud, using RTTOV, in case thermal inversion',\
                      '14: Spatial smoothing (gap filling in semi-transparent cloud field)','15: Spare for not yet defined methods'),\
                     ('0: No results (Non-processed, cloud free, no reliable method)','1: Good quality','2: Poor quality'))
        return

    def _CTTH_TEMPER(self):
        self.attr['CTTH_TEMPER']={}
        data = self.h5.get_node(self.swathnode, 'CTTH_TEMPER')
        var= data.read()
        var=(var&0xFF).astype(int) 
        gain=1 #K/count
        intercept=150 #K
        self.var['CTTH_TEMPER']=np.ma.array(gain*var+intercept)
        self.var['CTTH_TEMPER'].__setmask__(self.mask)
        self.var['CTTH_TEMPER']._sharedmask=False
        self.attr['CTTH_TEMPER']['units']='K'
        data = self.h5.get_node(self.swathnode, '03-PALETTE')
        self.attr['CTTH_TEMPER']['PALETTE']=data.read()
        self.attr['CTTH_TEMPER']['gain'] = gain
        self.attr['CTTH_TEMPER']['intercept'] = intercept
        return 
        
    def _CTTH_HEIGHT(self):
        self.attr['CTTH_HEIGHT']={}
        data = self.h5.get_node(self.swathnode, 'CTTH_HEIGHT')
        var = data.read()
        var=(var&0x7F).astype(int)
        gain=200 #m/count
        intercept=-2000 #m
        self.var['CTTH_HEIGHT']=np.ma.array(gain*var+intercept)
        self.var['CTTH_HEIGHT'].__setmask__(self.mask)
        self.var['CTTH_HEIGHT']._sharedmask=False
        self.attr['CTTH_HEIGHT']['units']='m'
        data = self.h5.get_node(self.swathnode, '02-PALETTE')
        self.attr['CTTH_HEIGHT']['PALETTE']=data.read()       
        self.attr['CTTH_HEIGHT']['gain'] = gain
        self.attr['CTTH_HEIGHT']['intercept'] = intercept
        return 
        
    def _CTTH_EFFECT(self):
        self.attr['CTTH_EFFECT']={}
        data = self.h5.get_node(self.swathnode, 'CTTH_EFFECT')
        var = data.read()
        var=(var&0x1F).astype(int)
        gain=5 #%/count
        intercept=-50 #%
        self.var['CTTH_EFFECT']=np.ma.array(gain*var+intercept)
        self.var['CTTH_EFFECT'].__setmask__(self.mask)
        self.var['CTTH_EFFECT']._sharedmask=False
        self.attr['CTTH_EFFECT']['units']='%'
        data = self.h5.get_node(self.swathnode, '04-PALETTE')
        self.attr['CTTH_EFFECT']['PALETTE']=data.read()
        self.attr['CTTH_EFFECT']['gain'] = gain
        self.attr['CTTH_EFFECT']['intercept'] = intercept
        return
        
   
    #Example;
    #To have the method (index 4) used for CTTH in the grid point (300,300)
    #qf[4][q[4,3100,300]]
    #'Opaqe cloud, using RTTOV'
    
