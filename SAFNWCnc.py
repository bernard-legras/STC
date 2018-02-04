import numpy as np
import geosat
import os
import matplotlib.pyplot as plt
import re
from netCDF4 import Dataset



class SAFNWC_nc(geosat.GeoSat):
    
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
           
        filename='S_NWC_'+typ+nam+'_globeI-VISIR_'+date.strftime("%Y%m%d")+'T'+date.strftime("%H%M")+'00Z.nc'      
        print(filename)
        try:
            self.time = re.search(typ+'(.+?)'+'_globeI-VISIR',filename).group(1)             
        except AttributeError:
                # AAA, ZZZ not found in the original string
            self.time = '' # apply your error handling
        fullname = os.path.join(geosat.root_dir,sat,'safnwc',date.strftime("%Y"),
                            date.strftime("%Y_%m_%d"),filename)
        print (fullname)
        try: 
            self.ncid = Dataset(fullname, mode='r')
        except:
            print('No file ',filename)
            return
        self.swathnode = '/'
        self.mask=masksat
                
   def close(self):
        '''
        close an opened object.
        '''
        self.ncid.close()
                  

   def CT(self):
        if 'CT' not in self.var.keys():
            self._CT()
        return self
         

   def CTTH_PRESS(self):
        if 'CTTH_PRESS' not in self.var.keys():
            self._CTTH_PRESS(self)
        return self
        

   def CTTH_TEMPER(self):
        if 'CTTH_TEMPER' not in self.var.keys():
            self._CTTH_TEMPER(self)
        return self

    
   def _merge(self,other):
       for var in other.var.keys():
           self.var[var] = other.var[var]
       for atr in other.attr.keys():
           self.attr[atr] = other.attr[atr]



class SAFNWC_CT_nc(SAFNWC_nc):
    '''
    Class to read SAF NWC Cloud Type products
    '''

    def __init__(self,date,sat):
        typ='CT_'
        SAFNWC_nc.__init__(self,date,sat,typ)  
        self.var={}
        self.attr={}   
        self.sat=sat
        
        
    def _CT(self):
        '''
        returns numpy array containing the Cloud Type product from SAF NWC file
        '''
        self.attr['CT']={}
        data = self.ncid.variables['ct'][:]
        self.var['CT']=np.ma.array(data)
        self.var['CT'].__setmask__(self.mask)
        self.attr['CT']['units']=('1: Cloud free land','2: Cloud free sea','3: Snow over land','4: Snow over sea','5: Very low clouds',\
        '6: Low Clouds','7: Mid-level clouds','8: High opaque clouds','9: Very high opaque clouds',\
        '10: Fractional clouds','11: High semitransparent thin clouds','12: High semitransparent meanly thick clouds','13: High semitransparent thick clouds',\
        '14: High semitransparent above low or medium clouds','15: High semitransparent above snow/ice ')
    
        data = self.ncid.variables['ct_pal'][:]
        self.attr['CT']['PALETTE']=data
        return

 
        



class SAFNWC_CTTH_nc(SAFNWC_nc):
    '''
    Class to read SAF NWC Cloud Type products
    '''

    def __init__(self,date,sat):
        typ='CTTH_'
        SAFNWC_nc.__init__(self,date,sat,typ)  
        self.var={}
        self.attr={}   
        self.sat=sat
                
    def _CTTH_PRESS(self):
        self.attr['CTTH_PRESS']={}
        data = self.ncid.variables['ctth_pres'][:]
        print("MAX press",np.amax(data))
        #gain=25 #hPa/count
        #intercept=-250 #hPa
        gain=10 #hPa/count
        intercept=0 #hPa
        fillvalue=65535.
        #self.var['CTTH_PRESS']=np.ma.array(gain*var+intercept)
        #self.var['CTTH_PRESS']=np.ma.array(gain*data+intercept)
        self.var['CTTH_PRESS']=np.ma.array(data/100)
        self.var['CTTH_PRESS'].__setmask__(self.mask)
        self.var['CTTH_PRESS']._sharedmask=False
        self.attr['CTTH_PRESS']['units']='hPa'
        data = self.ncid.variables['ctth_pres_pal'][:]        

        self.attr['CTTH_PRESS']['PALETTE']=data
        self.attr['CTTH_PRESS']['gain'] = gain
        self.attr['CTTH_PRESS']['intercept'] = intercept
        self.attr['CTTH_PRESS']['FillValue']=fillvalue

        return


    def _CTTH_TEMPER(self):
        self.attr['CTTH_TEMPER']={}
        data = self.ncid.variables['ctth_tempe'][:]
        print("MAX temp",np.amax(data))
        #gain=1 #K/count
        #intercept=150 #K
        gain=0.01 #K/count
        intercept=304 #K
        fillvalue=65535.
        self.var['CTTH_TEMPER']=np.ma.array(data)
        self.var['CTTH_TEMPER'].__setmask__(self.mask)
        self.var['CTTH_TEMPER']._sharedmask=False
        self.attr['CTTH_TEMPER']['units']='K'
        self.attr['CTTH_TEMPER']['FillValue']=-1
        data = self.ncid.variables['ctth_tempe_pal'][:]

        self.attr['CTTH_TEMPER']['PALETTE']=data
        self.attr['CTTH_TEMPER']['gain'] = gain
        self.attr['CTTH_TEMPER']['intercept'] = intercept
        return 
        
 
        
 
           
        
