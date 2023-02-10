import numpy as np
import geosat
import os
import re
from netCDF4 import Dataset

class SAFNWC(geosat.PureSat):
    
    def __init__(self,date,sat,typ,BBname=None,version=None):
     
        '''
        filename : name of the SAFNWC file
        '''
        self.date = date
        geosat.PureSat.__init__(self,sat)
        # switch the rootdir if SAFBox
        safnwc = 'safnwc'
        root_dir = geosat.root_dir
            
        if 'hima' in sat :
            sat = 'himawari'
            nam='HIMA08'
            region='globeJ'
            geosat.read_mask_himawari()
            masksat=geosat.mask_sat['himawari']
            VISIR = 'NR'
        elif sat=='msg3':
            nam='MSG3'
            region='globeM'
            geosat.read_mask_MSG()        
            masksat=geosat.mask_sat['msg']
            VISIR = 'VISIR'
        elif sat=='msg1':
            nam='MSG1'
            region='globeI'
            VISIR = 'VISIR'
            geosat.read_mask_MSG()        
            masksat=geosat.mask_sat['msg']
        if BBname == 'SAFBox':
            safnwc = 'safnwc-SAFBox'
            root_dir = geosat.alt_root_dir
            region = 'FULLAMA'
            # Meaning of the following lines
            if ('hima' in sat) & (version is None):
                nam='HIMAWARI08'
        if version is not None:
            safnwc = safnwc + '-' + version
        
        filename='S_NWC_'+typ+'_'+nam+'_'+region+'-'+VISIR+'_'+date.strftime("%Y%m%d")+'T'+date.strftime("%H%M")+'00Z.nc'      
        
        #print(filename)
        try:
            self.time = re.search(typ+'(.+?)'+'_globeI-VISIR',filename).group(1)             
        except AttributeError:
                # AAA, ZZZ not found in the original string
            self.time = '' # apply your error handling
                         
        #fullname = os.path.join(root_dir,sat,safnwc,'netcdf',date.strftime("%Y"),
        #                    date.strftime("%Y_%m_%d"),filename)
        fullname = os.path.join(root_dir,sat,safnwc,date.strftime("%Y"),
                            date.strftime("%Y_%m_%d"),filename)
        #print (fullname)
        try: 
            self.ncid = Dataset(fullname, mode='r')
        except AttributeError:
            print('No file ',filename)
            return
        # The mask must be truncated when BB is not None
        self.mask=masksat
        if BBname == 'SAFBox':
            BB = {'msg1':[[341,307],[1857,3351]],'himawari':[[475,445],[2751,3794]]} [sat]
            self.mask = masksat[BB[0][0]:BB[1][0]+1,BB[0][1]:BB[1][1]+1]
        self.nx = self.ncid.dimensions['nx'].size
        self.ny = self.ncid.dimensions['ny'].size
                
    def close(self):
        '''
        close an opened object.
        '''
        self.ncid.close()
        
    def _get_var(self,var):
        if var not in self.ncid.variables.keys():
            print('variable ',var,' is not available')
            return
        self.var[var] = np.ma.array(data=self.ncid.variables[var][:])
        self.var[var].__setmask__(self.mask)
        self.var[var]._sharedmask=False
        # This masking is used for the transformation to the latxlon grid
        # Therefore masking the filled pixels must be done at later stage
        self.attr[var] = {}
        self.attr[var]['FillValue'] = self.ncid.variables[var]._FillValue
        # Indicate that filled value have not been masked
        self.attr[var]['MaskedFill'] = False
        try:
            self.attr[var]['gain'] = self.ncid.variables[var].scale_factor
            self.attr[var]['intercept'] = self.ncid.variables[var].add_offset
        except:
            pass
        try:
            self.attr[var]['units'] = self.ncid.variables[var].units
        except:
            self.attr[var]['units'] = None
        try:
            self.attr[var]['PALETTE'] = self.ncid.variables[var+'_pal'][:]
        except:
            pass
    
    def _merge(self,other):
       for var in other.var.keys():
           self.var[var] = other.var[var]
       for atr in other.attr.keys():
           self.attr[atr] = other.attr[atr]

class SAFNWC_CMa(SAFNWC):
    '''
    Class to read SAF NWC Cloud mask products
    '''
    
    def __init__(self,date,sat,BBname=None,version=None):
        typ='CMA'  
        SAFNWC.__init__(self,date,sat,typ,BBname=BBname,version=version)
        self.typ = typ
        self.var={}
        self.attr={}   
        self.sat=sat
        
    def _CMa(self):
        '''
        returns numpy array containing the product from SAF NWC file
        '''
        self.attr['CMa']={}
        self.var['CMa'] = np.ma.array(self.ncid.variables['cma'][:])
        self.var['CMa'].__setmask__(self.mask)
        self.var['CMa']._sharedmask=False
        self.attr['CMa']['units']=('0: Non processed','1: Cloud free','2: Cloud contaminated','3: Cloud filled','4: Snow/Ice contaminated','5: Undefined')
        self.attr['CMa']['PALETTE'] = self.ncid.variables['cma_pal'][:]      
        return #self.var['CMa']
        
    def _CMa_QUALITY(self):
        '''
        returns numpy array containing the quality information on the product from SAF NWC file
        '''
        self.attr['CMa_QUALITY']={}
        var = self.ncid.variables['cma_quality'][:]
        q=np.zeros((6,self.ny,self.nx))
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
        var = self.ncid.variables['cma_test'][:]
        t=np.zeros((16,self.ny,self.nx))
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
        data = self.ncid.variables['cma_dust'][:]
        self.var['CMa_DUST']=np.ma.array(data.read())
        self.var['CMa_DUST'].__setmask__(self.mask)
        self.var['CMa_DUST']._sharedmask=False
        self.attr['CMa_DUST']['units']=('0: Non processed','1: Dust','2: Non dust','3: Undefined (separability problems)')
        self.attr['CMa_DUST']['PALETTE'] = self.ncid.variables['cma_dust-pal'][:]
        
        return 
    def _CMa_VOLCANIC(self):
        self.attr['CMa_VOLCANIC']={}
        self.var['CMa_VOLCANIC'] = np.ma.array(self.ncid.variables['cma_volcanic'][:])
        self.var['CMa_VOLCANIC'].__setmask__(self.mask)
        self.var['CMa_VOLCANIC']._sharedmask=False
        self.attr['CMa_VOLCANIC']['units']=('0: Non processed','1: Volcanic plume','2: Non volcanic plume','3: Undefined (separability problems)')
        self.attr['CMa_VOLCANIC']['PALETTE'] = self.ncid.variables['cma_volcanic_pal'][:]
        return 
        

class SAFNWC_CT(SAFNWC):
    '''
    Class to read SAF NWC Cloud Type products
    '''

    def __init__(self,date,sat,BBname=None,version=None):
        typ='CT'
        SAFNWC.__init__(self,date,sat,typ,BBname=BBname,version=version)
        self.typ = typ
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

    def _CT_QUALITY(self):
        self.attr['CT_QUALITY']={}
        var = self.ncid.variables['CT_QUALITY']
        q=np.zeros((5,self.ny,self.nx))
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
        self.var['CT_PHASE']=np.ma.array(self.ncid.variables['ct_phase'][:])
        self.var['CT_PHASE'].__setmask__(self.mask)
        self.var['CT_PHASE']._sharedmask=False
        self.attr['CT_PHASE']['units']=('0: Non processed or no cloud','1: water cloud','2: ice cloud','3: undefined')
        self.attr['CT_PHASE']['PALETTE'] = self.ncid.variables['ct_phase_pal']        

class SAFNWC_CTTH(SAFNWC):
    '''
    Class to read SAF NWC Cloud Type products
    '''

    def __init__(self,date,sat,BBname=None,version=None):
        typ='CTTH'
        SAFNWC.__init__(self,date,sat,typ,BBname=BBname,version=version)
        self.typ = typ
        self.var={}
        self.attr={}   
        self.sat=sat
                
    def _CTTH_PRESS(self):
        ''' Reads CTTH pressure 
            Beware that the pressure values are converted to hPa so that the Fillvalue
            63350 is now outside the range of eligible values '''
        self.attr['CTTH_PRESS']={}
        self.var['CTTH_PRESS'] = np.ma.array(data = self.ncid.variables['ctth_pres'][:]/100)
        # mask is replaced by the standard sat mask in order to allow sat_togrid conversion   
        self.var['CTTH_PRESS'].__setmask__(self.mask)
        self.var['CTTH_PRESS']._sharedmask=False
        self.attr['CTTH_PRESS']['units']='hPa'
        self.attr['CTTH_PRESS']['PALETTE'] = self.ncid.variables['ctth_pres_pal'][:]
        self.attr['CTTH_PRESS']['gain'] = self.ncid.variables['ctth_pres'].scale_factor
        self.attr['CTTH_PRESS']['intercept'] = self.ncid.variables['ctth_pres'].add_offset
        self.attr['CTTH_PRESS']['FillValue'] = self.ncid.variables['ctth_pres']._FillValue
        self.attr['CTTH_PRESS']['MaskedFill'] = False
        return

    def _CTTH_TEMPER(self):
        self.attr['CTTH_TEMPER']={}
        self.var['CTTH_TEMPER'] = np.ma.array(data = self.ncid.variables['ctth_tempe'][:])
        self.var['CTTH_TEMPER'].__setmask__(self.mask)
        self.var['CTTH_TEMPER']._sharedmask=False
        self.attr['CTTH_TEMPER']['units'] = self.ncid.variables['ctth_tempe'].units
        self.attr['CTTH_TEMPER']['PALETTE'] = self.ncid.variables['ctth_tempe_pal'][:]
        self.attr['CTTH_TEMPER']['gain'] = self.ncid.variables['ctth_tempe'].scale_factor
        self.attr['CTTH_TEMPER']['intercept'] = self.ncid.variables['ctth_tempe'].add_offset
        self.attr['CTTH_TEMPER']['FillValue'] = self.ncid.variables['ctth_tempe']._FillValue
        self.attr['CTTH_TEMPER']['MaskedFill'] = False
        return 
    
    def _CTTH_METHOD(self):
        self.attr['CTTH_METHOD']={}
        self.var['CTTH_METHOD'] = np.ma.array(data = self.ncid.variables['ctth_method'][:].astype(np.uint16))
        self.var['CTTH_METHOD'].__setmask__(self.mask)
        self.var['CTTH_METHOD']._sharedmask=False
        self.attr['CTTH_METHOD']['flag_mask'] = self.ncid.variables['ctth_method'].flag_mask
        self.attr['CTTH_METHOD']['flag_values'] = self.ncid.variables['ctth_method'].flag_values
        self.attr['CTTH_METHOD']['flag_meanings'] = self.ncid.variables['ctth_method'].flag_meanings.split()
        return
    
    def _CTTH_STATUS(self):
        self.attr['CTTH_STATUS']={}
        self.var['CTTH_STATUS'] = np.ma.array(data = self.ncid.variables['ctth_status_flag'][:].astype(np.uint16))
        self.var['CTTH_STATUS'].__setmask__(self.mask)
        self.var['CTTH_STATUS']._sharedmask=False
        self.attr['CTTH_STATUS']['flag_mask'] = self.ncid.variables['ctth_status_flag'].flag_mask
        self.attr['CTTH_STATUS']['flag_values'] = self.ncid.variables['ctth_status_flag'].flag_values
        self.attr['CTTH_STATUS']['flag_meanings'] = self.ncid.variables['ctth_status_flag'].flag_meanings.split()
        return
    
    def _CTTH_QUALITY(self):
        self.attr['CTTH_QUALITY']={}
        self.var['CTTH_QUALITY'] = np.ma.array(data = self.ncid.variables['ctth_quality'][:].astype(np.uint16))
        self.var['CTTH_QUALITY'].__setmask__(self.mask)
        self.var['CTTH_QUALITY']._sharedmask=False
        self.attr['CTTH_QUALITY']['flag_mask'] = self.ncid.variables['ctth_quality'].flag_mask
        self.attr['CTTH_QUALITY']['flag_values'] = self.ncid.variables['ctth_quality'].flag_values
        self.attr['CTTH_QUALITY']['flag_meanings'] = self.ncid.variables['ctth_quality'].flag_meanings.split()
        return
    
    def _CTTH_CONDITIONS(self):
        self.attr['CTTH_CONDITIONS']={}
        self.var['CTTH_CONDITIONS'] = np.ma.array(data = self.ncid.variables['ctth_conditions'][:].astype(np.uint16))
        self.var['CTTH_CONDITIONS'].__setmask__(self.mask)
        self.var['CTTH_CONDITIONS']._sharedmask=False
        self.attr['CTTH_CONDITIONS']['flag_mask'] = self.ncid.variables['ctth_conditions'].flag_mask
        self.attr['CTTH_CONDITIONS']['flag_values'] = self.ncid.variables['ctth_conditions'].flag_values
        self.attr['CTTH_CONDITIONS']['flag_meanings'] = self.ncid.variables['ctth_conditions'].flag_meanings.split()
        return

    def _CTTH_HEIGHT(self):       
        self.attr['CTTH_HEIGHT']={}
        self.var['CTTH_HEIGHT']=np.ma.array(data=self.ncid.variables['ctth_alti'][:])
        self.var['CTTH_HEIGHT'].__setmask__(self.mask)
        self.var['CTTH_HEIGHT']._sharedmask=False
        self.attr['CTTH_HEIGHT']['units'] = self.ncid.variables['ctth_alti'].units     
        self.attr['CTTH_HEIGHT']['PALETTE'] = self.ncid.variables['ctth_alti_pal'][:]
        self.attr['CTTH_HEIGHT']['gain'] = self.ncid.variables['ctth_alti'].scale_factor
        self.attr['CTTH_HEIGHT']['intercept'] = self.ncid.variables['ctth_alti'].add_offset
        self.attr['CTTH_HEIGHT']['FillValue'] = self.ncid.variables['ctth_alti']._FillValue
        self.attr['CTTH_HEIGHT']['MaskedFill'] = False
        return 
        
    def _CTTH_EFFECT(self):
        self.attr['CTTH_EFFECT']={}
        self.var['CTTH_EFFECT']=np.ma.array(self.ncid.variables[ 'ctth_effectiv'][:])
        self.var['CTTH_EFFECT'].__setmask__(self.mask)
        self.var['CTTH_EFFECT']._sharedmask=False
        self.attr['CTTH_EFFECT']['units']='%'
        self.attr['CTTH_EFFECT']['PALETTE'] = self.ncid.variables['ctth_effectiv_pal'][:]
        self.attr['CTTH_EFFECT']['gain'] = self.ncid.variables['ctth_effectiv'].scale_factor
        self.attr['CTTH_EFFECT']['intercept'] = self.ncid.variables['ctth_effectiv'].add_offset
        self.attr['CTTH_EFFECT']['FillValue'] = self.ncid.variables['ctth_effectiv']._Fillvalue
        return
