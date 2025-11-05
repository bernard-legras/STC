import numpy as np
import geosat
import os
import re
import glob
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
            nam='HIMA0?'
            region='globeJ'
            geosat.read_mask_himawari()
            masksat=geosat.mask_sat['himawari']
            VISIR = 'NR'
        elif sat=='msg0':
            nam='MSG0'
            region='globeM'
            geosat.read_mask_MSG()        
            masksat=geosat.mask_sat['msg']
            VISIR = 'VISIR'
        elif sat=='msg3':
            nam='MSG3'
            region='globeM'
            geosat.read_mask_MSG()        
            masksat=geosat.mask_sat['msg']
            VISIR = 'VISIR'    
        elif sat=='msg4':
            nam='MSG4'
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
        elif sat=='msg2':
            nam='MSG2'
            region='globeI'
            VISIR = 'VISIR'
            geosat.read_mask_MSG()
            masksat=geosat.mask_sat['msg']
        elif sat=='goesw':
            nam='GOES17'
            region='globeW'
            VISIR = 'NR'
            geosat.read_mask_GOES()
            masksat=geosat.mask_sat['goes']
        elif sat=='goese':
            nam='GOES16'
            region='globeE'
            VISIR = 'NR'
            geosat.read_mask_GOES()
            masksat=geosat.mask_sat['goes']
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
        temp_fullname = os.path.join(root_dir,sat,safnwc,date.strftime("%Y"),
                            date.strftime("%Y_%m_%d"),filename)
        try:
            fullname = glob.glob(temp_fullname)[0]
            self.ncid = Dataset(fullname, mode='r')
        except IndexError:
            print('NOT FOUND ',temp_fullname)
            self.ncid = None
            return
        except AttributeError:
            print('Cannot open',filename)
            self.ncid = None
            return
        except:
            print('Other error',filename)
            self.ncid = None
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
        self.attr['CMa']['flag_values'] = (0, 1)
        self.attr['CMa']['flag_meaning'] = ('Cloud_free', 'Cloudy')      
        return
        
    def _CMa_QUALITY(self):
        '''
        returns numpy array containing the quality information on the CMa product from SAF NWC file
        '''
        self.attr['CMa_QUALITY']={}
        self.var['CMa_QUALITY']=np.ma.array(self.ncid.variables['cma_quality'][:])
        self.var['CMa_QUALITY'].__setmask__(self.mask)
        self.var['CMa_QUALITY']._sharedmask=False
        self.attr['CMa_QUALITY']['flag_mask']   = (1, 2, 4, 56, 56, 56, 56)
        self.attr['CMa_QUALITY']['flag_values'] = (1, 2, 4,  8, 16, 24, 32)
        self.attr['CMa_QUALITY']['flag_meanings']=('nodata', 'internal_consistency', 'temporal_consistency',\
                                 'good', 'questionable',' bad',' interpolated')
        return
        
    def _CMa_STATUS(self):
        '''
        returns numpy array containing the status information on the CMa product from SAF NWC file
        '''
        self.attr['CMa_STATUS']={}
        self.var['CMa_STATUS']=np.ma.array(self.ncid.variables['cma_status_flag'][:])
        self.var['CMa_STATUS'].__setmask__(self.mask)
        self.var['CMa_STATUS']._sharedmask=False
        self.attr['CMa_STATUS']['flag_values'] = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
        self.attr['CMa_STATUS']['flag_meanings'] = ('Low_level_thermal_inversion_in_NWP_field',\
                                                  'Cold_snowy_ground_suspected',\
                                                  'Temporal_algorithm_passed',\
                                                  'High_resolution_satellite_data_used',\
                                                  'RTTOV_on-line_used', 'SST_analysis_available',\
                                                  'Snow_map_available', 'Sea_ice_map_available',\
                                                  'No_method_for_dust',\
                                                  ' No_method_for_volcanic_plume', 'No_method_for_smoke')
        return
                
    def _CMa_DUST(self):
        self.attr['CMa_DUST']={}
        self.var['CMa_DUST']=np.ma.array(self.ncid.variables['cma_dust'][:])
        self.var['CMa_DUST'].__setmask__(self.mask)
        self.var['CMa_DUST']._sharedmask=False
        self.attr['CMa_DUST']['flag_values'] = (0,1,2)
        self.attr['CMa_DUST']['flag_meanings'] = ('No_dust', 'Dust', 'Undefined_separability_problems')
        return 
    def _CMa_VOLCANIC(self):
        self.attr['CMa_VOLCANIC']={}
        self.var['CMa_VOLCANIC'] = np.ma.array(self.ncid.variables['cma_volcanic'][:])
        self.var['CMa_VOLCANIC'].__setmask__(self.mask)
        self.var['CMa_VOLCANIC']._sharedmask=False
        self.attr['CMa_VOLCANIC']['flag_values'] = (0, 1, 2)
        self.attr['CMa_VOLCANIC']['flag_meanings'] = ('No_volcanic_plume', 'Volcanic_plume', 'Undefined_separability_problems')
        return 
    def _CMa_SMOKE(self):
        self.attr['CMa_SMOKE']={}
        self.var['CMa_SMOKE']=np.ma.array(self.ncid.variables['cma_smoke'][:])
        self.var['CMa_SMOKE'].__setmask__(self.mask)
        self.var['CMa_SMOKE']._sharedmask=False
        self.attr['CMa_SMOKE']['flag_values'] = (0, 1, 2)
        self.attr['CMa_SMOKE']['flag_meanings'] = ('No_smoke', 'Smoke', 'Undefined_separability_problems')
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
        self.attr['CT']['units']=('1: Cloud free land','2: Cloud free sea','3: Snow over land','4: Snow over sea',\
                                  '5: Very low clouds','6: Low Clouds','7: Mid-level clouds','8: High opaque clouds',\
                                  '9: Very high opaque clouds','10: Fractional clouds','11: High semitransparent thin clouds',\
                                  '12: High semitransparent meanly thick clouds','13: High semitransparent thick clouds',\
                                  '14: High semitransparent above low or medium clouds','15: High semitransparent above snow/ice')
        self.attr['flag_values'] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
        self.attr['flag_meanings'] = ('Cloud-free_land', 'Cloud-free_sea', 'Snow_over_land', 'Sea_ice',\
                                      'Very_low_clouds', 'Low_clouds', 'Mid-level_clouds', 'High_opaque_clouds',\
                                      'Very_high_opaque_clouds', 'Fractional_clouds',\
                                      'High_semitransparent_thin_clouds', 'High_semitransparent_moderately_thick_clouds',\
                                      'High_semitransparent_thick_clouds',  'High_semitransparent_above_low_or_medium_clouds',\
                                      'High_semitransparent_above_snow_ice')
        self.attr['CT']['PALETTE'] = self.ncid.variables['ct_pal'][:]
        return

    def _CT_QUALITY(self):
        '''
        returns numpy array containing the quality information on the CT product from SAF NWC file
        '''
        self.attr['CT_QUALITY']={}
        self.var['CT_QUALITY']=np.ma.array(self.ncid.variables['ct_quality'][:])
        self.var['CT_QUALITY'].__setmask__(self.mask)
        self.var['CT_QUALITY']._sharedmask=False
        self.attr['CT_QUALITY']['flag_mask']   = (1, 2, 4, 56, 56, 56, 56)
        self.attr['CT_QUALITY']['flag_values'] = (1, 2, 4,  8, 16, 24, 32)
        self.attr['CT_QUALITY']['flag_meanings']=('nodata','internal_consistency','temporal_consistency',\
                                 'good','questionable','bad','interpolated')
        return
        
    def _CT_STATUS(self):
        '''
        returns numpy array containing the status information on the CT product from SAF NWC file
        '''
        self.attr['CT_STATUS']={}
        self.var['CT_STATUS']=np.ma.array(self.ncid.variables['ct_status_flag'][:])
        self.var['CT_STATUS'].__setmask__(self.mask)
        self.var['CT_STATUS']._sharedmask=False
        self.attr['CT_STATUS']['flag_values']=(1, 2, 4, 8, 16, 32)
        self.attr['CT_STATUS']['flag_meanings']=('Low_level_thermal_inversion_in_NWP_field',\
                                                 'Tropopause_temperature_available_from_NWP',\
                                                 '138um_used_for_cirrus_identification',\
                                                 'High_resolution_satellite_data_used',\
                                                 'No_method_for_stratiform_cumuliform_separation',\
                                                 'No_method_for_multilayer')
        return
     
    def _CT_CONDITIONS(self):
        self.attr['CT_CONDITIONS']={}
        self.var['CT_CONDITIONS']=np.ma.array(self.ncid.variables['ct_conditions'][:])
        self.var['CT_CONDITIONS'].__setmask__(self.mask)
        self.var['CT_CONDITIONS']._sharedmask=False
        self.attr['CT_CONDITIONS']['flag_mask']  = (1, 6, 6, 6, 8, 48, 48, 48, 64, 128, 768, 768, 768, 3072, 3072, 3072, 
                                               12288, 12288, 12288, 49152, 49152, 49152)
        self.attr['CT_CONDITIONS']['flag_values']= (1, 2, 4, 6, 8, 16, 32, 48, 64, 128, 256, 512, 768, 1024, 2048, 3072,
                                                4096,  8192, 12288, 16384, 32768, 49152)
        self.attr['CT_CONDITIONS']['flag_meanings'] = ('space', 'night', 'day', 'twilight', 'sunglint', 'land', 'sea', 'coast',\
                                   'not_used', 'not_used',\
                                   'all_satellite_channels_available', 'useful_satellite_channels_missing', 'mandatory_satellite_channels_missing',\
                                   'all_NWP_fields_available', 'useful_NWP_fields_missing', 'mandatory_NWP_fields_missing',\
                                   'all_product_data_available', 'useful_product_data_missing', 'mandatory_product_data_missing',\
                                   'all_auxiliary_data_available', 'useful_auxiliary_data_missing', 'mandatory_auxiliary_data_missing')
        return
        
    def _CT_CUMULIFORM(self):
        self.attr['CT_CUMULIFORM']={}
        self.var['CT_CUMULIFORM']=np.ma.array(self.ncid.variables['ct_cumuliform'][:])
        self.var['CT_CUMULIFORM'].__setmask__(self.mask)
        self.var['CT_CUMULIFORM']._sharedmask=False
        self.attr['CT_CUMULIFORM']['flag_values']= (1, 2, 3, 4, 5)
        self.attr['CT_CUMULIFORM']['flag_meanings'] = ('Stratiform_status', 'Cumuliform_status', 'Mixed_status', 'Cloud-free',\
                                   'Undefined_separability_problems')
        return
        
    def _CT_MULTILAYER(self):
        self.attr['CT_MULTILAYER']={}
        self.var['CT_MULTILAYER']=np.ma.array(self.ncid.variables['ct_multilayer'][:])
        self.var['CT_MULTILAYER'].__setmask__(self.mask)
        self.var['CT_MULTILAYER']._sharedmask=False
        self.attr['CT_MULTILAYER']['flag_values']= (0, 1, 2, 3)
        self.attr['CT_MULTILAYER']['flag_meanings'] = ('No_multilayer_detected', 'Multilayer_detected',\
                                                       'Cloud_free', 'Undefined_separability_problems')
        return

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
