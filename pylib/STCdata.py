# -*- coding: utf-8 -*-
"""

AMES1001 reader for StratoClim

Created on Sun Oct 15 10:43:24 2017
Works in python 3

@author: Bernard Legras
"""

from datetime import datetime
from os.path import basename, join
import numpy as np
import socket

if socket.gethostname() == 'Mentat':
    STC_database = 'C:\\cygwin64\\home\\berna\\data\\STC\\STC-M55\\database'
elif 'spirit' in socket.gethostname():
    STC_database ='/home/legras/STC/STC-M55/database'
elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
    STC_database = '/home/legras//STC/STC-M55/database'
elif socket.gethostname() == 'satie':
    STC_database = '/limbo/data/STC/STC-M55/database'
elif socket.gethostname() == ['zappa','couperin','coltrane','puccini','grapelli']:
    STC_database ='/net/grapelli/limbo/data/STC/STC-M55/database'
elif socket.gethostname() == 'gort':
    STC_database = '/dkol/data/STC/STC-M55/database'
else:
     print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')

""" General readAMES1001 reader, can be used beyond StratoClim """

class readAMES1001(object):
    def __init__(self,file):
        try:
            f = open(file)
        except:
            print('cannot open ',file)
            return
        #print(file)
        self.filename = basename(file)
        # lambdas to read and split line into integers of floats
        #intlb=lambda ll,n: list(int(x) for x in ll[:n])
        floatlb = lambda ll,n: list(float(x) for x in ll[:n])
        intl    = lambda n: list(int(x) for x in f.readline().split()[:n])
        floatl  = lambda n: list(float(x) for x in (f.readline()).split()[:n])
        """ Read header """
        # 1st line
        self.nlhead,self.ffi = intl(2)
        if self.ffi != 1001:
            print('incorrect format index',self.ffi)
            return
        # 2nd to 5th line: Author(s), org, sname, mname
        self.author = f.readline()
        self.org = f.readline()
        self.sname = f.readline()
        self.mname = f.readline()
        self.ivol,self.nvol = intl(2)
        # dates
        #patch to skip hagar
        try:
            _=self.inst
            if (self.inst == 'hagar') and (self.version =='1'):
                dd1,dd2 = intl(2)
                y1 = int(dd1/10000)
                y2 = int(dd2/10000)
                m1 = int(dd1/100-100*y1)
                m2 = int(dd2/100-100*y2)
                d1 = int(dd1-10000*y1-100*m1)
                d2 = int(dd2-10000*y2-100*m2)
            else:
                y1,m1,d1,y2,m2,d2 = intl(6)
        except:
            y1,m1,d1,y2,m2,d2 = intl(6)
        self.date = datetime(y1,m1,d1)
        self.rdate = datetime(y2,m2,d2)
        # dx1
        self.dx1 = floatl(1)[0]
        self.full_xname = f.readline()
        self.nv = intl(1)[0]
        # patches for erroneous files
        #try:
        #    if self.inst == 'ccp-cdp':
        #        self.nv = 7
        #except: pass
        self.vscale = floatl(self.nv)
        self.vmiss = floatl(self.nv)
        self.full_vname=[]
        for i in range(self.nv):
            self.full_vname.append(f.readline())
        self.nscoml = intl(1)[0]
        self.scoml=[]
        for i in range(self.nscoml):
            self.scoml.append(f.readline())
        self.nncoml = intl(1)[0]
        # fix for flash which has wrong nncoml
        # check that the total number of heading lines complies to the sum
        suml = 14 + self.nv + self.nscoml + self.nncoml
        if self.nlhead != suml:
            print('this file heading does not comply to rules ',self.nlhead,suml)
            if self.nlhead > suml:
                print('attempt correcting by reading more ncoml lines')
                self.nncoml += self.nlhead-suml
            #else:
                #print('attempt correcting by reading less ncoml lines' )
                #self.nncoml = max(0,self.nncoml-self.nlhead+suml)
        self.ncoml=[]
        for i in range(self.nncoml):
            self.ncoml.append(f.readline())

        # read data (as a single buffer)
        buffer1 = f.readlines()
        buffer2 = np.array(list(floatlb(ll.split(),1+self.nv) for ll in buffer1))
        del buffer1
        self.x = buffer2[:,0]
        self.var={}
        for i in range(self.nv):
            self.var[i] = buffer2[:,i+1]

""" Specialized StratoClim reader
Is the vname format really special to StratoClim
The last part of the format for the string between {} decoded  into long_name """

class STCdata(readAMES1001):

    def __init__(self,file):
        readAMES1001.__init__(self,file)
        self.xname = 'UTC'
        # extract vname and unit from full_vname
        # consider all the cases with and without unit and long_name
        self.vname = []
        self.vunit = []
        self.long_name = []
        for i in range(self.nv):
            # get name and unit
            try:
                v,u = self.full_vname[i].split('(')
                self.vname.append(v.rstrip())
                self.vunit.append(u.split(')')[0])
            # case without unit, separate name from long_name if any
            except ValueError:
                v = self.full_vname[i].split('{')
                self.vname.append(v[0].rstrip())
                self.vunit.append('None')
            # get long_name
            try:
                v,u = self.full_vname[i].split('{')
                self.long_name.append(u.split('}')[0])
            except ValueError:
                self.long_name.append('')

""" Handy class that reads the file according to instrument name and date """

class STCinst(STCdata):

    def __init__(self,inst,date,version='1',fv=False):
        self.STCdates = {'ktm1':datetime(2017,7,27),'ktm2':datetime(2017,7,29),\
                         'ktm3':datetime(2017,7,31),'ktm4':datetime(2017,8,2),\
                         'ktm5':datetime(2017,8,4),'ktm6':datetime(2017,8,6),\
                         'ktm7':datetime(2017,8,8),'ktm8':datetime(2017,8,10),\
                         'klm1':datetime(2016,8,30),'klm2':datetime(2016,9,1),\
                         'klm3':datetime(2016,9,6)}
        # TRACZILLA run start and end in UTC
        self.STCsegmt = {'ktm1':(7601,29100,36700),'ktm2':(13501,11600,25100),\
                         'ktm3':(11656,11625,23280),'ktm4':(13051,31000,44050),\
                         'ktm5':(12986,12425,25410),'ktm6':(11991,27190,39180),\
                         'ktm7':(9991,13850,23840),'ktm8':(12051,31900,43950)}
        # detection of the date
        if date in self.STCdates.keys():
            self.flight = date
            datett = self.STCdates[date]
        else:
            if not isinstance(date,datetime):
                print('date must be flight id or datetime')
                return
            for name, dd in self.STCdates.items():
                if dd == date:
                    self.flight = name
                    datett = date
            try:
                _ = self.flight
            except:
                print('date does not match known flights')
                return

        # Empty dictionary that may be enriched at later stage
        self.instruments = {'ucse':[],'tdc':[],'flash':[],'fish':[],'chiwis':[],\
                            'amica':[],'cip':[],'hagar':[],\
                            'uhsas':[],'copas':[],'ccp-cdp':[],'ccp-cipgs':[],\
                            'pip':[],'fozan':[],'funmass':[],'stratomas':[],\
                            'erica':[],'nixecaps':[],'sid3':[],'was':[],\
                            'cold':[],'hapaco':[],'phips':[],'haloholo':[],
                            'mas':[]}

        """ other instruments:
            mal non described
            gloria in netcdf """

        if inst not in self.instruments.keys():
            print ('unknown instrument')
            return
        self.inst = inst
        if fv:
            version = self.flight
        file = join(STC_database,inst,datett.strftime('%y%m%d_')+version+'_'+inst+'.nas')
        self.version = version
        STCdata.__init__(self,file)

        """ Generate masked data using the given mask value (or nan).
        The data are also rescaled according to vscale value. Done at the end in order not
        to change the missing value.
        The masked data are contained in a new dictionary with entries named after vnames
        instead of a number.
        We assume that time (x) does not need to be masked. """
        self.varmk = {}
        for iv in self.var.keys():
            # first fill nan with the missing value as several files are not consistent
            # in this respect
            self.var[iv][np.isnan(self.var[iv])] = self.vmiss[iv]
            self.varmk[self.vname[iv]] = np.ma.masked_equal(self.var[iv],self.vmiss[iv])
            self.varmk[self.vname[iv]] *= self.vscale[iv]
        # Finally rename the keys to match vnames
        i = 0
        for var in self.vname:
            self.var[var] = self.var.pop(i)
            i += 1
