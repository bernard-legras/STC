# -*- coding: utf-8 -*-
"""
Test of SAFNWCnc for 2021 2022 and the v2018.1-HVR version

Created on 25/01/2023

@author: Bernard Legras
"""

from datetime import datetime
import SAFNWCnc
import geosat
import matplotlib.pyplot as plt

#date = datetime(2022,8,31,15)
date = datetime(2021,8,31,15)
version = 'v2018.1-HVR'

sat1 = SAFNWCnc.SAFNWC_CTTH(date,'msg1',BBname='SAFBox',version=version)

sat2 = SAFNWCnc.SAFNWC_CTTH(date,'msg1')

sat1._CTTH_PRESS()
sat2._CTTH_PRESS()
sat1.show('CTTH_PRESS',clim=(50,80000))
sat2.show('CTTH_PRESS',clim=(50,80000))
sat1.var['CTTH_PRESS_REF'] = sat2.var['CTTH_PRESS'][341:1858,307:3352]
sat1.show('CTTH_PRESS_REF',clim=(50,80000))
pp1 = sat1.var['CTTH_PRESS'][:].compressed()
pp2 = sat1.var['CTTH_PRESS_REF'][:].compressed()

plt.hist(pp1[pp1<1000],bins=100)
plt.show()
plt.hist(pp2[pp2<1000],bins=100)
plt.show()
plt.hist(pp1[pp1<200],bins=100)
plt.show()
plt.hist(pp2[pp2<200],bins=100)
plt.show()
