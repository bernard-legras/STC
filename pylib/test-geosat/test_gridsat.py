# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 19:36:59 2019

@author: berna
"""

from datetime import datetime
from geosat import GridSat

date = datetime(2017,8,11,6)
gs = GridSat(date)
gs._get_IR0()
gs.chart('IR0')
#%%
gs._get_var('IR2')
gs.chart('IR2')
#%%
gs._get_WV67()
gs.chart('WV67',clim=(210,270))
#%%