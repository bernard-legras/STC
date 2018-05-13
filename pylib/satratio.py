#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mass saturation ratio of water vapour with respect to ice
Using formula from K. Emanuel book

Created on Sat Mar 31 19:34:58 2018

@author: Bernard Legras
"""
import numpy as np
# Calculation of the saturation mixing ratio from actual temperature and pressure
def satratio(p,T):
    """ Calculate the mass saturation ratio from pressure (in Pa) and temperature 
    (in K). Output in kg/kg """
    estar = 1.0008*np.exp(23.33086-(6111.72784/T)+0.15215*np.log(T))
    satr = 0.622 * estar/(0.01*p-estar)
    return satr