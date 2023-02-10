#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mass saturation ratio of water vapour with respect to ice
Using formula from K. Emanuel book

Created on Sat Mar 31 19:34:58 2018
Liquid formula added on 12/12/2022

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
def liquid_satratio(p,T):
    """ Calculate the mass saturation ratio from pressure (in Pa) and temperature 
    (in K). Output in kg/kg """
    estar = 6.112 * np.exp(17.67*(T-273.15)/(T-29.65))
    satr = 0.622 * estar/(0.01*p-estar)
    return satr
def liquid_satratio_2(p,T):
    """ Alternate formula from Vaissala """
    a = -6096.9385
    b = 21.2409642
    c = -2.711193e-2
    d = 1.673952e-5
    e = 2.433502
    vs = np.exp( (a/T) + b + (c*T) + (d*T**2) + e*np.log(T))
    satr = 0.622 * vs/(p-vs)
    return satr
