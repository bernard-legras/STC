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
# Calculation of the saturation mass mixing ratio from actual temperature and pressure
def satratio(p,T):
    """ Calculate the mass saturation ratio from pressure (in Pa) and temperature
    (in K). Output in kg/kg.Given by Emanuel to be valid up to 0.3% over the range 190K-273K
    """
    estar = 1.0008*np.exp(23.33086-(6111.72784/T)+0.15215*np.log(T))
    satr = 1.0021 * 0.622 * estar/(0.01*p-estar)
    return satr, estar*100
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
def esati_murphy(p, T):
    """ei in Pa saturation vapor pressure with respect to hexagonal (most stable) ice
       Murphy and Koop 2005, QJRMS
       Input in Pa and K
       Output in kg/kg """
    lnP = 9.550426-5723.265/T+3.53068*np.log(T)-0.00728332*T
    estar = np.exp(lnP)
    esati_murphy = estar / (p-estar)
    #: Convert to kg/kg
    esati_murphy = esati_murphy * 0.622
    #ppmv
    #kg
    return esati_murphy, estar

#%% Tests in the main program
""" The values of the tests are taken from the table C1 of Murphy & Koop, 2005,
    QJRMS. """
if __name__ == '__main__':
    Tl = (150, 180, 210, 240, 273.15, 273.16)
    pl = 100
    for T in Tl:
        rvs1,estar = esati_murphy(pl, T)
        print('MK {:.2f} {:.4E}'.format(T, estar))
        rvs2,estar = satratio(pl,T)
        print('KE {:.2f} {:.4E}'.format(T, estar))
