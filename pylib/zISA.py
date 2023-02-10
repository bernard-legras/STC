#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 17:56:12 2018

Barometric altitude following the standard atmospheric profile

hbaro uses a simplified form which is usable up to 11 km

zISA can be used up to the stratopause

The barometric altitude can be signicantly smaller than the true altitude in the tropics

@author: Bernard Legras
"""

import numpy as np

# barometric altitude in km from p in Pa
hbaro=lambda p: (1-(p/101325.)**0.1903)*(0.28825/0.0065)

# ISA altitude from p in hPa
# this is probably the geopotential altitude
z1 = lambda p: 44.3308*(1-0.267993*p**(0.190263))
z2 = lambda p: 11. + 6.34162 * np.log(226.32/p)
z3 = lambda p: 20. + 216.650 *(-1 + 1.12431*p**(-0.0292712))
z4 = lambda p: 32. + 81.6607 *(-1 + 1.19377*p**(-0.0819595))
z5 = lambda p: 47. + 7.92226 * np.log(1.10906/p)

def zISA(p):
    if p>22632.:
        return z1(0.01*p)
    elif p>5474.88:
        return z2(0.01*p)
    elif p>868.016:
        return z3(0.01*p)
    elif p>110.906:
        return z4(0.01*p)
    else:
        return z5(0.01*p)
zISA=np.vectorize(zISA)

# Derivation of the pressure from the geopotential altitude, reverting the formulas z1-z5
p1 = lambda z: ((1-z/44.3308)/0.267993)**(1/0.190263)
p2 = lambda z: 226.32 * np.exp(-(z-11.)/6.34162)
p3 = lambda z: ((1+(z-20.)/216.650)/1.12431)**(-1./0.0292712)
p4 = lambda z: ((1+(z-32.)/81.6607)/1.19377)**(-1./0.0819595)
p5 = lambda z: 1.10906 * np.exp(-(z-47.)/7.92226)

def pISA(z):
    if z < 11.:
        return 100 * p1(z)
    elif z < 20.:
        return 100 * p2(z)
    elif z < 32.:
        return 100 * p3(z)
    elif z < 47.:
        return 100 * p4(z)
    else:
        return 100 * p5(z)
pISA = np.vectorize(pISA)

T0 = 273.15
T1 = lambda z: T0 + 15.0265- 6.5 * z
T2 = lambda z: T0 - 56.5
T3 = lambda z: T0 - 56.5 + (z-20.)
T4 = lambda z: T0 - 44.5 + 2.8 * (z-32.)
T5 = lambda z: T0 - 2.5
T6 = lambda z: T0 - 2.5 - 2.8 * (z-51.)
T7 = lambda z: T0 - 58.5 - 2. * (z-71.)
T8 = lambda z: T0 - 86.28

# ISA temperature as a function of geopotentiel altitude in km

def TISA(z):
    if z < 11.:
        return T1(z)
    elif z < 20.:
        return T2(z)
    elif z < 32.:
        return T3(z)
    elif z < 47.:
        return T4(z)
    elif z < 51.:
        return T5(z)
    elif z < 71.:
        return T6(z)
    elif z < 84.852:
        return T7(z)
    else:
        return T8(z)
        return T
TISA = np.vectorize(TISA)
         