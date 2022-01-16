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
z1=lambda p: 44.3308*(1-0.267993*p**(0.190263))
z2=lambda p: 11. + 6.34162 * np.log(226.32/p)
z3=lambda p: 20. + 216.650 *(-1 + 1.12431*p**(-0.0292712))
z4=lambda p: 32. + 81.6607 *(-1 + 1.19377*p**(-0.0819595))
z5=lambda p: 47. + 7.92226 * np.log(1.10906/p)

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