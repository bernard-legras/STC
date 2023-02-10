#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 14:54:26 2018
Physical constants and reference values

@author: Bernard Legras
"""
# Standard thermodynamical constants fot the atmosphere
Ra = 8.3144626181
R = 287.048
Rv = 461.500
Cp = 1005.7
Cpv = 1878.
kappa = R/Cp
g = 9.81
p0 = 100000.
pref = 101325.
M = 28.9647
REarth = 6372000.
Omega = 2*3.141592653589793/86164
Zero_Celsius = 273.15
Na = 6.02214076e23
epsilon = R/Rv
# Formula from the Wikipedia latent heat page
Ll = lambda T: 1000*(2500.8 - 2.36 * (T-273.15) + 0.016 * (T-273.15)**2 \
                    - 0.00006 *(T-273.15)**3)
Li = lambda T: 1000*(2834.1 - 0.29 * (T-273.15) - 0.004 *(T-273.15)**2)
Cl = 4190
Ci = 2106
L0 = Ll(Zero_Celsius)
