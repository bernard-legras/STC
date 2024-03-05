#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 11:12:41 2017

Official StratoClim color map

@author: Bernard Legras
"""
import matplotlib.colors as colors
# color list with 10 colors
#listcolors=['#253494','#2c7fb8','#41b6c4','#a1dab4','#ffffcc','#fef0d9',
#            '#fdcc8a','#fc8d59','#e34a33','#b30000']
# color list with 20 colors  
listcolors20 = ['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']            
mymap20 = colors.ListedColormap(listcolors20)
listcolors16 = ['#161d58','#253494','#2850a6','#2c7fb8','#41b6c4',
            '#71c8bc','#a1dab4','#ffffcc','#fef0d9',
            '#fdcc8a','#fdac72','#fc8d59','#e34a33','#cb251a',
            '#b30000','#7f0000']            
mymap16 = colors.ListedColormap(listcolors16)
mymap = mymap20
