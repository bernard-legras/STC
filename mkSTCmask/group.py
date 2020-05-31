#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Common definition of the groups of regions

Created on Thu Oct 24 12:28:45 2019

@author: Bernard Legras
"""

group = {}
group['Land'] = ['IndianSub','SouthChina','Pen','Pakistan','Bangladesh']
group['AsiaLand']  = group['Land']+['NorthChina','JapanKorea']
group['Seas'] = ['BoB','SCSPhi','Philippines']
group['Ocean'] = group['Seas'] + ['IndianOcean','IndoMalaysia','WestPacific','MidPacific']
group['Tibet'] = ['TibetanPlateau',]
group['AfricaArabia'] = ['CentralAfrica','GuineaGulf','NorthAfrica','RedSea','MiddleEast']
group['Asia'] = group['Land']+group['Ocean']+['TibetanPlateau']
group['BigAsia'] = group['Ocean']+group['AsiaLand']+['NorthPacific','NorthAsia','WestAsia','Caspian'] + ['TibetanPlateau',]
group['All'] = group['Asia'] + group['AfricaArabia'] + ['Europe','Atlantic','Mediterranea'] 
