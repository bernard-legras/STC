# -*- coding: utf-8 -*-
"""
Tile images for the STC 2019 presentation

Created on Mon 13 May 2019

@author: Bernard Legras
"""
#import numpy as np
from PIL import Image
from os.path import join


# tile images from forward and backward calculations

forw_dir = join('.','figs-Box','All-theta','sh-1728','FullAMA')
forw_dirg = join('.','figs-Box','All-theta','sh-1728','global')
back_dir = join('..','STC-back','figs')
forw_dir2 = join('.','figs-Box')


# %% Tile 1 impact low levels
fname1 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-340K-All-h1728-sh.png')
fname3 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-350K-All-h1728-sh.png')
fname5 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-360K-All-h1728-sh.png')
fname7 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-370K-All-h1728-sh.png')
fname2 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-340K-All-h1728-sh.png')
fname4 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-350K-All-h1728-sh.png')
fname6 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-360K-All-h1728-sh.png')
fname8 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-370K-All-h1728-sh.png')
im = []
for chart in [fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8]:
    aa=Image.open(chart)
    im.append(aa)
    #if (chart == fname7) | (chart == fname8) :
    #    im.append(aa.crop([70,80,935,410]))
    #else:
    #    im.append(aa.crop([70,80,935,351]))
[w,h]=im[0].size
[w2,h2]=im[6].size
bb=Image.new('RGB',(2*w,3*h+h2))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.paste(im[6],(0,3*h,w,3*h+h2))
bb.paste(im[7],(w,3*h,2*w,3*h+h2))
bb.save('figs-Box/tile1.png')

#%% Tile 2 impact high levels
fname1 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-380K-All-h1728-sh.png')
fname3 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-390K-All-h1728-sh.png')
fname5 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-400K-All-h1728-sh.png')
fname7 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-420K-All-h1728-sh.png')
fname2 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-380K-All-h1728-sh.png')
fname4 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-390K-All-h1728-sh.png')
fname6 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-400K-All-h1728-sh.png')
fname8 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-420K-All-h1728-sh.png')
im = []
for chart in [fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8]:
    aa=Image.open(chart)
    im.append(aa.crop([0,0,4115,1272]))
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,4*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.paste(im[6],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.save('figs-Box/tile2.png')

#%% Tile 3 global impact
fname1 = join(forw_dirg,'target','chart-EID-FULL-Box-theta-global-target-380K-All-h1728-sh.png')
fname2 = join(forw_dirg,'target','chart-EID-FULL-Box-theta-global-target-400K-All-h1728-sh.png')
fname3 = join(forw_dirg,'source','chart-EID-FULL-Box-theta-global-source-380K-All-h1728-sh.png')
fname4 = join(forw_dirg,'source','chart-EID-FULL-Box-theta-global-source-400K-All-h1728-sh.png')
fname5 = join(forw_dir,'target','chart-EID-FULL-Box-theta-fullAMA-target-380K-All-h1728-sh.png')
fname6 = join(forw_dir,'target','chart-EID-FULL-Box-theta-fullAMA-target-400K-All-h1728-sh.png')
im = []
for chart in [fname1,fname2,fname3,fname4,fname5,fname6]:
    aa=Image.open(chart)
    if (chart==fname2) | (chart==fname1):
        im.append(aa.crop([0,0,3800,1679]))
    else:
        im.append(aa.crop([0,0,4115,1272]))
[w,h]=im[0].size
[w2,h2]=im[3].size
w3 = int((w2-w)/2)
bb=Image.new('RGB',size=(2*w2,h+2*h2),color=(255,255,255))
bb.paste(im[4],(0,0,w2,h2))
bb.paste(im[5],(w2,0,2*w2,h2))
bb.paste(im[0],(w3,h2,w+w3,h+h2))
bb.paste(im[1],(w2+w3,h2,w2+w3+w,h+h2))
bb.paste(im[2],(0,h+h2,w2,h+2*h2))
bb.paste(im[3],(w2,h+h2,2*w2,h+2*h2))

bb.save('figs-Box/tile3.png')


#%% Tile 4 mean age in the target and in the source domain
fname1 = join(forw_dir,'mage-target','chart-EAD-Box-theta-FullAMA-mage-target-360K-All-h1728-sh.png')
fname2 = join(forw_dir,'mage-target','chart-EAD-Box-theta-FullAMA-mage-target-380K-All-h1728-sh.png')
fname3 = join(forw_dir,'mage-source','chart-EAD-Box-theta-FullAMA-mage-source-360K-All-h1728-sh.png')
fname4 = join(forw_dir,'mage-source','chart-EAD-Box-theta-FullAMA-mage-source-380K-All-h1728-sh.png')
im = []
for chart in [fname1,fname2,fname3,fname4]:
    aa=Image.open(chart)
    im.append(aa.crop([0,0,4115,1272]))
[w,h]=im[0].size
bb=Image.new('RGB',size=(2*w,2*h),color=(255,255,255))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.save('figs-Box/tile4.png')

#%% Tile 5 comparison impact / back percentage
fname1 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-360K-All-h1728-sh.png')
fname2 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-370K-All-h1728-sh.png')
fname3 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-380K-All-h1728-sh.png')
fname4 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-400K-All-h1728-sh.png')
fname5 = join(back_dir,'chart-EAD-percentage-hits-360K.png')
fname6 = join(back_dir,'chart-EAD-percentage-hits-370K.png')
fname7 = join(back_dir,'chart-EAD-percentage-hits-380K.png')
fname8 = join(back_dir,'chart-EAD-percentage-hits-400K.png')
im = []
for chart in [fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8]:
    aa=Image.open(chart)
    im.append(aa.crop([0,0,4115,1272]))
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,4*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[4],(w,0,2*w,h))
bb.paste(im[1],(0,h,w,2*h))
bb.paste(im[5],(w,h,2*w,2*h))
bb.paste(im[2],(0,2*h,w,3*h))
bb.paste(im[6],(w,2*h,2*w,3*h))
bb.paste(im[3],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.save('figs-Box/tile5.png')

#%% Tile 6 comparison sources
fname1 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-360K-All-h1728-sh.png')
fname2 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-370K-All-h1728-sh.png')
fname3 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-380K-All-h1728-sh.png')
fname4 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-400K-All-h1728-sh.png')
fname5 = join(back_dir,'chart-EAD-distrib-sources-360K.png')
fname6 = join(back_dir,'chart-EAD-distrib-sources-370K.png')
fname7 = join(back_dir,'chart-EAD-distrib-sources-380K.png')
fname8 = join(back_dir,'chart-EAD-distrib-sources-400K.png')
im = []
for chart in [fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8]:
    aa=Image.open(chart)
    im.append(aa.crop([0,0,4115,1272]))
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,4*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[4],(w,0,2*w,h))
bb.paste(im[1],(0,h,w,2*h))
bb.paste(im[5],(w,h,2*w,2*h))
bb.paste(im[2],(0,2*h,w,3*h))
bb.paste(im[6],(w,2*h,2*w,3*h))
bb.paste(im[3],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.save('figs-Box/tile6.png')

#%% Tile 7 comparison confinement
fname1 = join(forw_dir2,'target-target-theta','chart-EAD-Box-theta-FullAMA-Jul-11-target-370K-sh.png')
fname2 = join(forw_dir2,'target-target-theta','chart-EAD-Box-theta-FullAMA-Jul-21-target-370K-sh.png')
fname3 = join(forw_dir2,'target-target-theta','chart-EAD-Box-theta-FullAMA-Aug-01-target-370K-sh.png')
fname4 = join(forw_dir2,'target-target-theta','chart-EAD-Box-theta-FullAMA-Aug-11-target-370K-sh.png')
fname5 = join(forw_dir2,'target-target-theta','chart-EAD-Box-theta-FullAMA-Aug-21-target-370K-sh.png')
fname6 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-370K-All-h1728-sh.png')
im = []
for chart in [fname1,fname2,fname3,fname4,fname5,fname6]:
    aa=Image.open(chart)
    im.append(aa.crop([0,0,4115,1272]))
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,3*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.save('figs-Box/tile7.png')

