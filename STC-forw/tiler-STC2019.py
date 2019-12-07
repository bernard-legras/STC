# -*- coding: utf-8 -*-
"""
Tile images for the STC 2019 presentation

Created on Mon 13 May 2019

@author: Bernard Legras
"""
# Initialization sequence to be run as a cell
#import numpy as np
from PIL import Image, ImageDraw, ImageFont
from os.path import join
from os import chdir

try:
    chdir('data')
    chdir('STC')
    chdir('STC-forw')
except:pass
# tile images from forward and backward calculations

forw_dir = join('.','figs-Box','All-theta','sh-1728','FullAMA')
forw_dirg = join('.','figs-Box','All-theta','sh-1728','global')
back_dir = join('..','STC-back','figs')
forw_dir2 = join('.','figs-Box')
forw_target = join('.','figs-Box','target-target-theta')


# %% Tile 1 impact low levels
fname1 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-340K-All-h1728-sh.png')
fname3 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-350K-All-h1728-sh.png')
fname5 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-360K-All-h1728-sh.png')
fname7 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-370K-All-h1728-sh.png')
fname2 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-340K-All-h1728-sh.png')
fname4 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-350K-All-h1728-sh.png')
fname6 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-360K-All-h1728-sh.png')
fname8 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-370K-All-h1728-sh.png')
im = []
for chart in [fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8]:
    aa=Image.open(chart)
    cc = Image.new('RGB',size=(3610,1132),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
    #im.append(aa)
    #if (chart == fname7) | (chart == fname8) :
    #    im.append(aa.crop([70,80,935,410]))
    #else:
    #    im.append(aa.crop([70,80,935,351]))
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

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((300,150),u'target 340 K',fill='yellow',font=font)
draw.text((300,h+150),u'target 350 K',fill='yellow',font=font)
draw.text((300,2*h+150),u'target 360 K',fill='yellow',font=font)
draw.text((300,3*h+150),u'target 370 K',fill='yellow',font=font)
draw.text((w+300,150),u'source 340 K',fill='yellow',font=font)
draw.text((w+300,h+150),u'source 350 K',fill='yellow',font=font)
draw.text((w+300,2*h+150),u'source 360 K',fill='yellow',font=font)
draw.text((w+300,3*h+150),u'source 370 K',fill='yellow',font=font)
bb.save('figs-Box/tile1.png',dpi=(300,300))
bb.show()

#%% Tile 2 impact high levels
fname0 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-380K-All-h1728-sh.png')
fname2 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-390K-All-h1728-sh.png')
fname4 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-400K-All-h1728-sh.png')
fname6 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-420K-All-h1728-sh.png')
fname1 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-380K-All-h1728-sh.png')
fname3 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-390K-All-h1728-sh.png')
fname5 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-400K-All-h1728-sh.png')
fname7 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-420K-All-h1728-sh.png')
im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7]:
    aa=Image.open(chart)
    if chart == fname5:
        im.append(aa.crop([0,0,3610,1132]))
    else:
        cc = Image.new('RGB',size=(3610,1132),color=(255,255,255))
        cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
        im.append(cc)
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

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((300,150),u'target 380 K',fill='yellow',font=font)
draw.text((300,h+150),u'target 390 K',fill='yellow',font=font)
draw.text((300,2*h+150),u'target 400 K',fill='yellow',font=font)
draw.text((300,3*h+150),u'target 420 K',fill='yellow',font=font)
draw.text((w+300,150),u'source 380 K',fill='yellow',font=font)
draw.text((w+300,h+150),u'source 390 K',fill='yellow',font=font)
draw.text((w+300,2*h+150),u'source 400 K',fill='yellow',font=font)
draw.text((w+300,3*h+150),u'source 420 K',fill='yellow',font=font)

bb.save('figs-Box/tile2.png',dpi=(300,300))
bb.show()

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

#%% Tile 3b global impact
fname0 = join(forw_dir2,'SPE','chart-EID-FULL-Box-theta-global-target-SPE-380K-All-h1728-sh.png')
fname1 = join(forw_dir2,'SPE','chart-EID-FULL-Box-theta-FullAMA-target-SPE-380K-All-h1728-sh.png')
fname2 = join(forw_dir2,'SPE','chart-EID-FULL-Box-theta-global-source-SPE-380K-All-h1728-sh.png')
fname3 = join(forw_dir2,'SPE','chart-EID-FULL-Box-theta-FullAMA-source-SPE-380K-All-h1728-sh.png')
fname4 = join(forw_dir2,'SPE','chart-EAD-Box-theta-FullAMA-target-SPE-380K-All-h1728-sh.png')
fname5 = join(forw_dir2,'SPE','chart-EAD-Box-theta-FullAMA-source-SPE-380K-All-h1728-sh.png')
im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5]:
    aa=Image.open(chart)
   #im.append(aa.crop([0,0,4115,1679]))
    if chart==fname0 :
        hh = 1673
    else:
        hh = 1141
    cc = Image.new('RGB',size=(3644,hh),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
[w4,h4]=im[5].size
[w2,h2]=im[2].size
[w0,h0]=im[0].size
wi = int((w2-w0)/2)
hi = int((h0-h2)/2)
bb=Image.new('RGB',size=(2*w2,h2+h0+h4),color=(255,255,255))
bb.paste(im[1],(0,0,w2,h2))
bb.paste(im[3],(w2,0,2*w2,h2))
bb.paste(im[0],(wi,h2,w0+wi,h0+h2))
bb.paste(im[2],(w2,h2+hi,2*w2,2*h2+hi))
bb.paste(im[4],(0,h0+h2,w2,h0+h2+h4))
bb.paste(im[5],(w2,h0+h2,2*w2,h0+h2+h4))

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h2),u'(c)',fill='blue',font=font)
draw.text((0,h2+h0),u'(e)',fill='blue',font=font)
draw.text((w2,0),u'(b)',fill='blue',font=font)
draw.text((w2,h2+hi),u'(d)',fill='blue',font=font)
draw.text((w2,h2+h0),u'(f)',fill='blue',font=font)

bb.save('figs-Box/tile3b.png',dpi=(300,300))
bb.show()

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

##%% Tile 5 comparison impact / back percentage
#fname1 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-360K-All-h1728-sh.png')
#fname2 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-370K-All-h1728-sh.png')
#fname3 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-380K-All-h1728-sh.png')
#fname4 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-400K-All-h1728-sh.png')
#fname5 = join(back_dir,'chart-EAD-percentage-hits-Jul-Aug-Sep-360K.png')
#fname6 = join(back_dir,'chart-EAD-percentage-hits-Jul-Aug-Sep-370K.png')
#fname7 = join(back_dir,'chart-EAD-percentage-hits-Jul-Aug-Sep-380K.png')
#fname8 = join(back_dir,'chart-EAD-percentage-hits-Jul-Aug-Sep-400K.png')
#im = []
#for chart in [fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8]:
#    aa=Image.open(chart)
#    #im.append(aa.crop([0,0,4115,1272]))
#    im.append(aa)
#[w,h]=im[0].size
#bb=Image.new('RGB',(2*w,4*h))
#bb.paste(im[0],(0,0,w,h))
#bb.paste(im[4],(w,0,2*w,h))
#bb.paste(im[1],(0,h,w,2*h))
#bb.paste(im[5],(w,h,2*w,2*h))
#bb.paste(im[2],(0,2*h,w,3*h))
#bb.paste(im[6],(w,2*h,2*w,3*h))
#bb.paste(im[3],(0,3*h,w,4*h))
#bb.paste(im[7],(w,3*h,2*w,4*h))
#bb.save('figs-Box/tile5.png')

##%% Tile 6 comparison sources
#fname1 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-360K-All-h1728-sh.png')
#fname2 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-370K-All-h1728-sh.png')
#fname3 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-380K-All-h1728-sh.png')
#fname4 = join(forw_dir,'source','chart-EAD-Box-theta-FullAMA-source-400K-All-h1728-sh.png')
#fname5 = join(back_dir,'chart-EAD-distrib-sources-360K.png')
#fname6 = join(back_dir,'chart-EAD-distrib-sources-370K.png')
#fname7 = join(back_dir,'chart-EAD-distrib-sources-380K.png')
#fname8 = join(back_dir,'chart-EAD-distrib-sources-400K.png')
#im = []
#for chart in [fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8]:
#    aa=Image.open(chart)
#    im.append(aa.crop([0,0,4115,1272]))
#[w,h]=im[0].size
#bb=Image.new('RGB',(2*w,4*h))
#bb.paste(im[0],(0,0,w,h))
#bb.paste(im[4],(w,0,2*w,h))
#bb.paste(im[1],(0,h,w,2*h))
#bb.paste(im[5],(w,h,2*w,2*h))
#bb.paste(im[2],(0,2*h,w,3*h))
#bb.paste(im[6],(w,2*h,2*w,3*h))
#bb.paste(im[3],(0,3*h,w,4*h))
#bb.paste(im[7],(w,3*h,2*w,4*h))
#bb.save('figs-Box/tile6.png')

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

#%% Tile 8 tiling of images for the backward calculations in the same way as 
# tile 2 and 3 for the forward calculations plus the pdf of the sources

fname0 = join(back_dir,'EAD-JAS','chart-EAD-percentage-hits-Jul-Aug-Sep-360K.png')
fname1 = join(back_dir,'EAD-JAS','chart-EAD-percentage-hits-Jul-Aug-Sep-370K.png')
fname2 = join(back_dir,'EAD-JAS','chart-EAD-percentage-hits-Jul-Aug-Sep-380K.png')
fname3 = join(back_dir,'EAD-JAS','chart-EAD-percentage-hits-Jul-Aug-Sep-400K.png')
fname4 = join(back_dir,'EAD-JAS','chart-EAD-distrib-sources-Jul-Aug-Sep-360K.png')
fname5 = join(back_dir,'EAD-JAS','chart-EAD-distrib-sources-Jul-Aug-Sep-370K.png')
fname6 = join(back_dir,'EAD-JAS','chart-EAD-distrib-sources-Jul-Aug-Sep-380K.png')
fname7 = join(back_dir,'EAD-JAS','chart-EAD-distrib-sources-Jul-Aug-Sep-400K.png')
fname8 = join(back_dir,'EAD-JAS','plot-source-theta-EAD-Jul-Aug-Sep-360K.png')
fname9 = join(back_dir,'EAD-JAS','plot-source-theta-EAD-Jul-Aug-Sep-370K.png')
fname10 = join(back_dir,'EAD-JAS','plot-source-theta-EAD-Jul-Aug-Sep-380K.png')
fname11 = join(back_dir,'EAD-JAS','plot-source-theta-EAD-Jul-Aug-Sep-400K.png')

im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,fname9,fname10,fname11]:
    aa=Image.open(chart)
    if chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7]:
        print(chart)
        cc = Image.new('RGB',size=(3644,1173),color=(255,255,255))
        cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
        im.append(cc)
    elif chart in (fname8,fname9,fname10,fname11):
        im.append(aa.crop((0,0,1219,1173)))
    else:
        im.append(aa)
       
[w,h] = im[0].size
[w2,h2] = im[8].size
bb=Image.new('RGB',(2*w+w2,4*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[4],(w,0,2*w,h))
bb.paste(im[1],(0,h,w,2*h))
bb.paste(im[5],(w,h,2*w,2*h))
bb.paste(im[2],(0,2*h,w,3*h))
bb.paste(im[6],(w,2*h,2*w,3*h))
bb.paste(im[3],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.paste(im[8],(2*w,0,2*w+w2,h))
bb.paste(im[9],(2*w,h,2*w+w2,2*h))
bb.paste(im[10],(2*w,2*h,2*w+w2,3*h))
bb.paste(im[11],(2*w,3*h,2*w+w2,4*h))

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(d)',fill='blue',font=font)
draw.text((0,2*h),u'(g)',fill='blue',font=font)
draw.text((0,3*h),u'(j)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(e)',fill='blue',font=font)
draw.text((w,2*h),u'(h)',fill='blue',font=font)
draw.text((w,3*h),u'(k)',fill='blue',font=font)
draw.text((2*w,0),u'(c)',fill='blue',font=font)
draw.text((2*w,h),u'(f)',fill='blue',font=font)
draw.text((2*w,2*h),u'(i)',fill='blue',font=font)
draw.text((2*w,3*h),u'(l)',fill='blue',font=font)
draw.text((300,150),u'target 360 K',fill='yellow',font=font)
draw.text((300,h+150),u'target 370 K',fill='yellow',font=font)
draw.text((300,2*h+150),u'target 380 K',fill='yellow',font=font)
draw.text((300,3*h+150),u'target 400 K',fill='yellow',font=font)
draw.text((w+300,150),u'source 360 K',fill='yellow',font=font)
draw.text((w+300,h+150),u'source 370 K',fill='yellow',font=font)
draw.text((w+300,2*h+150),u'source 380 K',fill='yellow',font=font)
draw.text((w+300,3*h+150),u'source 400 K',fill='yellow',font=font)

bb.save('figs-Box/tile8.png',dpi=(300,300))
bb.show()


#%% Tile 9 impact of EAD against EAZ
fname0 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-340K-All-h1728-sh.png')
fname2 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-350K-All-h1728-sh.png')
fname4 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-360K-All-h1728-sh.png')
fname6 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-370K-All-h1728-sh.png')
fname8 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-380K-All-h1728-sh.png')
fname10 = join(forw_dir,'target-norm','chart-EAD-Box-theta-FullAMA-target-norm-400K-All-h1728-sh.png')
fname1 = join(forw_dir,'target-norm','chart-EAZ-Box-theta-FullAMA-target-norm-340K-All-h1728-sh.png')
fname3 = join(forw_dir,'target-norm','chart-EAZ-Box-theta-FullAMA-target-norm-350K-All-h1728-sh.png')
fname5 = join(forw_dir,'target-norm','chart-EAZ-Box-theta-FullAMA-target-norm-360K-All-h1728-sh.png')
fname7 = join(forw_dir,'target-norm','chart-EAZ-Box-theta-FullAMA-target-norm-370K-All-h1728-sh.png')
fname9 = join(forw_dir,'target-norm','chart-EAZ-Box-theta-FullAMA-target-norm-380K-All-h1728-sh.png')
fname11 = join(forw_dir,'target-norm','chart-EAZ-Box-theta-FullAMA-target-norm-400K-All-h1728-sh.png')
im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,fname9,fname10,fname11]:
    aa=Image.open(chart)
#    if chart == fname5:
#        im.append(aa.crop([0,0,3610,1132]))
#    else:
    cc = Image.new('RGB',size=(3610,1132),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,6*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.paste(im[6],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.paste(im[8],(0,4*h,w,5*h))
bb.paste(im[9],(w,4*h,2*w,5*h))
bb.paste(im[10],(0,5*h,w,6*h))
bb.paste(im[11],(w,5*h,2*w,6*h))

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arial.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((300,150),u'target 340 K',fill='yellow',font=font)
draw.text((300,h+150),u'target 350 K',fill='yellow',font=font)
draw.text((300,2*h+150),u'target 360 K',fill='yellow',font=font)
draw.text((300,3*h+150),u'target 370 K',fill='yellow',font=font)
draw.text((300,4*h+150),u'target 380 K',fill='yellow',font=font)
draw.text((300,5*h+150),u'target 400 K',fill='yellow',font=font)
draw.text((w+300,150),u'target 340 K',fill='yellow',font=font)
draw.text((w+300,h+150),u'target 350 K',fill='yellow',font=font)
draw.text((w+300,2*h+150),u'target 360 K',fill='yellow',font=font)
draw.text((w+300,3*h+150),u'target 370 K',fill='yellow',font=font)
draw.text((w+300,4*h+150),u'target 380 K',fill='yellow',font=font)
draw.text((w+300,5*h+150),u'target 400 K',fill='yellow',font=font)

bb.save('figs-Box/tile9.png',dpi=(300,300))
bb.show()

#%% Tile 10 sources of EAD against EAZ
fname0 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-340K-All-h1728-sh.png')
fname2 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-350K-All-h1728-sh.png')
fname4 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-360K-All-h1728-sh.png')
fname6 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-370K-All-h1728-sh.png')
fname8 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-380K-All-h1728-sh.png')
fname10 = join(forw_dir,'source-norm','chart-EAD-Box-theta-FullAMA-source-norm-400K-All-h1728-sh.png')
fname1 = join(forw_dir,'source-norm','chart-EAZ-Box-theta-FullAMA-source-norm-340K-All-h1728-sh.png')
fname3 = join(forw_dir,'source-norm','chart-EAZ-Box-theta-FullAMA-source-norm-350K-All-h1728-sh.png')
fname5 = join(forw_dir,'source-norm','chart-EAZ-Box-theta-FullAMA-source-norm-360K-All-h1728-sh.png')
fname7 = join(forw_dir,'source-norm','chart-EAZ-Box-theta-FullAMA-source-norm-370K-All-h1728-sh.png')
fname9 = join(forw_dir,'source-norm','chart-EAZ-Box-theta-FullAMA-source-norm-380K-All-h1728-sh.png')
fname11 = join(forw_dir,'source-norm','chart-EAZ-Box-theta-FullAMA-source-norm-400K-All-h1728-sh.png')
im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,fname9,fname10,fname11]:
    aa=Image.open(chart)
#    if chart == fname5:
#        im.append(aa.crop([0,0,3610,1132]))
#    else:
    cc = Image.new('RGB',size=(3610,1132),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,6*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.paste(im[6],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.paste(im[8],(0,4*h,w,5*h))
bb.paste(im[9],(w,4*h,2*w,5*h))
bb.paste(im[10],(0,5*h,w,6*h))
bb.paste(im[11],(w,5*h,2*w,6*h))

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((300,150),u'source 340 K',fill='yellow',font=font)
draw.text((300,h+150),u'source 350 K',fill='yellow',font=font)
draw.text((300,2*h+150),u'source 360 K',fill='yellow',font=font)
draw.text((300,3*h+150),u'source 370 K',fill='yellow',font=font)
draw.text((300,4*h+150),u'source 380 K',fill='yellow',font=font)
draw.text((300,5*h+150),u'source 400 K',fill='yellow',font=font)
draw.text((w+300,150),u'source 340 K',fill='yellow',font=font)
draw.text((w+300,h+150),u'source 350 K',fill='yellow',font=font)
draw.text((w+300,2*h+150),u'source 360 K',fill='yellow',font=font)
draw.text((w+300,3*h+150),u'source 370 K',fill='yellow',font=font)
draw.text((w+300,4*h+150),u'source 380 K',fill='yellow',font=font)
draw.text((w+300,5*h+150),u'source 400 K',fill='yellow',font=font)

bb.save('figs-Box/tile10.png',dpi=(300,300))
bb.show()

#%% Tile 11 impact of EID-FULL against EIZ-FULL
fname0 = join(forw_dirg,'target-norm','chart-EID-FULL-Box-theta-global-target-norm-340K-All-h1728-sh.png')
fname2 = join(forw_dirg,'target-norm','chart-EID-FULL-Box-theta-global-target-norm-350K-All-h1728-sh.png')
fname4 = join(forw_dirg,'target-norm','chart-EID-FULL-Box-theta-global-target-norm-360K-All-h1728-sh.png')
fname6 = join(forw_dirg,'target-norm','chart-EID-FULL-Box-theta-global-target-norm-370K-All-h1728-sh.png')
fname8 = join(forw_dirg,'target-norm','chart-EID-FULL-Box-theta-global-target-norm-380K-All-h1728-sh.png')
fname10 = join(forw_dirg,'target-norm','chart-EID-FULL-Box-theta-global-target-norm-400K-All-h1728-sh.png')
fname1 = join(forw_dirg,'target-norm','chart-EIZ-FULL-Box-theta-global-target-norm-340K-All-h1728-sh.png')
fname3 = join(forw_dirg,'target-norm','chart-EIZ-FULL-Box-theta-global-target-norm-350K-All-h1728-sh.png')
fname5 = join(forw_dirg,'target-norm','chart-EIZ-FULL-Box-theta-global-target-norm-360K-All-h1728-sh.png')
fname7 = join(forw_dirg,'target-norm','chart-EIZ-FULL-Box-theta-global-target-norm-370K-All-h1728-sh.png')
fname9 = join(forw_dirg,'target-norm','chart-EIZ-FULL-Box-theta-global-target-norm-380K-All-h1728-sh.png')
fname11 = join(forw_dirg,'target-norm','chart-EIZ-FULL-Box-theta-global-target-norm-400K-All-h1728-sh.png')
im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,fname9,fname10,fname11]:
    aa=Image.open(chart)
    h = aa.size[1]
    #if chart == fname5:
    #    im.append(aa.crop([0,0,3446,h]))
    #else:
    cc = Image.new('RGB',size=(3608,h),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,6*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.paste(im[6],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.paste(im[8],(0,4*h,w,5*h))
bb.paste(im[9],(w,4*h,2*w,5*h))
bb.paste(im[10],(0,5*h,w,6*h))
bb.paste(im[11],(w,5*h,2*w,6*h))

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((300,150),u'target 340 K',fill='yellow',font=font)
draw.text((300,h+150),u'target 350 K',fill='yellow',font=font)
draw.text((300,2*h+150),u'target 360 K',fill='yellow',font=font)
draw.text((300,3*h+150),u'target 370 K',fill='yellow',font=font)
draw.text((300,4*h+150),u'target 380 K',fill='yellow',font=font)
draw.text((300,5*h+150),u'target 400 K',fill='yellow',font=font)
draw.text((w+300,150),u'target 340 K',fill='yellow',font=font)
draw.text((w+300,h+150),u'target 350 K',fill='yellow',font=font)
draw.text((w+300,2*h+150),u'target 360 K',fill='yellow',font=font)
draw.text((w+300,3*h+150),u'target 370 K',fill='yellow',font=font)
draw.text((w+300,4*h+150),u'target 380 K',fill='yellow',font=font)
draw.text((w+300,5*h+150),u'target 400 K',fill='yellow',font=font)

bb.save('figs-Box/tile11.png',dpi=(300,300))
bb.show()

#%% Tile 12 sources of EID-FULL against EIZ-FULL
fname0 = join(forw_dirg,'source-norm','chart-EID-FULL-Box-theta-global-source-norm-340K-All-h1728-sh.png')
fname2 = join(forw_dirg,'source-norm','chart-EID-FULL-Box-theta-global-source-norm-350K-All-h1728-sh.png')
fname4 = join(forw_dirg,'source-norm','chart-EID-FULL-Box-theta-global-source-norm-360K-All-h1728-sh.png')
fname6 = join(forw_dirg,'source-norm','chart-EID-FULL-Box-theta-global-source-norm-370K-All-h1728-sh.png')
fname8 = join(forw_dirg,'source-norm','chart-EID-FULL-Box-theta-global-source-norm-380K-All-h1728-sh.png')
fname10 = join(forw_dirg,'source-norm','chart-EID-FULL-Box-theta-global-source-norm-400K-All-h1728-sh.png')
fname1 = join(forw_dirg,'source-norm','chart-EIZ-FULL-Box-theta-global-source-norm-340K-All-h1728-sh.png')
fname3 = join(forw_dirg,'source-norm','chart-EIZ-FULL-Box-theta-global-source-norm-350K-All-h1728-sh.png')
fname5 = join(forw_dirg,'source-norm','chart-EIZ-FULL-Box-theta-global-source-norm-360K-All-h1728-sh.png')
fname7 = join(forw_dirg,'source-norm','chart-EIZ-FULL-Box-theta-global-source-norm-370K-All-h1728-sh.png')
fname9 = join(forw_dirg,'source-norm','chart-EIZ-FULL-Box-theta-global-source-norm-380K-All-h1728-sh.png')
fname11 = join(forw_dirg,'source-norm','chart-EIZ-FULL-Box-theta-global-source-norm-400K-All-h1728-sh.png')
im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,fname9,fname10,fname11]:
    aa=Image.open(chart)
#    if chart == fname5:
#        im.append(aa.crop([0,0,3446,1664]))
#    else:
    cc = Image.new('RGB',size=(3610,1132),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,6*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.paste(im[6],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.paste(im[8],(0,4*h,w,5*h))
bb.paste(im[9],(w,4*h,2*w,5*h))
bb.paste(im[10],(0,5*h,w,6*h))
bb.paste(im[11],(w,5*h,2*w,6*h))

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((300,150),u'source 340 K',fill='yellow',font=font)
draw.text((300,h+150),u'source 350 K',fill='yellow',font=font)
draw.text((300,2*h+150),u'source 360 K',fill='yellow',font=font)
draw.text((300,3*h+150),u'source 370 K',fill='yellow',font=font)
draw.text((300,4*h+150),u'source 380 K',fill='yellow',font=font)
draw.text((300,5*h+150),u'source 400 K',fill='yellow',font=font)
draw.text((w+300,150),u'source 340 K',fill='yellow',font=font)
draw.text((w+300,h+150),u'source 350 K',fill='yellow',font=font)
draw.text((w+300,2*h+150),u'source 360 K',fill='yellow',font=font)
draw.text((w+300,3*h+150),u'source 370 K',fill='yellow',font=font)
draw.text((w+300,4*h+150),u'source 380 K',fill='yellow',font=font)
draw.text((w+300,5*h+150),u'source 400 K',fill='yellow',font=font)

bb.save('figs-Box/tile12.png',dpi=(300,300))
bb.show()

#%% Tile 13 ages of EAD in the target and source domain
fname0 = join(forw_dir,'mage-target','chart-EAD-Box-theta-FullAMA-mage-target-340K-All-h1728-sh.png')
fname2 = join(forw_dir,'mage-target','chart-EAD-Box-theta-FullAMA-mage-target-350K-All-h1728-sh.png')
fname4 = join(forw_dir,'mage-target','chart-EAD-Box-theta-FullAMA-mage-target-360K-All-h1728-sh.png')
fname6 = join(forw_dir,'mage-target','chart-EAD-Box-theta-FullAMA-mage-target-380K-All-h1728-sh.png')
fname8 = join(forw_dir,'mage-target','chart-EAD-Box-theta-FullAMA-mage-target-400K-All-h1728-sh.png')
fname1 = join(forw_dir,'mage-source','chart-EAD-Box-theta-FullAMA-mage-source-340K-All-h1728-sh.png')
fname3 = join(forw_dir,'mage-source','chart-EAD-Box-theta-FullAMA-mage-source-350K-All-h1728-sh.png')
fname5 = join(forw_dir,'mage-source','chart-EAD-Box-theta-FullAMA-mage-source-360K-All-h1728-sh.png')
fname7 = join(forw_dir,'mage-source','chart-EAD-Box-theta-FullAMA-mage-source-380K-All-h1728-sh.png')
fname9 = join(forw_dir,'mage-source','chart-EAD-Box-theta-FullAMA-mage-source-400K-All-h1728-sh.png')
im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,fname9]:
    aa=Image.open(chart)
    #im.append(aa)
    cc = Image.new('RGB',size=(3610,1136),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
[w,h]=im[0].size
bb=Image.new('RGB',(2*w,5*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))
bb.paste(im[6],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))
bb.paste(im[8],(0,4*h,w,5*h))
bb.paste(im[9],(w,4*h,2*w,5*h))

# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((0,4*h),u'(i)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((w,5*h),u'(j)',fill='blue',font=font)
draw.text((1500,150),u'target 340 K',fill='black',font=font)
draw.text((1500,h+150),u'target 350 K',fill='black',font=font)
draw.text((2650,2*h+700),u'target 360 K',fill='black',font=font)
draw.text((2650,3*h+700),u'target 380 K',fill='black',font=font)
draw.text((2650,4*h+700),u'target 400 K',fill='black',font=font)
draw.text((w+300,500),u'source 340 K',fill='black',font=font)
draw.text((w+300,h+500),u'source 350 K',fill='black',font=font)
draw.text((w+300,2*h+500),u'source 360 K',fill='black',font=font)
draw.text((w+300,3*h+500),u'source 380 K',fill='black',font=font)
draw.text((w+300,4*h+500),u'source 400 K',fill='black',font=font)

bb.save('figs-Box/tile13.png',dpi=(300,300))
bb.show()

#%% Tile 14 tiling of images for the backward calculations in the same way as 
# tile 13 (minus 340K) for the forward calculations 

fname0 = join(back_dir,'EAD-JAS','chart-EAD-age-target-Jul-Aug-Sep-350K.png')
fname1 = join(back_dir,'EAD-JAS','chart-EAD-age-target-Jul-Aug-Sep-360K.png')
fname2 = join(back_dir,'EAD-JAS','chart-EAD-age-target-Jul-Aug-Sep-380K.png')
fname3 = join(back_dir,'EAD-JAS','chart-EAD-age-target-Jul-Aug-Sep-400K.png')
fname4 = join(back_dir,'EAD-JAS','chart-EAD-age-source-Jul-Aug-Sep-350K.png')
fname5 = join(back_dir,'EAD-JAS','chart-EAD-age-source-Jul-Aug-Sep-360K.png')

im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5,fname6,fname7]:
    aa=Image.open(chart)
    cc = Image.new('RGB',size=(3644,1136),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(cc)
       
[w,h] = im[0].size
bb=Image.new('RGB',(2*w,4*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[4],(w,0,2*w,h))
bb.paste(im[1],(0,h,w,2*h))
bb.paste(im[5],(w,h,2*w,2*h))
bb.paste(im[2],(0,2*h,w,3*h))
bb.paste(im[6],(w,2*h,2*w,3*h))
bb.paste(im[3],(0,3*h,w,4*h))
bb.paste(im[7],(w,3*h,2*w,4*h))

# add labels
# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((0,3*h),u'(g)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((w,3*h),u'(h)',fill='blue',font=font)
draw.text((1500,150),u'target 350 K',fill='black',font=font)
draw.text((1500,h+150),u'target 360 K',fill='black',font=font)
draw.text((2650,2*h+700),u'target 380 K',fill='black',font=font)
draw.text((2650,3*h+700),u'target 400 K',fill='black',font=font)
draw.text((w+300,500),u'source 350 K',fill='black',font=font)
draw.text((w+300,h+500),u'source 360 K',fill='black',font=font)
draw.text((w+300,2*h+500),u'source 380 K',fill='black',font=font)
draw.text((w+300,3*h+500),u'source 400 K',fill='black',font=font)

bb.save('figs-Box/tile14.png',dpi=(300,300))
bb.show()

#%% Tile 15 tiling of images for the variability study

fname0 = join(forw_target,'chart-EAD-Box-theta-FullAMA-Jul-11-target-380K-sh.png')
fname1 = join(forw_target,'chart-EAD-Box-theta-FullAMA-Jul-21-target-380K-sh.png')
fname2 = join(forw_target,'chart-EAD-Box-theta-FullAMA-Aug-01-target-380K-sh.png')
fname3 = join(forw_target,'chart-EAD-Box-theta-FullAMA-Aug-11-target-380K-sh.png')
fname4 = join(forw_target,'chart-EAD-Box-theta-FullAMA-Aug-21-target-380K-sh.png')
fname5 = join(forw_dir,'target','chart-EAD-Box-theta-FullAMA-target-370K-All-h1728-sh.png')


im = []
for chart in [fname0,fname1,fname2,fname3,fname4,fname5]:
    aa=Image.open(chart)
    cc = Image.new('RGB',size=(3644,1136),color=(255,255,255))
    cc.paste(aa,(0,0,aa.size[0],aa.size[1]))
    im.append(aa)
       
[w,h] = im[0].size
bb=Image.new('RGB',(2*w,3*h))
bb.paste(im[0],(0,0,w,h))
bb.paste(im[1],(w,0,2*w,h))
bb.paste(im[2],(0,h,w,2*h))
bb.paste(im[3],(w,h,2*w,2*h))
bb.paste(im[4],(0,2*h,w,3*h))
bb.paste(im[5],(w,2*h,2*w,3*h))

# add labels
# add labels
draw = ImageDraw.Draw(bb)
font  = ImageFont.truetype("arialbd.ttf",96)
draw.text((0,0),u'(a)',fill='blue',font=font)
draw.text((0,h),u'(c)',fill='blue',font=font)
draw.text((0,2*h),u'(e)',fill='blue',font=font)
draw.text((w,0),u'(b)',fill='blue',font=font)
draw.text((w,h),u'(d)',fill='blue',font=font)
draw.text((w,2*h),u'(f)',fill='blue',font=font)
draw.text((2650,800),u'Jul 11-20',fill='yellow',font=font)
draw.text((2650,h+800),u'Aug 01-10',fill='yellow',font=font)
draw.text((2650,2*h+800),u'Aug 21-31',fill='yellow',font=font)
draw.text((w+2650,800),u'Jul 21-31',fill='yellow',font=font)
draw.text((w+2650,h+800),u'Aug 11-20',fill='yellow',font=font)
draw.text((w+2650,2*h+800),u'All summer',fill='yellow',font=font)

bb.save('figs-Box/tile15.png',dpi=(300,300))
bb.show()