# -*- coding: utf-8 -*-
"""

Script to generate the contour of the Tibetan Plateau

Created on 12 May 2020

@author: Bernard Legras
"""
import pickle,gzip
import matplotlib.pyplot as plt

#We open the existing mask
mask0=pickle.load(gzip.open('MaskCartopy2-STCfine.pkl','rb'))
mask=mask0['mask']
lonc=mask0['lons']
latc=mask0['lats']
# make an intermediate mask that contains only the Tibetan Plateau
tibet = mask == mask0['regcode']['TibetanPlateau']
# make a contour plot from this new mask
cs=plt.contour(lonc,latc,tibet,[1,])
plt.show()
# extract the vertices of the main contour
tibetc = cs.collections[0].get_paths()[0].vertices
plt.plot(tibetc[:,0],tibetc[:,1])
plt.show()
# save the vertices
pickle.dump(tibetc,open('TibetContour.pkl','wb'))
