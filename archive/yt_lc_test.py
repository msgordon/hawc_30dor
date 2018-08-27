#! /usr/bin/env python
import numpy as np
from astropy.io import fits
from astropy.table import Table
#from yt.utilities.lib.line_integral_convolution import line_integral_convolution_2d
from liccy import line_integral_convolution_2d
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from astropy.nddata import block_replicate
from scipy.interpolate import griddata


def upsample(data,block, method='cubic'):
    points,values = zip(*np.ndenumerate(data))
    points = np.array(points)
    values = np.array(values)

    grid_x,grid_y = np.mgrid[0:data.shape[0]:data.shape[0]*block*1j,
                             0:data.shape[1]:data.shape[1]*block*1j]

    real_idx = np.array([~np.isnan(v) for v in values])

    points = points[real_idx]
    values = values[real_idx]

    pf = griddata(points,values,(grid_x,grid_y),method=method)
    return pf
    

hdu = fits.open('sofia_data/F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits')

block = 5

stokes_i = hdu['STOKES I'].data

#stokes_i = block_replicate(stokes_i,block)

p = hdu['DEBIASED PERCENT POL'].data
theta = hdu['ROTATED POL ANGLE'].data

#plt.imshow(theta,origin='lower')

stokes_if = upsample(stokes_i,block)
pf = upsample(p, block)
tf = upsample(theta,block)
#plt.figure()
#plt.imshow(tf,origin='lower')
#plt.show()

u = pf*np.cos(tf*np.pi/180)
v = pf*np.sin(tf*np.pi/180)

vectors = np.concatenate((u[...,np.newaxis],
                          v[...,np.newaxis]),axis=2)

kernellen = 50
texture = np.random.rand(pf.shape[0],pf.shape[1]).astype(np.double)

kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)

lic_data = line_integral_convolution_2d(vectors, texture, kernel)
lic_data = lic_data / np.nanmax(lic_data)

lim = (0,1)
alpha = 0.75
cmap = 'binary'

cmapf = plt.get_cmap(cmap)
new_cmap = colors.LinearSegmentedColormap.from_list('trunc(%s,0.3,0.8)'%
                                                    cmap,cmapf(np.linspace(0.3,1.0,100)))

lic_data_clip = np.clip(lic_data,lim[0],lim[1])

lic_data_rgba = cm.ScalarMappable(norm=None,cmap=cmap).to_rgba(lic_data_clip)
lic_data_clip_rescale = (lic_data_clip - lim[0]) / (lim[1] - lim[0])
lic_data_rgba[...,3] = lic_data_clip_rescale * alpha

plt.imshow(stokes_if,origin='lower')
plt.imshow(lic_data_rgba,origin='lower',cmap=cmap,alpha=alpha)

plt.show()
exit()

plt.figure()
#plt.imshow(lic_data_clip,origin='lower')
plt.imshow(stokes_if,origin='lower')

#plt.figure()
from skimage.feature import canny
from skimage.measure import find_contours

#lines = canny(lic_data_rgba)
lines = canny(lic_data_clip)
contours = find_contours(lines,0.8, fully_connected='high')
for contour in contours:
    plt.plot(contour[:,1], contour[:,0],alpha=0.9,color='w',lw=0.5)
#plt.imshow(lines,origin='lower')


plt.show()

exit()

#print(np.indices(p.shape))


exit()
#xx = np.arange(0,stokes_i.shape[0])
#yy = np.arange(0,stokes_i.shape[1])

#xx,yy = np.meshgrid(xx,yy)
yy,xx = np.mgrid[:p.shape[1],:p.shape[0]]

pf = interp2d(xx,yy,p,kind='linear')

XX = np.arange(0,stokes_i.shape[0],0.5)
YY = np.arange(0,stokes_i.shape[1],0.5)

pnew = pf(XX,YY)
plt.figure()
plt.imshow(pnew,origin='lower')
plt.show()
exit()


#p = block_replicate(p,block)
#theta = block_replicate(theta,block)

u = p*np.cos(theta*np.pi/180)
v = p*np.sin(theta*np.pi/180)

vectors = np.concatenate((u[...,np.newaxis],
                          v[...,np.newaxis]),axis=2)

kernellen = 50
texture = np.random.rand(stokes_i.shape[0],stokes_i.shape[1]).astype(np.double)

kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)

lic_data = line_integral_convolution_2d(vectors, texture, kernel)
lic_data = lic_data / np.nanmax(lic_data)

lim = (0.2,0.6)
alpha = 0.8
cmap = 'binary_r'

lic_data_clip = np.clip(lic_data,lim[0],lim[1])

lic_data_rgba = cm.ScalarMappable(norm=None,cmap=cmap).to_rgba(lic_data_clip)
lic_data_clip_rescale = (lic_data_clip - lim[0]) / (lim[1] - lim[0])
lic_data_rgba[...,3] = lic_data_clip_rescale * alpha

plt.imshow(stokes_i,origin='lower')
plt.imshow(lic_data_rgba,origin='lower',cmap=cmap,alpha=alpha)
plt.show()
