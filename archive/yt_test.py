#! /usr/bin/env python
import numpy as np
from astropy.io import fits
from astropy.table import Table
import yt
from astropy.wcs import WCS

hdu = fits.open('sofia_data/F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits')

stokes_i = hdu['STOKES I']
density = np.expand_dims(stokes_i.data,axis=-1)

p = hdu['DEBIASED PERCENT POL'].data
theta = hdu['ROTATED POL ANGLE'].data

u = p*np.cos(theta*np.pi/180)
v = p*np.sin(theta*np.pi/180)
w = np.zeros(p.shape)

'''
data = {'particle_velocity_x':u,
        'particle_velocity_y':v,
        'particle_velocity_z':w,
        'density':density}
'''

#density = np.expand_dims(p,axis=-1)

data = {'x-velocity':u,
        'y-velocity':v,
        'z-velocity':w,
        'density':density}


#ds = yt.visualization.fits_image.FITSImageData(data,wcs=WCS(stokes_i.header))

#s = yt.FITSSlice(ds,'z','density')
#s.show()


ds = yt.load_uniform_grid(data,(p.shape[0],p.shape[1],1))

print(dir(ds.fields.stream));exit()
print(ds.field_list)
print(ds.derived_field_list)

s = yt.SlicePlot(ds,'z','density')
#s.annotate_streamlines('density_gradient_x','density_gradient_y',factor=16)
s.annotate_streamlines('x-velocity','y-velocity',factor=16)

#s.set_origin('lower-left-domain')
#print(dir(s))
#print(s)
#print(ds.field_list)
#print(ds.derived_field_list)
#exit()
#s.annotate_line_integral_convolution('x', 'y', alpha=0.4)#,  lim=(0.2,0.8),alpha=0.4)
s.save()
