#! /usr/bin/env python
from astropy.io import fits
import matplotlib.pyplot as plt
from aplpy import FITSFigure

filename = '../sofia_data/F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits'
hawc = fits.open(filename)

stokes_i = hawc['STOKES I']
fig = FITSFigure(stokes_i)
fig.show_colorscale()
plt.show()
