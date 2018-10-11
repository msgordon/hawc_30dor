
# coding: utf-8

# # HAWC+ Data Analysis Techniques
# 
# In this jupyter cookbook, we will explore the [HAWC+](https://www.sofia.usra.edu/science/instruments/hawc) data cube and describe some of the basic analysis techniques involving imaging polarimetry data.
# 
# This cookbook follows the SOFIA press release of 30 Doradus observations: [SOFIA Reveals Never-Before-Seen Magnetic Field Details](https://www.sofia.usra.edu/multimedia/science-results-archive/sofia-reveals-never-seen-magnetic-field-details).
# 
# The Level 4 reduced data from this program has been released immediately to the public and is available on the [SOFIA Data Cycle System (DCS)](https://dcs.sofia.usra.edu/).  This notebook will guide the reader through downloading the 30 Doradus data with a walkthrough of basic analysis techniques with `python`.

# # Downloading HAWC+ Data
# 
# - If you do not yet have a DCS account, register for one at [https://dcs.sofia.usra.edu/userSupport/registration.jsp](https://dcs.sofia.usra.edu/userSupport/registration.jsp)
# - Log into DCS: [https://dcs.sofia.usra.edu](https://dcs.sofia.usra.edu)
# - Go to [Search Science Archive](https://dcs.sofia.usra.edu/dataRetrieval/SearchScienceArchiveInfoBasic.jsp)
# - Fill in:
#   - Instrument: `HAWC_PLUS` from drop-down menu
#   - Processing State: `LEVEL_4` from drop-down menu
#   - Target: 30Dor
#   - Click the `SIMBAD Position` button
#   - Change `Spatial Search Radius` to `600` arcsec
#   - Click the `Submit` button
# - After the results load, select the checkboxes next to each of the six rows of the table.
# - Click `Get Selected Data in Current Page`
# - Click `Request Data Bundle`
# - After a few minutes, an email with a download link will be sent to the email address associated with your DCS account.
# - For more information, consult the HAWC+ Data Handbook accessible at [https://www.sofia.usra.edu/science/proposing-and-observing/data-products/data-resources](https://www.sofia.usra.edu/science/proposing-and-observing/data-products/data-resources)

# # SOFIA Data Organization
# After downloading the SOFIA DCS bundle to your working directory you will want to unzip it, which will produce a directory structure like this:

# ```console
# .
# └── sofia_data
#     ├── level4
#     │   └── p5813
#     │       └── F0484_HA_POL_7600018_HAWCHWPC_PMP_022-114.fits
#     └── missions
#         ├── 2018-07-05_HA_F481
#         │   └── p5827
#         │       └── F0481_HA_POL_7600012_HAWDHWPD_PMP_050-083.fits
#         ├── 2018-07-07_HA_F483
#         │   └── p5646
#         │       └── F0483_HA_POL_7600014_HAWCHWPC_PMP_022-065.fits
#         ├── 2018-07-11_HA_F484
#         │   └── p5648
#         │       └── F0484_HA_POL_7600017_HAWCHWPC_PMP_065-114.fits
#         └── 2018-07-12_HA_F485
#             └── p5658
#                 ├── g1
#                 │   └── F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits
#                 └── g2
#                     └── F0485_HA_POL_7600019_HAWEHWPE_PMP_055-075.fits
# ```

# Note that each file represents observations with a different filter.  However, two observations were made with the same filter (HAWC C, $89\,\mathrm{\mu m}$).  These files, `F0483_HA_POL_7600014_HAWCHWPC_PMP_022-065.fits` and `F0484_HA_POL_7600017_HAWCHWPC_PMP_065-114.fits`, were combined into one: `level4->p5813->F0484_HA_POL_7600018_HAWCHWPC_PMP_022-114.fits`.

# You can choose to keep the `fits` files nested, or copy them into one directory.

# For the purpose of this basic analysis, though, let us dump all the files into one `sofia_data` directory:

# ```console
# .
# └── sofia_data
#     ├── F0481_HA_POL_7600012_HAWDHWPD_PMP_050-083.fits
#     ├── F0483_HA_POL_7600014_HAWCHWPC_PMP_022-065.fits
#     ├── F0484_HA_POL_7600017_HAWCHWPC_PMP_065-114.fits
#     ├── F0484_HA_POL_7600018_HAWCHWPC_PMP_022-114.fits
#     ├── F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits
#     └── F0485_HA_POL_7600019_HAWEHWPE_PMP_055-075.fits
# ```

# # Data Structure
# 
# For this analysis, we require the standard [numpy/scipy/matplotlib stack](https://scipy.org/install.html) as well the [astropy](http://docs.astropy.org/en/stable/) and [aplpy](https://aplpy.readthedocs.io/en/stable/index.html) modules.
# 
# With just a few lines of code, we can explore the HAWC+ `fits` data cubes and plot the images.

# In[1]:


from astropy.io import fits

efile = 'sofia_data/F0485_HA_POL_7600019_HAWEHWPE_PMP_055-075.fits'
dfile = 'sofia_data/F0481_HA_POL_7600012_HAWDHWPD_PMP_050-083.fits'
cfile = 'sofia_data/F0484_HA_POL_7600018_HAWCHWPC_PMP_022-114.fits'


afile = 'sofia_data/F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits'
hawc = fits.open(afile)
hawc.info()


# We can see above the data structure of the multi-extension `fits` files.  Each file contains 19 extensions which encapsulates all of the Stokes parameters in a single package.

# ## Stokes I
# Stokes $I$---the zeroth extension in the `fits` file---represents the total intensity of the image, where $I^2 = Q^2 + U^2$.
# 
# 
# Let us go ahead and plot this extension:

# In[2]:


import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'notebook')
# ^jupyter magic for inline plots
from aplpy import FITSFigure

# set colormap for all plots
cmap = 'rainbow'

stokes_i = hawc['STOKES I']               # or hawc[0]. Note the extension is from the hawc.info() table above

fig = plt.figure(figsize=(7,7))

axs = FITSFigure(stokes_i, figure=fig)    # load HDU into aplpy figure
axs.show_colorscale(cmap=cmap)            # display the data with WCS projection and chosen colormap

# FORMATTING
axs.set_tick_labels_font(size='small')
axs.set_axis_labels_font(size='small')

# Add colorbar
axs.add_colorbar()
axs.colorbar.set_axis_label_text('Flux (Jy/pix)')


# ## Stokes Q and U
# Similarly, we can plot the Stokes Q and Stokes U images:

# In[3]:


stokes_q = hawc['STOKES Q']
stokes_u = hawc['STOKES U']

axq = FITSFigure(stokes_q, subplot=(1,2,1))  # generate FITSFigure as subplot to have two axes together
axq.show_colorscale(cmap=cmap)               # show Q


axu = FITSFigure(stokes_u, subplot=(1,2,2),
                 figure=plt.gcf())
axu.show_colorscale(cmap=cmap)               # show U

# FORMATTING
axq.set_title('Stokes Q')
axu.set_title('Stokes U')
#axu.axis_labels.hide_y()                     # hide axis ticklabels for U figure
#axu.tick_labels.hide_y()
axu.axis_labels.set_yposition('right')
axu.tick_labels.set_yposition('right')
axq.set_tick_labels_font(size='small')
axq.set_axis_labels_font(size='small')
axu.set_tick_labels_font(size='small')
axu.set_axis_labels_font(size='small')


# We can additionally plot the associated error maps for each extension.

# In[4]:


stokes_q = hawc['STOKES Q']
error_q = hawc['ERROR Q']

axq = FITSFigure(stokes_q, subplot=(1,2,1))  # generate FITSFigure as subplot to have two axes together
axq.show_colorscale(cmap=cmap)               # show Q


axe = FITSFigure(error_q, subplot=(1,2,2), figure=plt.gcf())
axe.show_colorscale(cmap=cmap)               # show error

# FORMATTING
axq.set_title('Stokes Q')
axe.set_title('Error Q')
axq.axis_labels.hide()                       # hide axis/tick labels
axe.axis_labels.hide()
axq.tick_labels.hide()
axe.tick_labels.hide()


# ## Polarized Intensity $I_p$
# Level 4 HAWC+ additionally provides extensions with the polarization percentage ($p$), angle ($\theta$), and their associated errors ($\sigma$).
# 
# Percent polarization ($p$) and error ($\sigma_p$) are calculated as:
# 
# $p = 100\sqrt{\left(\frac{Q}{I}\right)^2+\left(\frac{U}{I}\right)^2}$
# 
# $\sigma_p = \frac{100}{I}\sqrt{\frac{1}{(Q^2+U^2)}\left[(Q\,\sigma_Q)^2+(U\,\sigma_U)^2+2QU\,\sigma_{QU}\right]+\left[\left(\frac{Q}{I}\right)^2+\left(\frac{U}{I}\right)^2\right]\sigma_I^2-2\frac{Q}{I}\sigma_{QI}-2\frac{U}{I}\sigma_{UI}}$ .
# 
# Note that $p$ here represents the **percent** polarization as opposed to the more typical convention for $p$ as the **fractional** polarization.
# 
# Maps of these data are found in extensions 7 (PERCENT POL) and 9 (ERROR PERCENT POL).
# 
# Polarized intensity, $I_p$, can then be calculated as $I_p = \frac{I\times p}{100}$, which is included in extension 13 (POL FLUX).
# 
# Also included is the debiased polarization percentage ($p^\prime$) calculated as:
# 
# $p^\prime=\sqrt{p^2-\sigma_p^2}$, found in extension 8 (DEBIASED PERCENT POL).
# 
# We similarly define the debiased polarized intensity as $I_{p^\prime} = \frac{I\times p^\prime}{100}$, which is included in extension 15 (DEBIASED POL FLUX).

# In[5]:


stokes_ip = hawc['DEBIASED POL FLUX']

axi = FITSFigure(stokes_i, subplot=(1,2,1))
axi.show_colorscale(cmap=cmap)               # show I


axp = FITSFigure(stokes_ip, subplot=(1,2,2), figure=plt.gcf())
axp.show_colorscale(cmap=cmap)               # show Ip

# FORMATTING
axi.set_title(r'$I$')
axp.set_title(r'$I_{p^\prime}$')
#axi.axis_labels.hide()                       # hide axis/tick labels
#axp.axis_labels.hide()
#axi.tick_labels.hide()
#axp.tick_labels.hide()
axp.axis_labels.set_yposition('right')
axp.tick_labels.set_yposition('right')
axi.set_tick_labels_font(size='small')
axi.set_axis_labels_font(size='small')
axp.set_tick_labels_font(size='small')
axp.set_axis_labels_font(size='small')


# # Plotting Polarization Vectors
# 
# From the $Q$ and $U$ maps, the polarization angle $\theta$ is calculated in the standard way:
# 
# $\theta = \frac{90}{\pi}\,\mathrm{tan^{-1}}\left(\frac{U}{Q}\right)$
# 
# with associated error:
# 
# $\sigma_\theta = \frac{90}{\pi\left(Q^2+U^2\right)}\sqrt{\left(Q\sigma_Q\right)^2+\left(U\sigma_U\right)^2-2QU\sigma_{QU}}$
# 
# The angle map is stored in extension 10 (POL ANGLE), with its error in extension 12 (ERROR POL ANGLE).  
# 
# However, these angles are relative to detector coordinates.  The angle we are more interested in is the angle on the sky.  As part of the HAWC+ reduction pipeline, $\theta$ is corrected for the vertical position angle of the instrument on the sky, the angle of the HWP plate, as well as an offset angle that is calibrated to each filter configuration.  This correction angle is applied to $\theta\rightarrow\theta^\prime$ and is saved to extension 11 (ROTATED POL ANGLE).  This also affects the Stokes $Q$ and $U$ parameters, and so this factor has already been rolled into the STOKES Q and STOKES U extensions (and their corresponding error maps) in the HAWC+ data cube.
# 
# We can now use the $p^\prime$ and $\theta^\prime$ maps to plot the polarization vectors.  First, however, let us make a quality cut.  Rather than defining a $\sigma$ cut on the polarization vectors themselves, it is more useful to define a signal-to-noise cut on total intensity, $I$, the measured the quantity.
# 
# Returning to the definition of $p$ for the moment:
# 
# $p = \frac{100\sqrt{Q^2+U^2}}{I}$
# 
# Let's assume the errors in $Q$ and $U$ are comparable such that there are no covariant (cross) terms in the error expansion.  Essentially, we define a quantity $x\equiv\sqrt{Q^2+U^2}$ so that:
# 
# \begin{align}
#     p & = \frac{100\sqrt{Q^2+U^2}}{I} = \frac{x}{I} \\
#     \left(\frac{\sigma_p}{p}\right)^2 & = \left(\frac{\sigma_x}{x}\right)^2 + \left(\frac{\sigma_I}{I}\right)^2
# \end{align}
# 
# If $p\sim 1\%\times I$, then
# 
# \begin{align}
#     p & = 0.01 = \frac{x}{I} \\
#     x & = 0.01\,I \\
#     \sigma_x & = 0.01\,\sigma_I
# \end{align}
#     
# \begin{equation*}
#     \Rightarrow\frac{\sigma_x}{x} = \frac{\sigma_I}{I}
# \end{equation*}
# 
# Therefore,
# \begin{equation*}
#     \left(\frac{\sigma_p}{p}\right)^2 \sim 2\,\left(\frac{\sigma_I}{I}\right)^2
# \end{equation*}
# Inverting this, where $\frac{\sigma_x}{x}$ is the S/N of that quantity,
# \begin{align*}
#     \left(\mathrm{S/N}\right)_p & \sim \frac{1}{\sqrt{2}}\,\left(\mathrm{S/N}\right)_I \\
#     \left(\mathrm{S/N}\right)_I & \sim \sqrt{2}\left(\mathrm{S/N}\right)_p \\
#     & \sim \sqrt{2}\left(\frac{p}{\sigma_p}\right)
# \end{align*}
# 
# One way to think about how to proceed from here is to ask, if a fractional polarization is measured $\sim1\%$, what is the maximum error we would tolerate in that quantity? For an error of $0.5\%$ we have:
# \begin{align}
#     \left(\mathrm{S/N}\right)_I & \sim \sqrt{2}\left(\frac{p}{\sigma_p}\right) \sim \sqrt{2}\frac{1}{0.5\%} \\
#     & \sim \frac{\sqrt{2}}{0.005} \sim 283
# \end{align}
# 
# So, therefore if we desire an accuracy of $\sigma_p\sim0.5\%$, we require a S/N in total intensity $I$ of $\sim283$. 
# 
# We perform the following steps:
# 1.  use the Stokes $I$ image as the background for the vector plot
# 2.  perform a quality cut on Stokes $I/\sigma_I > 100$ to make a mask
# 3.  mask out low S/N vectors
# 4.  plot remaining polarization vectors
# 5.  add contours to better visualize changes in flux across the map

# In[6]:


from astropy.io import fits
import numpy as np
from aplpy import FITSFigure

def make_polmap(filename, title=None, figure=None, subplot=(1,1,1)):
    hawc = fits.open(filename)
    p = hawc['DEBIASED PERCENT POL']    # %
    theta = hawc['ROTATED POL ANGLE']   # deg
    stokes_i = hawc['STOKES I']         # I
    error_i = hawc['ERROR I']           # error I

    # 1. plot Stokes I
    #  convert from Jy/pix to Jy/sq. arcsec
    pxscale = stokes_i.header['CDELT2']*3600  # map scale in arcsec/pix
    stokes_i.data /= pxscale**2
    error_i.data /= pxscale**2

    fig = FITSFigure(stokes_i, figure=figure, subplot=subplot)

    # 2. perform S/N cut on I/\sigma_I
    err_lim = 100
    mask = np.where(stokes_i.data/error_i.data < err_lim)

    # 3. mask out low S/N vectors by setting masked indices to NaN
    p.data[mask] = np.nan

    # 4. plot vectors
    scalevec = 0.4  # 1pix = scalevec * 1% pol          # scale vectors to make it easier to see 
    fig.show_vectors(p, theta, scale=scalevec, step=2)  # step size = display every 'step' vectors
                                                        #   step size of 2 is effectively Nyquist sampling
                                                        #   --close to the beam size

    # 5. plot contours
    ncontours = 30
    fig.show_contour(stokes_i, cmap=cmap, levels=ncontours,
                     filled=True, smooth=1, kernel='box')
    fig.show_contour(stokes_i, colors='gray', levels=ncontours,
                     smooth=1, kernel='box', linewidths=0.3)

    # Show image
    fig.show_colorscale(cmap=cmap)
    
    # If title, set it
    if title:
        fig.set_title(title)

    # Add colorbar
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')

    # Add beam indicator
    fig.add_beam(facecolor='red', edgecolor='black',
                 linewidth=2, pad=1, corner='bottom left')
    fig.add_label(0.02, 0.02, 'Beam FWHM',
                  horizontalalignment='left', weight='bold',
                  relative=True, size='small')

    # Add vector scale
    #   polarization vectors are displayed such that 'scalevec' * 1% pol is 1 pix long
    #   must translate pixel size to angular degrees since the 'add_scalebar' function assumes a physical scale
    vectscale = scalevec * pxscale/3600
    fig.add_scalebar(5 * vectscale, "p = 5%",corner='top right',frame=True)
    
    # FORMATTING
    fig.set_tick_labels_font(size='small')
    fig.set_axis_labels_font(size='small')
    
    return stokes_i, p, mask, fig


# In[7]:


stokes_i, p, mask, fig = make_polmap(afile, title='A')


# In[8]:


fig = FITSFigure(p)

# Show image
fig.show_colorscale(cmap=cmap)

# Plot contours
ncontours = 30
fig.show_contour(stokes_i, colors='gray', levels=ncontours,
                 smooth=1, kernel='box', linewidths=0.3)

# Add colorbar
fig.add_colorbar()
fig.colorbar.set_axis_label_text('$p$ (%)')


# In[9]:


files = [afile,cfile,dfile,efile]
titles = ['A','C','D','E']

for file, title in zip(files,titles):
    make_polmap(file,title)


# In[10]:


# reproject
from reproject import reproject_exact

a_orig = fits.open(afile)['STOKES I']
hawc_a_header = a_orig.header

new_c, footprint = reproject_exact(cfile,output_projection=hawc_a_header,hdu_in='STOKES I')

c_repr = fits.PrimaryHDU(new_c,header=hawc_a_header)

afig = FITSFigure(a_orig,subplot=(1,2,1))
afig.show_colorscale(cmap=cmap)

cfig = FITSFigure(c_repr, subplot=(1,2,2), figure=plt.gcf())
cfig.show_colorscale(cmap=cmap)

# FORMATTING
afig.set_title('A')
cfig.set_title('C')
cfig.axis_labels.set_yposition('right')
cfig.tick_labels.set_yposition('right')
afig.set_tick_labels_font(size='small')
afig.set_axis_labels_font(size='small')
cfig.set_tick_labels_font(size='small')
cfig.set_axis_labels_font(size='small')


# In[11]:


new_d, footprint = reproject_exact(dfile,output_projection=hawc_a_header,hdu_in='STOKES I')
new_e, footprint = reproject_exact(efile,output_projection=hawc_a_header,hdu_in='STOKES I')

d_repr = fits.PrimaryHDU(new_d,header=hawc_a_header)
e_repr = fits.PrimaryHDU(new_e,header=hawc_a_header)

dfig = FITSFigure(d_repr, subplot=(1,2,1))
dfig.show_colorscale(cmap=cmap)

efig = FITSFigure(e_repr, subplot=(1,2,2),figure=plt.gcf())
efig.show_colorscale(cmap=cmap)

# FORMATTING
dfig.set_title('D')
efig.set_title('E')
efig.axis_labels.set_yposition('right')
efig.tick_labels.set_yposition('right')
dfig.set_tick_labels_font(size='small')
dfig.set_axis_labels_font(size='small')
efig.set_tick_labels_font(size='small')
efig.set_axis_labels_font(size='small')


# # Temperature map
# 
# In this section, we demonstrate the use of `astropy` modeling to fit a graybody function through all four HAWC+ images. We define an emissivity function $\varepsilon_\lambda$ which we represent as a powerlaw with spectral index $\beta \leq 0$. This function then adds a wavelength-dependent scaling factor to the Planck function:
# \begin{align}
#     \varepsilon_\lambda & = \varepsilon_0\left(\frac{\lambda}{\lambda_0}\right)^{\beta} \\
#     I_{\lambda,T} & = \varepsilon_\lambda\times B_{\lambda}(T)
# \end{align}
# 
# Ideally, $\varepsilon$ would also have a temperature dependence, which we neglect for this model.  Note that for $\beta=0$ the wavelength-dependence disappears, and $I_{\lambda,T}$ is simply the Planck function (modulo some scaling factor $\varepsilon_0$ which we set to 1 in the following code).
# 
# To construct our model, we use the `astropy` `PowerLaw1D` and `BlackBody1D` model objects.  Note that the $\alpha$ quantity in `PowerLaw1D` is $-\beta$ in our formulation above.  We set the amplitude, $\varepsilon_0$, to 1 and fix all values upon instantiation.  This means that these quantities are not left as free parameters for fitting later.

# In[22]:


#from astropy.modeling.models import PowerLaw1D, BlackBody1D
from astropy.modeling.blackbody import blackbody_lambda
import astropy.units as u
from astropy.visualization import quantity_support
import astropy.constants as const

def GrayBody(waves, temperature, bolometric_flux=1*u.erg/(u.cm**2*u.s), wave_0=53*u.micron, beta=-1.5,
             return_units=u.erg/(u.s*u.cm**2)):

    eps = (waves/wave_0)**beta

    # We normalize the returned blackbody so that the integral would be
    # unity, and we then multiply by the bolometric flux. A normalized
    # blackbody has f_lam = pi * B_nu / (sigma * T^4), which is what we
    # calculate here.
    bb = (np.pi * u.sr * blackbody_lambda(waves, temperature) /
          const.sigma_sb / temperature ** 4)

    flambda = bolometric_flux * eps * bb

    return flambda.to(return_units, u.spectral_density(waves))

'''
exit()
# Generate a number of epsilon models with different alpha (-beta) values to visualize
model0 = PowerLaw1D(amplitude=1, x_0=53*u.micron, alpha=0,
                    fixed={'ampitude':True,'x_0':True, 'alpha':True})
model1 = PowerLaw1D(amplitude=1, x_0=53*u.micron, alpha=1,
                    fixed={'ampitude':True,'x_0':True, 'alpha':True})
model2 = PowerLaw1D(amplitude=1, x_0=53*u.micron, alpha=2,
                    fixed={'ampitude':True,'x_0':True, 'alpha':True})


# Generate blackbody model and multiply to form compound model
BB = BlackBody1D(temperature=500*u.K)

model0 *= BB
model1 *= BB
model2 *= BB
'''
# Generate a number of epsilon models with different -beta values to visualize
# Apply models to wavelengths and plot
waves = np.logspace(0,3)*u.micron

model0 = GrayBody(waves, 500*u.K, beta=0)
model1 = GrayBody(waves, 500*u.K, beta=-1)
model2 = GrayBody(waves, 500*u.K, beta=-2)
plt.close()
with quantity_support():
    plt.figure()
    plt.loglog(waves, model0, label=r'$\beta$=0')
    plt.loglog(waves, model1, label=r'$\beta$=-1')
    plt.loglog(waves, model2, label=r'$\beta$=-2')
    plt.legend()






    
# For this demonstration, we generate an initial model with $\beta=-1.5$ and a blackbody temperature of 200 K.  The two fittable parameters for this model are then the resulting graybody temperature as well as the bolometric flux of the resulting spectral energy distribution.

# In[44]:
exit()

# Make initial graybody model with beta = -1.5
'''
GrayBody = PowerLaw1D(amplitude=1, x_0=53*u.micron, alpha=1.5,
                      fixed={'amplitude':True,'x_0':True,'alpha':True})
GrayBody *= BlackBody1D(temperature=200*u.K)
#GrayBody = BlackBody1D() * PowerLaw1D()
GrayBody.name = 'GrayBody'
'''

class GrayBody(PowerLaw1D * BlackBody1D):
    '''CompoundModel of an emissivity-modified BlackBody'''

gb_model = GrayBody(amplitude_0=1, x_0_0=53*u.micron, alpha_0=1.5,
                    temperature_1=200*u.K,
                    fixed={'amplitude_0':True,'x_0_0':True,'alpha_0':True})


'''
def unit_dict(input_units,output_units):
    return {'temperature_1':u.K,'bolometric_flux_1':FNU,
            'x_0_0':u.micron}
GrayBody._parameter_units_for_data_units = unit_dict
GrayBody = GrayBody.with_units_from_data(x=u.micron,
                                         y=FNU)
#GrayBody.return_units = u.Jy
#GrayBody._supports_unit_fitting = True
print(GrayBody)
print(GrayBody._supports_unit_fitting)
#exit()
'''


# In[52]:


# Fit graybody
from astropy.modeling.fitting import LevMarLSQFitter, SLSQPLSQFitter
from astropy.utils.console import ProgressBar
from scipy.optimize import curve_fit

#fitter = LevMarLSQFitter()
fitter = SLSQPLSQFitter()

waves = [53, 89, 155, 216]*u.micron

# dataiter holds 4-tuples of each pixel in all 4 wavelengths
dataiter = list(zip(a_orig.data.flat,c_repr.data.flat,d_repr.data.flat,e_repr.data.flat))

# This is the function that we apply to every pixel in the maps.
#   Note that we first convert the data values from Jy to erg/cm2/s to match the output from GrayBody
def fit_gb(y):
    y = (y*u.Jy).to(FNU, equivalencies=u.spectral_density(waves))
    print(y)
    print('wut')
    return fitter(GrayBody, waves, y, maxiter=50, acc=1e-5)

#print(GrayBody)
#print(gb_model(waves.value))
##fitter(gb_model, waves, [1,1,1,1], maxiter=50, acc=1e-5)
#print(fitter(GrayBody,waves,[1,1,1,1]))

#curve_fit(gb_model,waves,dataiter[1],p0=


exit()
#print(fit_gb([5,1,0.8,0.7]))
#results = ProgressBar.map(fit_gb, dataiter,
#                          ipython_widget=True,multiprocess=False)

###temps = [bb.temperature.value for bb in results]*u.K
#temps = temps.reshape(a_orig.data.shape)
#temps = fits.PrimaryHDU(temps,header=hawc_a_header)
#temps.writeto('temps_gray.fits',overwrite=True)


# In[13]:


def GreyBody(x, temperature, bolometric_flux, wave_0, beta):
    tau = (wave_0 / x)**beta
    
    bb = BlackBody1D(temperature,bolometric_flux)    

class GreyBody(Fittable1DModel):
    temperature = Parameter(default=200, min=0, unit=u.K)
    bolometric_flux = Parameter(default=1, unit=u.erg / u.cm ** 2 / u.s)
    wave_0 = Parameter(default=53, unit=u.micron, fixed=True)
    beta = Parameter(default=1, fixed=True)

    # We allow values without units to be passed when evaluating the model, and
    # in this case the input x values are assumed to be frequencies in Hz.
    _input_units_allow_dimensionless = True

    # We enable the spectral equivalency by default for the spectral axis
    input_units_equivalencies = {'x': u.spectral()}
    
    def evaluate(self, x, temperature, bolometric_flux, wave_0, beta):
        if isinstance(temperature, u.Quantity):
            temperature = temperature.to(u.K, equivalencies=u.temperature())
        else:
            temperature = u.Quantity(temperature, u.K)

        # We normalize the returned blackbody so that the integral would be
        # unity, and we then multiply by the bolometric flux. A normalized
        # blackbody has f_nu = pi * B_lambda / (sigma * T^4), which is what we
        # calculate here. We multiply by the bolometric
        # flux to get the normalization right.
        flambda = ((np.pi * u.sr * blackbody_nu(x, temperature) /
                    const.sigma_sb / temperature ** 4).to(1 / u.Hz) *
                    bolometric_flux)
        
        tau = (wave_0 / x)**beta

        fnu = np.exp(-tau) * fnu

        # If the bolometric_flux parameter has no unit, we should drop the 1/Hz
        # and return a unitless value. This occurs for instance during fitting,
        # since we drop the units temporarily.
        if hasattr(bolometric_flux, 'unit'):
            return flambda
        else:
            return flambda.value


# In[ ]:


from astropy.modeling.blackbody import blackbody_lambda
import astropy.units as u
from astropy.visualization import quantity_support

'''
def greybody_lambda(x, temperature, scale=1, wave_0=53*u.micron, beta=1.5):
    #Blackbody model modified to include an emissivity factor (beta)
    tau = (wave_0 / x)**beta
    
    bb = blackbody_lambda(x, temperature)
    gb = (1 - np.exp(-tau)) * bb
    return scale*gb
'''

def greybody_lambda(wave, temperature, scale=1, wave_0=53*u.micron, beta=-1.5):
    '''Blackbody model modified to include an emissivity factor (beta)'''
    bb = blackbody_lambda(wave, temperature)
    gb = (wave/wave_0)**beta * bb
    return scale*gb


# In[ ]:


waves = np.logspace(0,3)*u.micron

model0 = greybody_lambda(waves,temperature=500*u.K, beta=0)
model1 = greybody_lambda(waves,temperature=500*u.K, beta=1)
model2 = greybody_lambda(waves,temperature=500*u.K, beta=2)

with quantity_support():
    plt.figure()
    plt.loglog(waves,model0/np.max(model0),label=r'$\beta$=0')
    plt.loglog(waves,model1/np.max(model1),label=r'$\beta$=1')
    plt.loglog(waves,model2/np.max(model2),label=r'$\beta$=2')
    plt.legend()


# In[ ]:


# Fit GB_Ratio model to data
from scipy.optimize import brentq
from functools import partial
from astropy.utils.console import ProgressBar

waves = [53, 89]*u.micron
fit_temp = partial(GB_Ratio, x1=waves[0], x2=waves[1], beta=1.5)

def minimize(val,func=fit_temp):
    if val <= 0 or np.isnan(val):
        return np.nan
    try:
        temp = brentq(lambda T:val-func(T=T),a=10,b=1200,maxiter=50,xtol=1e-8)
        return temp
    except ValueError:
        return np.nan


# In[ ]:


# Define ratio of greybodies
def GB_Ratio(x1,x2, T, wave_0=53*u.micron, beta=1.5):
    '''Return ratios of greybodies (gb[t2]/gb[t1])'''
    gb1 = greybody_lambda(x1, T, wave_0, beta)
    gb2 = greybody_lambda(x2, T, wave_0, beta)
    return gb2/gb1


# In[ ]:


# Fit GB_Ratio model to data
from scipy.optimize import brentq
from functools import partial
from astropy.utils.console import ProgressBar

waves = [53, 89]*u.micron
fit_temp = partial(GB_Ratio, x1=waves[0], x2=waves[1], beta=1.5)

def minimize(val,func=fit_temp):
    if val <= 0 or np.isnan(val):
        return np.nan
    try:
        temp = brentq(lambda T:val-func(T=T),a=10,b=1200,maxiter=50,xtol=1e-8)
        return temp
    except ValueError:
        return np.nan


ratios = c_repr.data / a_orig.data

#results = ProgressBar.map(minimize, ratios.flat,
#                          ipython_widget=True,multiprocess=True)

#temps = np.array(results).reshape(a_orig.data.shape)
#temps = fits.PrimaryHDU(temps,header=hawc_a_header)
#temps.writeto('temps_AC.fits',overwrite=True)
temps = fits.open('temps_AC.fits')


# In[ ]:


fig = FITSFigure(temps)
fig.show_colorscale(cmap=cmap)

# Add colorbar
fig.add_colorbar()
fig.colorbar.set_axis_label_text('Temp (K)')


# In[ ]:


# Fit blackbody
from astropy.modeling.models import BlackBody1D
from astropy.modeling.blackbody import FNU
from astropy.modeling.fitting import LevMarLSQFitter
#import astropy.units as u
from astropy.utils.console import ProgressBar

fitter = LevMarLSQFitter()

#bb = BlackBody1D(temperature=200*u.K)

#dataiter = zip(a_orig.data.flat,c_repr.data.flat,d_repr.data.flat,e_repr.data.flat)
#waves = [53, 89, 155, 216]*u.micron
dataiter = zip(a_orig.data.flat,c_repr.data.flat)
waves = [53, 89]*u.micron
#dataiter = zip(d_repr.data.flatten()[3000:3010],e_repr.data.flatten()[3000:3010])
#waves = [155,216]*u.micron

#gb = GreyBody(beta=1.5)
#print(gb)

#def fit_bb(y):
#    y = (y*u.Jy).to(FNU, equivalencies=u.spectral_density(waves))
#    return fitter(BB,waves,y, maxiter=50, acc=1e-5)

#def fit_gb(y):
#    y *= u.Jy
    #return fitter(gb, waves, y, maxiter=50, acc=1e-5)

'''
from astropy.table import Table
t = Table.read('temps_2color.ecsv')
temps = np.array(t['temps']).reshape(a_orig.data.shape)
temps[temps==200] = np.nan
temps = fits.PrimaryHDU(temps,hawc_a_header)
fig = FITSFigure(temps)
fig.show_colorscale(cmap=cmap)
fig.add_colorbar()
'''


#fit_gb([1,0.8])

###results = ProgressBar.map(fit_gb,list(dataiter),ipython_widget=True,multiprocess=True,step=1000)

###temps = [bb.temperature.value for bb in results]*u.K
###temps = temps.reshape(a_orig.data.shape)
###temps = fits.PrimaryHDU(temps,header=hawc_a_header)
###temps.write('temps_AC.fits',overwrite=True)
#from astropy.table import Table,Column
#t = Table()
#t.add_column(Column(temps,name='temps'))
#t.pprint()
#t.write('temps_AC.ecsv',overwrite=True)
#'''


# In[ ]:


#from matplotlib import gridspec

#fig = plt.figure(figsize=(10,10))
#gs = gridspec.GridSpec(2,2)
#positions = [list(ax.get_position(fig).bounds) for ax in gs]

files = ['sofia_data/F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits',
         'sofia_data/F0484_HA_POL_7600018_HAWCHWPC_PMP_022-114.fits',
         'sofia_data/F0481_HA_POL_7600012_HAWDHWPD_PMP_050-083.fits',
         'sofia_data/F0485_HA_POL_7600019_HAWEHWPE_PMP_055-075.fits']
titles = ['A','C','D','E']
         
#files = [('sofia_data/F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits','HAWC A', pos[0]),
#         ('sofia_data/F0484_HA_POL_7600018_HAWCHWPC_PMP_022-114.fits', 'HAWC C', ),
#         ('sofia_data/F0481_HA_POL_7600012_HAWDHWPD_PMP_050-083.fits', 'HAWC D', gs[1,0]),
#         ('sofia_data/F0485_HA_POL_7600019_HAWEHWPE_PMP_055-075.fits', 'HAWC E', gs[1,1])]
#fig = plt.figure(figsize=(8,8))

#for file,title,pos in zip(files,titles,positions):
#    make_polmap(file,title,fig,subplot=pos)
#    f,t,g = tup
#    make_polmap(f,t,fig,g)


# In[ ]:


'''
from astropy.nddata import block_reduce, block_replicate

# Bin and reduce p by a factor of 2
block_size = 2
binned_p = block_reduce(p.data,block_size, func=np.nansum)

# Upsample back to original size
rebinned_p = block_replicate(binned_p,block_size)

# Apply the same quality cut mask from before, and mask any negative values
#rebinned_p[mask] = np.nan
rebinned_p[rebinned_p <= 0] = np.nan

# Convert HDU by copying header from the original p frame
rebinned_p = fits.PrimaryHDU(rebinned_p, header=p.header)
'''


# In[ ]:


#fig = FITSFigure(p)

# Show image
#fig.show_colorscale(cmap=cmap)

# Plot contours
#ncontours = 30
#fig.show_contour(stokes_i, colors='gray', levels=ncontours,
#                 smooth=1, kernel='box', linewidths=0.3)

# Add colorbar
#fig.add_colorbar()
#fig.colorbar.set_axis_label_text('$p$ (%)')


# In[ ]:


#from astropy.modeling import models, fitting

#binned_i = block_reduce(stokes_i.data, block_size, func=np.nanmean)

#plt.figure()
#plt.loglog(binned_i.flatten(),binned_p.flatten(),'ko')

#linear_mod = models.Linear1D()
#fitter = fitting.LevMarLSQFitter()

#fit = fitter(linear_mod,binned_i.flatten(), binned_p.flatten())
#print(fit)
#plt.show()

