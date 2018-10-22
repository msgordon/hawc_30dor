
# coding: utf-8

# # HAWC+ Data Cookbook / 30 Doradus Data Release
# -------------------
# In this `jupyter` cookbook, we will explore the [HAWC+](https://www.sofia.usra.edu/science/instruments/hawc) data cube and describe some of the basic analysis techniques involving imaging polarimetry data.
# 
# This cookbook follows the SOFIA press release of 30 Doradus observations: [SOFIA Reveals Never-Before-Seen Magnetic Field Details](https://www.sofia.usra.edu/multimedia/science-results-archive/sofia-reveals-never-seen-magnetic-field-details).
# 
# The Level 4 reduced data from this program has been released immediately to the public and is available on the [SOFIA Data Cycle System (DCS)](https://dcs.sofia.usra.edu/).  This notebook will guide the reader through downloading the 30 Doradus data with a walkthrough of basic analysis techniques with `python`.

# ## Downloading HAWC+ Data
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

# ## SOFIA Data Organization
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

# ## Data Structure
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
# <hr>

# ### Stokes I
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
pxscale = stokes_i.header['CDELT2']*3600  # map scale in arcsec/pix
disp_i = stokes_i.copy()
disp_i.data /= pxscale**2
axs = FITSFigure(disp_i, figure=fig)    # load HDU into aplpy figure
axs.show_colorscale(cmap=cmap)            # display the data with WCS projection and chosen colormap

# FORMATTING
axs.set_tick_labels_font(size='small')
axs.set_axis_labels_font(size='small')

# Add colorbar
axs.add_colorbar()
axs.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')

plt.savefig('figs/Stokes_I.pdf',dpi=300)
plt.savefig('figs/Stokes_I.png',dpi=300)

# ### Stokes Q and U
# Similarly, we can plot the Stokes Q and Stokes U images:

# In[3]:


stokes_q = hawc['STOKES Q']
stokes_u = hawc['STOKES U']

fig = plt.figure(figsize=(12.8,9.6))
axq = FITSFigure(stokes_q, subplot=(1,2,1),figure=fig)  # generate FITSFigure as subplot to have two axes together
axq.show_colorscale(cmap=cmap)               # show Q
axq.add_colorbar()
axq.colorbar.set_axis_label_text('Flux (Jy/pix)')

axu = FITSFigure(stokes_u, subplot=(1,2,2),
                 figure=plt.gcf())
axu.show_colorscale(cmap=cmap)               # show U
axu.add_colorbar()
axu.colorbar.set_axis_label_text('Flux (Jy/pix)')

# FORMATTING
axq.set_title('Stokes Q')
axu.set_title('Stokes U')
axu.tick_labels.hide_y()
axu.axis_labels.hide_y()
#axq.ticks.set_xspacing(0.05)
#axu.ticks.set_xspacing(0.05)
axq.set_tick_labels_font(size='small')
axq.set_axis_labels_font(size='small')
axu.set_tick_labels_font(size='small')
axu.set_axis_labels_font(size='small')

plt.savefig('figs/Stokes_Q_U.pdf',dpi=300)
plt.savefig('figs/Stokes_Q_U.png',dpi=300)

# <hr>
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

plt.savefig('figs/Stokes_Q_error.pdf',dpi=300)
plt.savefig('figs/Stokes_Q_error.png',dpi=300)

# ------
# ### Polarized Intensity $I_p$
# Level 4 HAWC+ additionally provides extensions with the polarization percentage ($p$), angle ($\theta$), and their associated errors ($\sigma$).
# 
# Percent polarization ($p$) and error ($\sigma_p$) are calculated as:
# 
# \begin{align}
#     p & = 100\sqrt{\left(\frac{Q}{I}\right)^2+\left(\frac{U}{I}\right)^2} \\
#     \sigma_p & = \frac{100}{I}\sqrt{\frac{1}{(Q^2+U^2)}\left[(Q\,\sigma_Q)^2+(U\,\sigma_U)^2+2QU\,\sigma_{QU}\right]+\left[\left(\frac{Q}{I}\right)^2+\left(\frac{U}{I}\right)^2\right]\sigma_I^2-2\frac{Q}{I}\sigma_{QI}-2\frac{U}{I}\sigma_{UI}}
# \end{align}
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
pxscale = stokes_i.header['CDELT2']*3600
stokes_ip.data /= pxscale**2
stokes_i = hawc['STOKES I'].copy()
pxscale = stokes_i.header['CDELT2']*3600  # map scale in arcsec/pix
stokes_i.data /= pxscale**2

fig = plt.figure(figsize=(12.8,9.6))
axi = FITSFigure(stokes_i, subplot=(1,2,1),figure=fig)
axi.show_colorscale(cmap=cmap)               # show I
axi.add_colorbar()
axi.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')

axp = FITSFigure(stokes_ip, subplot=(1,2,2), figure=plt.gcf())
axp.show_colorscale(cmap=cmap)               # show Ip
axp.add_colorbar()
axp.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')


# FORMATTING
axi.set_title(r'$I$')
axp.set_title(r'$I_{p^\prime}$')
axp.tick_labels.hide_y()
axp.axis_labels.hide_y()
axi.set_tick_labels_font(size='small')
axi.set_axis_labels_font(size='small')
axp.set_tick_labels_font(size='small')
axp.set_axis_labels_font(size='small')

plt.savefig('figs/Stokes_Ip.pdf',dpi=300)
plt.savefig('figs/Stokes_Ip.png',dpi=300)


# ----
# ## Plotting Polarization Vectors
# 
# From the $Q$ and $U$ maps, the polarization angle $\theta$ is calculated in the standard way:
# 
# $\theta = \frac{90}{\pi}\,\mathrm{tan}^{-1}\left(\frac{U}{Q}\right)$
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
        fig.set_title(title,fontsize=14)

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


stokes_i, p, mask, fig = make_polmap(afile)
plt.savefig('figs/A_polmap.png',dpi=300)
plt.savefig('figs/A_polmap.pdf',dpi=300)

# ----
# ### Plotting Polarization Fraction
# 
# We can also plot the polarization fraction $p$ to better visualize the structure of 30 Doradus.  We plot the same contours from total intensity $I$ in the background.

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

plt.savefig('figs/A_pmap.pdf',dpi=300)
plt.savefig('figs/A_pmap.png',dpi=300)

# ----
# ## HAWC+ Polarization Maps
# 
# Finally, using the function defined above, we plot all four HAWC+ observations of 30 Doradus.

# In[9]:


files = [afile,cfile,dfile,efile]
titles = ['A (53 $\mu m$)',
          'C (89 $\mu m$)',
          'D (154 $\mu m$)',
          'E (214 $\mu m$)']

rows = []
for file, title in zip(files,titles):
    stokes_i,p,mask,fig = make_polmap(file,title)
    plt.savefig('figs/%s_polmap.pdf'%title[0],dpi=300)
    plt.savefig('figs/%s_polmap.png'%title[0],dpi=300)

    rows.append((stokes_i,p,mask,fig))


ncontours = [30,30,30,30]#[30,20,20,20]

for row,title,cont in zip(rows,titles,ncontours):
    stokes_i,p,_,_ = row

    fig = FITSFigure(p)
    fig.show_colorscale(cmap=cmap)
    fig.show_contour(stokes_i,colors='gray',levels=cont,
                     smooth=1,kernel='box',linewidths=0.3)
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('$p^\prime$ (%)')
    fig.set_title(title,fontsize=14)
    
    plt.savefig('figs/%s_pmap.pdf'%title[0],dpi=300)
    plt.savefig('figs/%s_pmap.png'%title[0],dpi=300)
