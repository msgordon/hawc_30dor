'''
Make better colorbars
'''
from mpl_toolkits.axes_grid1 import make_axes_locatable

def add_colorbar(fig, ax, im, side='right',size='5%',pad=0.05,label=None):
    # just a wrapper for a series of matplotlib functions
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side,size=size,pad=pad)
    bar = fig.colorbar(im,cax=cax)
    if label:
        bar.set_label(label)
    return bar
