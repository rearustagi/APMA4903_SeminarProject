# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 22:19:39 2021

@author: rearu
"""

import warnings

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText

path='C:/Users/rearu/OneDrive/Desktop/APMA_Seminar_Project/Data/Decadal/'
mask_file = 'C:/Users/rearu/OneDrive/Desktop/APMA_Seminar_Project/regmask_1117_anl.nc'
mask=xr.open_dataset(mask_file)

def map_plot(mask_ds, vals, b, t):
    lats = mask_ds.variables['lat'][:]
    lons = mask_ds.variables['lon'][:]

    ax = plt.axes(projection=ccrs.PlateCarree())
    
    plt.contourf(lons, lats, vals, levels=np.linspace(b[0],b[1],20),
                 transform=ccrs.PlateCarree(),cmap="BuPu", extend = "max")
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    states50 = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none',
            edgecolor='k')
    ax.add_feature(states50, zorder=1, linewidth=1)
    ax.set_extent([-87, -58, 37, 52], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', linestyle='--')
    gl.bottom_labels = False
    gl.right_labels = False
    
    max, min = np.round(vals.max().values, 2), np.round(vals.min().values, 2)
    mean = np.round(vals.mean().values, 2)
    #print(max, min)
    text = AnchoredText("min: {}, max: {}\nmean: {}".format(min, max, mean), loc=4, prop={'size': 12}, frameon=True)
    ax.add_artist(text)
    
    plt.colorbar(orientation="horizontal", label = '% of days')
    plt.tight_layout()
    plt.title(t, fontsize=20)
    plt.show()
    
def apply_mask(ds):
    layer = mask.regmask[0]
    mask_area = layer.where(layer == 10)
    mask_area = mask_area.where(mask_area != 10, 1)

    ne_usa_prec = ds.prec * mask_area.values
    ne_usa_prec = ne_usa_prec.where(ne_usa_prec >= 0)
    return ne_usa_prec, mask_area

def get_quantiles( ds, qt_values = (0.99, 0.999) ):
    qt_dims = ("lat", "lon", "time")
    ds_qt = ds.quantile(qt_values, dim=qt_dims)
    return qt_values, ds_qt.values

precs_1980s, m = apply_mask(xr.open_dataset(path+'1980s_prec.nc') / 24)
precs_2010s, m = apply_mask(xr.open_dataset(path+'2010s_prec.nc') / 24)

pct, prec_exts = get_quantiles(precs_1980s)
bounds = [[0,3],[0,0.4]]

t_steps = precs_1980s.time.count().values
compare = {"1980s" : precs_1980s, "2010s" : precs_2010s}


for p, e, b in zip(pct, prec_exts, bounds):
    for dec, data in compare.items():
        ext_count = data \
                    .where(data >= e).count(dim='time') / t_steps * 100 * m.values
        map_plot(mask, ext_count, b, \
                 "{} JJA, {}th Percentile = {} ".format(dec, p*100, round(e, 2))+r"$\frac{mm}{hr}$")