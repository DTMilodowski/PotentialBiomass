import numpy as np
import os
import sys

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

sys.path.append('/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/src')
import prepare_EOlab_layers as EO

sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomass/src')
import geospatial_utility_tools as geo

sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/')
import data_io as io

# Get perceptionally uniform colourmaps
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.set_cmap(cmaps.viridis)

#plt.figure(1, facecolor='White',figsize=[2, 1])
#plt.show()

SAVEDIR = '/home/dmilodow/DataStore_DTM/'
NetCDF_file = '/disk/scratch/local.2/southeast_asia_PFB/southeast_asia_PFB_mean_WorldClim2.nc'
ForestCover_file = '/home/dmilodow/DataStore_DTM/FOREST2020/PartnerCountries/Indonesia/ForestCover/Primary_and_intact_forest_and_loss/change/margono_indonesia_forestcover_and_loss_1km.tif'

agb_ds,geoTrans1 = EO.load_NetCDF(NetCDF_file,lat_var = 'lat', lon_var = 'lon')
forestclass, geoTrans2, coord_sys = io.load_GeoTIFF_band_and_georeferencing(File)
