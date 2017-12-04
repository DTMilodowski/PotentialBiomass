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

ds,geoTrans1 = EO.load_NetCDF(NetCDF_file,lat_var = 'lat', lon_var = 'lon')
forestclass, geoTrans2, coord_sys = io.load_GeoTIFF_band_and_georeferencing(ForestCover_file)

# Calculate reforestation sequestration potential
cell_areas = ds.variables['areas']

nodata_mask = ds.variables['AGB_mean']==-9999

agb_ds = ds.variables['AGB_mean'][:]*cell_areas/10.**4 # convert carbon density to total carbon. Note conversion from m2 to ha
agb_ds[nodata_mask] = np.nan
agbpot_ds = ds.variables['AGBpot_mean'][:]*cell_areas/10.**4 # convert carbon density to total carbon. Converting from m2 to ha
agbpot_ds[nodata_mask] = np.nan

forests_ds = ds.variables['forests'][:]
seqpot_ds = agbpot_ds-agb_ds

# tidy up seqpot
seqpot_ds[forests_ds==1] = 0.
seqpot_ds[seqpot_ds<0] = 0.

# Forest Classes:
# 0   - Out of area study
# 1   - No change of primary degraded forest from 2000-2012
# 2   - No change of primary intact forest from 2000-2012
# 3   - No change of non-primary from 2000-2012
# 4   - Primary intact, cleared 2005
# 5   - Primary intact, cleared 2010
# 6   - Primary intact, cleared 2012
# 7   - Primary intact, degraded 2005
# 8   - Primary intact, degraded 2010
# 9   - Primary intact, degraded 2012
# 10  - Primary degraded, cleared 2005
# 11  - Primary degraded, cleared 2010
# 12  - Primary degraded, cleared 2012
# 13  - Primary intact degraded 2005, cleared 2010
# 14  - Primary intact degraded 2005, cleared 2012
# 15  - Primary intact degraded 2010, cleared 2012

forestclass[forestclass==0]=np.nan
classes = np.unique(forestclass[~np.isnan(forestclass)])
nclass = classes.size

agb = np.zeros(15)
agbpot = np.zeros(15)
seqpot = np.zeros(15)

for cc in range(0,nclass):
    mask = forestclass==cc+1
    agb[cc] = np.sum(agb_ds[mask])
    agbpot[cc] = np.sum(agb_ds[mask])
    seqpot[cc] = np.sum(seqpot_ds[mask])
    
    
