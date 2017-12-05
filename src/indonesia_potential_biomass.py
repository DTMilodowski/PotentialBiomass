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
forestclass = forestclass[2500:,200:6000]
# Calculate reforestation sequestration potential
cell_areas = ds.variables['areas'][2500:,200:6000]

nodata_mask = ds.variables['AGB_mean'][2500:,200:6000]==-9999

agb_ds = ds.variables['AGB_mean'][2500:,200:6000]*cell_areas/10.**4 # convert carbon density to total carbon. Note conversion from m2 to ha
agb_ds[nodata_mask] = np.nan
agbpot_ds = ds.variables['AGBpot_mean'][2500:,200:6000]*cell_areas/10.**4 # convert carbon density to total carbon. Converting from m2 to ha
agbpot_ds[nodata_mask] = np.nan

agbdef_ds = agb_ds-agbpot_ds
agbratio_ds = agbdef_ds/agbpot_ds


forests_ds = ds.variables['forests'][2500:,200:6000]
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

# create land cover maps
forestclass2000 = forestclass.copy()
forestclass2005 = forestclass.copy()
forestclass2010 = forestclass.copy()
forestclass2012 = forestclass.copy()

# 2000
intact = np.any((forestclass==4,forestclass==5,forestclass==6,forestclass==7,forestclass==8,forestclass==9,forestclass==13,forestclass==14,forestclass==15),axis=0)
degraded = np.any((forestclass==10,forestclass==11,forestclass==12),axis=0)
forestclass2000[intact]=2
forestclass2000[degraded]=1

# 2005
cleared = np.any((forestclass==4,forestclass==10),axis=0)
intact = np.any((forestclass==5,forestclass==6,forestclass==8,forestclass==9,forestclass==15),axis=0)
degraded = np.any((forestclass==7,forestclass==11,forestclass==12,forestclass==13,forestclass==14),axis=0)
forestclass2005[intact]=2
forestclass2005[degraded]=1
forestclass2005[cleared]=3

# 2010
cleared = np.any((forestclass==4,forestclass==5,forestclass==10,forestclass==11,forestclass==13),axis=0)
intact = np.any((forestclass==6,forestclass==9),axis=0)
degraded = np.any((forestclass==7,forestclass==8,forestclass==12,forestclass==14,forestclass==15),axis=0)
forestclass2010[intact]=2
forestclass2010[degraded]=1
forestclass2010[cleared]=3

# 2012
cleared = np.any((forestclass==4,forestclass==5,forestclass==6,forestclass==10,forestclass==11,forestclass==12,forestclass==13,forestclass==14,forestclass==15),axis=0)
degraded = np.any((forestclass==7,forestclass==8,forestclass==9),axis=0)
forestclass2012[degraded]=1
forestclass2012[cleared]=3

print '====================================================================='
print '\tland cover class areas in 1000s of km'
print '---------------------------------------------------------------------'
print 'year,\t\tintact,\t\tdegraded,\tcleared,\ttotal'
print '2000,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(forestclass2000==2)/1000.,np.sum(forestclass2000==1)/1000.,np.sum(forestclass2000==3)/1000.,np.sum(~np.isnan(forestclass2000)/1000.))
print '2005,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(forestclass2005==2)/1000.,np.sum(forestclass2005==1)/1000.,np.sum(forestclass2005==3)/1000.,np.sum(~np.isnan(forestclass2005)/1000.))
print '2010,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(forestclass2010==2)/1000.,np.sum(forestclass2010==1)/1000.,np.sum(forestclass2010==3)/1000.,np.sum(~np.isnan(forestclass2010)/1000.))
print '2012,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(forestclass2012==2)/1000.,np.sum(forestclass2012==1)/1000.,np.sum(forestclass2012==3)/1000.,np.sum(~np.isnan(forestclass2012)/1000.))
print '====================================================================='

print '====================================================================='
print '\tpotential biomass within each class, in 10^6 Mg C'
print '---------------------------------------------------------------------'
print 'year,\t\tintact,\t\tdegraded,\tcleared,\ttotal'
print '2000,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(agbpot_ds[forestclass2000==2])/10.**6.,np.sum(agbpot_ds[forestclass2000==1])/10.**6.,np.sum(agbpot_ds[forestclass2000==3])/10.**6.,np.sum(agbpot_ds/10.**6))
print '2005,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(agbpot_ds[forestclass2005==2])/10.**6.,np.sum(agbpot_ds[forestclass2005==1])/10.**6.,np.sum(agbpot_ds[forestclass2005==3])/10.**6.,np.sum(agbpot_ds/10.**6.))
print '2010,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(agbpot_ds[forestclass2010==2])/10.**6.,np.sum(agbpot_ds[forestclass2010==1])/10.**6.,np.sum(agbpot_ds[forestclass2010==3])/10.**6.,np.sum(agbpot_ds/10.**6.))
print '2012,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(agbpot_ds[forestclass2012==2])/10.**6.,np.sum(agbpot_ds[forestclass2012==1])/10.**6.,np.sum(agbpot_ds[forestclass2012==3])/10.**6.,np.sum(agbpot_ds/10.**6.))
print '====================================================================='

print '====================================================================='
print ' AGB deficit within each class, in 10^6 Mg C (2005 only)'
print '---------------------------------------------------------------------'
print 'year,\t\tintact,\t\tdegraded,\tcleared,\ttotal'
print '2005,\t\t%.0f,\t\t%.0f,\t\t%.0f,\t\t%.0f' % (np.sum(agbdef_ds[forestclass2005==2])/10.**6.,np.sum(agbdef_ds[forestclass2005==1])/10.**6.,np.sum(agbdef_ds[forestclass2005==3])/10.**6.,np.sum(agbdef_ds/10.**6.))
print '====================================================================='

print '====================================================================='
print ' Average ratio of deficit to potential agb (2005 only)'
print '---------------------------------------------------------------------'
print 'year,\t\tintact,\t\tdegraded,\tcleared,\ttotal'
print '2005,\t\t%.3f,\t\t%.3f,\t\t%.3f,\t\t%.3f' % (np.mean(agbratio_ds[forestclass2012==2]),np.mean(agbratio_ds[forestclass2012==1]),np.mean(agbratio_ds[forestclass2012==3]),np.mean(agbratio_ds[~np.isnan(forestclass2012)]))
print '====================================================================='
