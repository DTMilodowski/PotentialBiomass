# This is a relatively simple program that calculate raster statistics for regions denoted by a multipolygon shapefile. 
# This is widely useful, but in this instance I am using it to collate national-level statistics for potential biomass
# simulations undertaken by J. Exbrayat et al.

import numpy as np
import sys

import rasterstats as rs
import fiona

sys.path.append('./src')
import geospatial_utility_tools as geo

import qgis.core as qgis
from qgis.analysis import QgsZonalStatistics

sys.path.append('/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/src')
import prepare_EOlab_layers as EO

# The data files and shapefile for calculating zonal stats
NetCDF_file = '/home/dmilodow/DataStore_GCEL/AGB/AGBregpot.nc'
# specify variables of interest - note this order will determine the order with which the stats are returned in the output file
vars = ['AGBobs','AGBpot','AGBreg']

# some other info
DATADIR = './data/'
REPORTDIR = './reports/'
national_boundaries = '/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/ne_50m_admin_0_tropical_countries_small_islands_removed.shp'

# Filename for output report
potAGB_report = REPORTDIR+'tropical_potAGB_0_25deg_national_summary.csv'

# first up - load the features (so that country names can be written into summary)
shapefile= fiona.open(national_boundaries)
features = list(shapefile)
N = len(features)

"""
# Alternative is to use qgis
#polygonLayer = qgis.QgsVectorLayer(national_boundaries, layer_name, 'ogr') 
polygonLayer = qgis.QgsVectorLayer(national_boundaries) 
"""

# Note that rasterstats does not weight observations so we calculate the mean values separately using total values and
# total area used in calculation.  That is, the input AGB arrays should be in units of Mg, not spatial densities i.e.:
# Mh/ha (or similar).
ds, geoTrans = EO.load_NetCDF(NetCDF_file,lat_var = 'lat', lon_var = 'lon')

# ordering of geoTrans [ XMinimum, DataResX, 0, YMinimum, 0, DataResY ]
rows, cols = ds.variables[vars[0]].shape
latitude = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3],geoTrans[5])
longitude =  np.arange(geoTrans[0],cols*geoTrans[1]+geoTrans[0],geoTrans[1])
areas = geo.calculate_cell_area_array(latitude,longitude, area_scalar = 1./10.**4)

# loop through the variables, multiplying by cell areas to give values in Mg
for vv in range(0,len(vars)):
    print vars[vv]
    file_prefix = DATADIR + 'tropics_' + vars[vv] + '_total'

    out_array = np.asarray(ds.variables[vars[vv]]) * areas
    out_array[np.asarray(ds.variables[vars[vv]])==-9999]=-9999  # not sure why np.asarray step is needed but gets the job done
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
areas_out = areas.copy()
areas_out[np.asarray(ds.variables[vars[0]])==-9999] = -9999
area_file_prefix = DATADIR + 'tropics_cell_areas'
EO.write_array_to_data_layer_GeoTiff(areas_out, geoTrans, area_file_prefix)


#------------------------------------------------------------------------------------------------------
# Now we can use the geotiffs we've just created with rasterstats to get the summaries that we need
# Let's start with the area variable before moving on
area_stats = rs.zonal_stats(national_boundaries, area_file_prefix+'.tif', stats="count sum")

"""
# QGIS alternative - not working at the moment!
area_stats = QgsZonalStatistics(polygonLayer, area_file_prefix+'.tif','zonal_stats_',1)
area_stats.calculateStatistics(None)
"""

# for other variables, want to be as flexible as possible, so I'll create a dictionary that can hold as many
# variables as required, according to what is available (and wanted) in the original netcdf file
out_stats = {}
for vv in range(0,len(vars)):
    raster = DATADIR + 'tropics_' + vars[vv] + '_total.tif'
    out_stats[vars[vv]] = rs.zonal_stats(national_boundaries, raster, stats="count sum")
    
# Write report to file
out = open(pot_AGB_report,'w')
# First set up header
out.write('Country, Area (Ha)')
for vv in range(0,len(vars)):
    out.write(', mean '+ vars[vv] + ', total ' + vars[vv])
out.write('\n')

# now loop through data and write each line of the file
for i in range(0,N):
    country = features[i]['properties']['name']
    print country
    out.write(country + ', ' + str(area_stats[i]['sum']))
    for vv in range(0,len(vars)):
        out.write(', ' + str(out_stats[vars[vv]][i]['sum']/area_stats[i]['sum']) + ', ' + str(out_stats[vars[vv]][i]['sum']))
    out.write('\n')
out.close()
