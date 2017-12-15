#===============================================================================
# sankey.py
#-------------------------------------------------------------------------------
# Code to produce swanky sankey diagrams
# Inputs:
# - axis = a matplotlib axis into which the sankey diagram will be plotted
# - A = array (dimensions N x T) containing the data to be sankey-tified. Should
#   be presented as the class for each data element in N for timestep in T
# - (optional) date_time = 1D array (dimensions T)  containing a date or time
#   stamp to distribute sankey bars according to time series. If blank, evenly
#   distributed with no time info
# - (optional) colours = 1D array (dimensions A)  containing colours to be used
#   for associated classes. If not provided, you risk the wrath of a potentially
#   random colour generator.
#-------------------------------------------------------------------------------
# @author D.T. Milodowski
# @date December 2017
#===============================================================================
import numpy as np
from matplotlib import pyplot as plt

# Get perceptionally uniform colourmaps
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.set_cmap(cmaps.viridis)

# This is the plotting script
def plot_sankey(axis,A,date_time = None,colours = None):

    classes = np.unique(A)
    
    
    return 0