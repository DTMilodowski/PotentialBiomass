#==============================================================================#
# sankey.py                                                                    #
#------------------------------------------------------------------------------#
# Code to produce swanky sankey diagrams                                       #
# Inputs:                                                                      #
# - axis = a matplotlib axis into which the sankey diagram will be plotted     #
# - A = array (dimensions N x T) containing the data to be sankey-tified.      #
#   Should be presented as the class for each data element in N for timestep   #
#   in T                                                                       #
# - (optional) date_time = 1D array (dimensions T)  containing a date or time  #
#   stamp to distribute sankey bars according to time series. If blank, evenly #
#   distributed with no time info                                              #
# - (optional) colours = 1D array (dimensions A)  containing colours to be     #
#   used for associated classes. If not provided, uses colour ramp             #
# - (optional) colourmap = string (colour map name). This picks colours from   #
#   specified colour ramp. Default is viridis if not specified and no colours  #
#   given                                                                      #
# - (optional) patch_width_fraction = float (between 0 and 1) specifying       #
#   fraction of "inter-column" distance occupied by each column. Default is    #
#   0.1                                                                        #
#------------------------------------------------------------------------------#
# @author D.T. Milodowski                                                      #
# @date December 2017                                                          #
#==============================================================================#
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# Get perceptionally uniform colourmaps
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.set_cmap(cmaps.viridis)



# This is the plotting script
def plot_sankey(axis,A,x_locs = None,colours = None,colourmap = None,patch_width=None):
    
    colours_specified = True # will use colourmap if colours not specified
    if colours not in locals():
        colours_specified = False
        if colourmap in locals():
            cmap = cm.get_cmap(colourmap)
        else:
            cmap = cm.get_cmap('viridis') # default colour scheme
    
    # find classes
    classes = np.unique(A)
    n_points,n_steps = A.shape
    n_class = classes.size
    class_abundance = np.zeros((n_class,n_steps))

    # there are n_class**2 possible paths to consider for plotting
    paths_i = np.zeros((n_class**2,n_steps-1))
    paths_f = np.zeros((n_class**2,n_steps-1))

    # get colours for classes from cmap if required
    scale = np.arange(0.,n_class)
    scale /=scale.max()
    if colours_specified == False:
        colours = cmap(scale)

    if x_locs not in locals():
        x_locs = np.arange(n_steps)+1.
    if patch_width_fraction not in locals():
        path_width_fraction=0.1
    patch_width = path_width_fraction*(x_locs[0]-x_locs[-1]/float(n_steps)) # default patch width is 10% of average interval
        
    # loop through timesteps and fill in class abundances and path details
    for cc in range(0,n_class):
        for tt in range(0,n_steps):
            class_abundance[cc,tt] = np.sum(A[:,tt]==classes[cc])
        for cc2 in range(0,n_class):
            for tt in range(0,n_tsteps-1):
                paths_i[cc*n_class+cc2,tt] = np.sum(np.all((A[:,tt]==classes[cc],A[:,tt+1]==classes[cc2])))
                paths_f[cc2*n_class+cc,tt] = np.sum(np.all((A[:,tt+1]==classes[cc2],A[:,tt]==classes[cc])))

    # Now get the cumsum of the classes and paths so that we can plot the sankey
    # diagram
    class_abundance_cum = np.cumsum(class_abundance,axis=0)
    paths_i_cum = np.cumsum(paths_i,axis=0)
    paths_f_cum = np.cumsum(paths_f,axis=0)

    # Now we have everything that we need to plot the diagram
    # Loop through the number of steps
    patches = []
    for tt in range(0,n_steps):
        # Loop through each class, and create a polygon patch for each class
        class_llim = np.zeros(n_class)
        class_llim[1:] = class_abundance_cum[:-1,tt]
        class_ulim = class_abundance_cum[:,tt]

        x_l = x_locs[tt]-patch_width/2.
        x_r = x_locs[tt]+patch_width/2.

        for cc in range(0, n_class):
            UL = [x_l,class_ulim]
            UR = [x_r,class_ulim]
            LR = [x_r,class_llim]
            LL = [x_l,class_llim]
            box = np.array([UL,UR,LR,LL])
            polygon = Polygon(box,True,ec='0.5',fc=colours[cc])
            patches.append(polygon)
            
    p_class = PatchCollection(patches, match_original=True)
    axis.add_collection(p_class)
    
    return 0
