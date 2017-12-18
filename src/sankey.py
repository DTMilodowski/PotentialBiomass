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
import sys
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
def plot_sankey_from_class_timeseries(axis,A, x_locs = np.array([]), colours = np.array([]), colourmap = 'viridis', patch_width_fraction=0.1):

    colours_specified = True # will use colourmap if colours not specified
    if colours.size==0:
        colours_specified = False
        cmap = cm.get_cmap(colourmap)
        
    # find classes
    classes = np.unique(A)
    n_points,n_steps = A.shape
    n_class = classes.size
    class_abundance = np.zeros((n_class,n_steps))

    # there are n_class**2 possible paths to consider for plotting
    paths = np.zeros((n_class,n_class,n_steps-1))

    # get colours for classes from cmap if required
    scale = np.arange(0.,n_class)
    scale /=scale.max()
    if colours_specified == False:
        colours = cmap(scale)

    if x_locs.size==0:
        x_locs = np.arange(n_steps)+1.
    patch_width = patch_width_fraction*(x_locs[0]-x_locs[-1])/float(n_steps)
        
    # loop through timesteps and fill in class abundances
    for cc in range(0,n_class):
        for tt in range(0,n_steps):
            class_abundance[cc,tt] = np.sum(A[:,tt]==classes[cc])

    # Now get the cumsum of the classes
    class_ulim = np.cumsum(class_abundance,axis=0)
    class_llim = np.zeros(class_ulim.shape)
    class_llim[1:,:] = class_ulim[:-1,:]
            
    # Now deal with the path details
    for cc1 in range(0,n_class):
        for cc2 in range(0,n_class):
            for tt in range(0,n_steps-1):
                paths[cc1,cc2,tt] = np.sum(np.all((A[:,tt]==classes[cc1],A[:,tt+1]==classes[cc2]),axis=0))

    # get starting ulim and llim for the paths
    paths_vec = paths.reshape(n_class**2,n_steps-1)
    paths_i_ulim = np.cumsum(paths_vec,axis=0)
    paths_i_llim = paths_i_ulim-paths_vec

    # now need to find the end ulim and llim
    paths_f_increment = np.cumsum(paths,axis=0)
    paths_f_u=np.zeros(paths.shape)
    paths_f_l=np.zeros(paths.shape)
    
    for cc1 in range(0,n_class):
        for cc2 in range(0,n_class):
            paths_f_u[cc1,cc2,:] = class_llim[cc2,1:]+paths_f_increment[cc1,cc2,:]
            paths_f_l[cc1,cc2,:] = paths_f_u[cc1,cc2,:]-paths[cc1,cc2,:]
            
    paths_f_ulim = paths_f_u.reshape(n_class**2,n_steps-1)
    paths_f_llim = paths_f_l.reshape(n_class**2,n_steps-1)

    # Now we have everything that we need to plot the diagram
    plot_sankey(axis, class_ulim, class_llim, paths_i_ulim, paths_i_llim, paths_f_ulim, paths_f_llim, x_locs, colours, patch_width)

# This is the plotting script, adapted to now take in a matrix of fluxes between classes
# rather than the raw class maps. This enables more generic use.
# The input matrix F should have 3 dimensions: n_class*n_class*n_steps
# The rows correspond with the initial class at each timestep
# The cols correspond with the new class at each timestep
# The value of the matric represents the flux from one class type to another in the
# timestep.
# So for a three class system, with classes i,j,k, at each timestep F is:
#
#       i:i   i:j   i:k
# F = ( j:i   j:j   j:k )
#       k:i   j:k   k:k

def plot_sankey_generic(axis,F, x_locs = np.array([]), colours = np.array([]), colourmap = 'viridis', patch_width_fraction=0.1):
    # basic info
    n_class,temp,n_fluxes = F.shape
    n_steps = n_fluxes+1 # fences and fence posts
    
    colours_specified = True # will use colourmap if colours not specified
    if colours.size==0:
        colours_specified = False
        cmap = cm.get_cmap(colourmap)
        
    # get colours for classes from cmap if required
    scale = np.arange(0.,n_class)
    scale /=scale.max()
    if colours_specified == False:
        colours = cmap(scale)

    if x_locs.size==0:
        x_locs = np.arange(n_steps)+1.
    patch_width = patch_width_fraction*(x_locs[0]-x_locs[-1])/float(n_steps)
        
    # get class sizes for each timestep
    class_abundance = np.zeros((n_class,n_steps))
    class_abundance[:,:n_fluxes]=np.sum(F,axis=1)
    class_abundance[:,-1]=np.sum(F[:,-1,:],axis=1)

    # Now get the cumsum of the classes
    class_ulim = np.cumsum(class_abundance,axis=0)
    class_llim = np.zeros(class_ulim.shape)
    class_llim[1:,:] = class_ulim[:-1,:]
            
    # Now deal with the path details
    # there are n_class**2 possible paths to consider for plotting
    paths = F.copy()
    
    # get starting ulim and llim for the paths
    paths_vec = paths.reshape(n_class**2,n_steps-1)
    paths_i_ulim = np.cumsum(paths_vec,axis=0)
    paths_i_llim = paths_i_ulim-paths_vec

    # now need to find the end ulim and llim
    paths_f_increment = np.cumsum(paths,axis=0)
    paths_f_u=np.zeros(paths.shape)
    paths_f_l=np.zeros(paths.shape)
    
    for cc1 in range(0,n_class):
        for cc2 in range(0,n_class):
            paths_f_u[cc1,cc2,:] = class_llim[cc2,1:]+paths_f_increment[cc1,cc2,:]
            paths_f_l[cc1,cc2,:] = paths_f_u[cc1,cc2,:]-paths[cc1,cc2,:]
            
    paths_f_ulim = paths_f_u.reshape(n_class**2,n_steps-1)
    paths_f_llim = paths_f_l.reshape(n_class**2,n_steps-1)

    # Now we have everything that we need to plot the diagram
    plot_sankey(axis, class_ulim, class_llim, paths_i_ulim, paths_i_llim, paths_f_ulim, paths_f_llim, x_locs, colours, patch_width)

    
# This is the basic plotting script
def plot_sankey (axis, class_ulim, class_llim, paths_i_ulim, paths_i_llim, paths_f_ulim, paths_f_llim, x_locs, colours, patch_width, k=12):
    n_class,n_steps = class_ulim.shape
    # First we add the paths. These are going to be curved because it looks a lot nicer
    # Will use a logistic function, which has a parameters:
    # L = vertical difference between coordinates to be joined
    # x0 = midpoint
    # k = steepness of curve, which we start with as k=12 as default parameter as it
    # gives nice curves
    for tt in range(0,n_steps-1):
        t1 = x_locs[tt]
        t2 = x_locs[tt+1]

        # loop through the paths.
        x = np.arange(t1,t2+(t2-t1)/1000.,(t2-t1)/1000.)
        x0 = (t2+t1)/2.
        y_prime = 1./(1.+np.exp(-k*(x-x0)))
        # now plot connectors
        for cc in range(0,n_class):    
          for cc2 in range(0,n_class):
            UL = paths_i_ulim[cc*n_class+cc2,tt]
            LL = paths_i_llim[cc*n_class+cc2,tt]
            
            UR = paths_f_ulim[cc*n_class+cc2,tt]
            LR = paths_f_llim[cc*n_class+cc2,tt]

            L_u = UR-UL
            upper_limit = UL + L_u*y_prime
            
            L_l = LR-LL
            lower_limit = LL + L_l*y_prime
            
            # only print paths if flux > 0
            if UL - LL > 0:
                axis.fill_between(x,lower_limit,upper_limit,color='0.5',alpha=0.25)

    # Next we add the columns
    # Loop through the number of steps
    patches = []
    for tt in range(0,n_steps):
        # Loop through each class, and create a polygon patch for each class
        x_l = x_locs[tt]-patch_width/2.
        x_r = x_locs[tt]+patch_width/2.

        for cc in range(0, n_class):
            UL = [x_l,class_ulim[cc,tt]]
            UR = [x_r,class_ulim[cc,tt]]
            LR = [x_r,class_llim[cc,tt]]
            LL = [x_l,class_llim[cc,tt]]
            box = np.array([UL,UR,LR,LL])
            polygon = Polygon(box,True,ec='0.5',fc=colours[cc])
            patches.append(polygon)
    
    p_class = PatchCollection(patches, match_original=True)
    axis.add_collection(p_class)
    
