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

def plot_sankey(axis,A,date_time = None,colours = None):

    return 0
