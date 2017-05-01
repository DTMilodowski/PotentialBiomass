# This is a relatively simple program that calculate raster statistics for regions denoted by a multipolygon shapefile. 
# This is widely useful, but in this instance I am using it to collate national-level statistics for potential biomass
# simulations undertaken by J. Exbrayat et al.

import numpy as np
import rasterstats as rs
