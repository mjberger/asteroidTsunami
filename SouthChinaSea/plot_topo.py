"""
Plot bathymetry

"""

import matplotlib.pyplot as plt
import numpy
import clawpack.geoclaw.topotools as tt

topo = tt.Topography(path='./path_to_bathy_file')
topo.plot()
plt.show()
