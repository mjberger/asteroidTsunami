"""
Plot the ocean and island topography from the files 
    ocean.topotype2
    island.topotype2
"""

from pylab import *

from clawpack.geoclaw import topotools

ocean = topotools.Topography()
ocean.read('ocean.topotype2', 2)

island = topotools.Topography()
island.read('island.topotype2', 2)

figure(figsize=(8,11))

ax = subplot(2,1,1)
#ocean.plot(axes=ax)
ocean.plot(axes=ax,contour_levels=[-4000,-3000,-2000,-1000,0],contour_kwargs={'colors':'r', 'linestyles':'dashed'})
ax.set_title('ocean topography')

ax = subplot(2,1,2)
island.plot(axes=ax, contour_levels=[-120,-90,-60,-30,0,30], \
        contour_kwargs={'colors':'k', 'linestyles':'solid'})
ax.set_title('zoom on island topography')

show()
fname = 'topography.png'
savefig(fname)
print "Created ",fname

