
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import os

from clawpack.geoclaw.data import Rearth  # radius of earth
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata

from mapper import latlong

import clawpack.geoclaw.surge.plot as surgeplot
import clawpack.geoclaw.surge.data as surgedata


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data
    #plotdata.format = 'binary'
    plotdata.format = 'ascii'

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir,'claw.data'))
    amrdata = amrclaw.AmrclawInputData(clawdata)
    amrdata.read(os.path.join(plotdata.outdir,'amr.data'))
    physics = geodata.GeoClawData()
    physics.read(os.path.join(plotdata.outdir,'geoclaw.data'))
    surge_data = surgedata.SurgeData()
    surge_data.read(os.path.join(plotdata.outdir,'surge.data'))

    probdata = clawutil.ClawData()
    #probdata.read('setprob.data', force=True)
    #theta_island = probdata.theta_island


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    #-----------------------------------------
    # Figure for pcolor plot for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    #plotitem.plot_var = 0/1/2 or plot that entry into q instead of a function
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    #plotitem.pcolor_cmin = -1.0
    #plotitem.pcolor_cmax = 1.0
    plotitem.pcolor_cmin = -.005
    plotitem.pcolor_cmax = .005
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1
    #plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    # plotitem.pcolor_cmap = colormaps.all_white  # to print better in B&W
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0]
    #plotitem.patchedges_show = 1
    plotitem.patchedges_show = 0
    #plotaxes.xlimits = [-85,-55]
    #plotaxes.ylimits = [13,45]
    plotaxes.xlimits = [-75,-70]
    plotaxes.ylimits = [38,43]


    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = True #False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [-3500,-2500,-1500, -500, 0, 500]
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'dashed','linewidths':2,'colors' : 'red' }
    plotitem.amr_contour_show = [0,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    #-----------------------------------------
    # Figure for speeds
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Speeds', figno=10)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Speed'
    plotaxes.scaled = True

    #def fixup(current_data):
    #    import pylab
    #    addgauges(current_data)
    #    t = current_data.t
    #    t = t / 3600.  # hours
    #    pylab.title('Speed at %4.2f hours' % t, fontsize=20)
    #    pylab.xticks(fontsize=15)
    #    pylab.yticks(fontsize=15)
    #plotaxes.afteraxes = fixup

    def speed(current_data):
        from pylab import where,sqrt
        q = current_data.q
        h = q[0,:]
        dry_tol = 0.001
        u = q[1,:] # this is just to initialize
        v = q[2,:] # to correct size
        s = 0*q[2,:] # to correct size

        nq = len(q[1,:])
        [n,m] = h.shape
        for ii in range(0,n):
           for jj in range(0,m):
             if h[ii,jj] > dry_tol:
                u[ii,jj] = q[1,ii,jj]/h[ii,jj]
                v[ii,jj] = q[2,ii,jj]/h[ii,jj]
                s[ii,jj] = sqrt(u[ii,jj]* u[ii,jj] + v[ii,jj]*v[ii,jj])
             else:
                u[ii,jj] = 0.
                v[ii,jj] = 0.
                s[ii,jj] = 0
        #print("max of u = " + str(max(u)))
              
        #u = where(h>dry_tol, q[1,:]/h, 0.)
        #v = where(h>dry_tol, q[2,:]/h, 0.)
        #s = sqrt(u**2 + v**2)
        #s = sqrt(u*2+v*v) #try not dividing or using where
        return s

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = speed
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmap = \
           colormaps.make_colormap({0:[1,1,1],0.5:[0.5,0.5,1],1:[1,0.3,0.3]})
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = .01
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [-75,-70]
    plotaxes.ylimits = [38,43]


    #-----------------------------------------
    # Figure for zoom around NYC
    #-----------------------------------------
    zoomWanted = True
    if zoomWanted:
        plotfigure = plotdata.new_plotfigure(name='Zoom1', figno=7)
        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes('zoom on nyc')
        plotaxes.title = 'Surface elevation'
        plotaxes.scaled = True
        manhattan_island=-73.5
        xisland,yisland = latlong(1600e3, manhattan_island, 40., Rearth)
        #plotaxes.xlimits = [xisland-0.6, xisland+0.6]
        plotaxes.xlimits = [manhattan_island-1, manhattan_island+1]
        plotaxes.ylimits = [40.5,41.5]
        #plotaxes.afteraxes = addgauges
        def bigfont(current_data):
            import pylab
            t = current_data.t
            pylab.title("Surface at t = %8.1f" % t, fontsize=20)
            pylab.xticks(fontsize=15)
            pylab.yticks(fontsize=15)
        plotaxes.afteraxes = bigfont

        # Water
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.show = True
        #plotitem.plot_var = geoplot.surface
        plotitem.plot_var = geoplot.surface_or_depth
        plotitem.pcolor_cmap = geoplot.tsunami_colormap
        plotitem.pcolor_cmin = -1.0
        plotitem.pcolor_cmax = 1.0
        plotitem.add_colorbar = False
        plotitem.amr_celledges_show = [0]
        plotitem.patchedges_show = 1

        # Land
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.show = True
        plotitem.plot_var = geoplot.land
        # plotitem.pcolor_cmap = colormaps.all_white  # to print better in B&W
        plotitem.pcolor_cmap = geoplot.land_colors
        plotitem.pcolor_cmin = 0.0
        plotitem.pcolor_cmax = 100.0
        plotitem.add_colorbar = False
        plotitem.amr_celledges_show = [1,0,0,0,0]
        plotitem.patchedges_show = 1

        # contour lines:
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.show = True
        plotitem.plot_var = geoplot.surface
        plotitem.contour_levels = [-0.8, -0.4, 0.4, 0.8]
        plotitem.amr_contour_colors = ['k']  # color on each level
        plotitem.kwargs = {'linewidths':2}
        plotitem.amr_contour_show = [0,0,0,1,1]  
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0

        # add contour lines of bathy if desired:
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.show = False
        plotitem.plot_var = geoplot.topo
        plotitem.contour_levels = linspace(-1000,-1000,1)
        plotitem.amr_contour_colors = ['g']  # color on each level
        plotitem.kwargs = {'linestyles':'dashed','linewidths':2}
        plotitem.amr_contour_show = [0,0,1]  
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0000, 7000]
    plotaxes.ylimits = [-.05, .05]
    #plotaxes.xlimits = 'auto'
    #plotaxes.ylimits = 'auto'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo

    def gaugedpress(current_data):
        q = current_data.q
        #dpress = (q[4,:] - 101300)/101300
        dpress = q[4,:]  # already output as relative:  dp/amb_pr
        return dpress

    def gs(current_data):
        q = current_data.q
        # different than speed function because q is function of time, not
        # x,y at the gauges.
        from numpy import where, sqrt
        h = q[0,:]
        #print('shape of h ' +  str(h.shape))
        dry_tol = 0.001
        u = where(h>dry_tol, q[1,:]/h, 0.)
        v = where(h>dry_tol, q[2,:]/h, 0.)
        ssq = sqrt(u*u+v*v)
        #s = sqrt(u**2 + v**2)
        s = sqrt(ssq)
        return ssq
        
    #plotitem.plot_var = gaugetopo
    #plotitem.plotstyle = 'g-'

    # Plot relative delta pressure as red curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = gaugedpress
    plotitem.plotstyle = 'r-'

    # add speed to this plot since cant get new one going
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = gs
    plotitem.plotstyle = 'g-'

    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
#        legend(('surface','topography','dp'),loc='lower left')
        legend(('surface','dp','speed'),loc='upper right')
        plot(t, 0*t, 'k')

    plotaxes.afteraxes = add_zeroline

    # -------------------------------
    # Figure for speed at gauges
    #--------------------------------

#    plotfigure = plotdata.new_plotfigure(name='Speed at gauges', figno=310, \
#                    type='each_gauge')
#

#
#    # Set up for axes in this figure:
#    plotaxes = plotfigure.new_plotaxes()
#    plotaxes.xlimits = 'auto'
#    plotaxes.ylimits = 'auto'
#    #plotaxes.xlimits = [0000, 7000]
#    #plotaxes.ylimits = 'auto'
#    plotaxes.title = 'Speed'
#
#    # Plot surface as blue curve:
#    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
#    #plotitem.plot_var = speed
#    plotitem.plot_var = gs  #gauge_speed
#    plotitem.plotstyle = 'b-'
#
#    plotfigure.show = True 
#

    #-----------------------------------------
    # Figure for grids alone
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='grids', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.title = 'grids'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = [1,1,0]   
    #plotitem.amr_patchedges_show = [1]     
    plotitem.amr_patchedges_show = [0]     

    
    # Pressure field
    plotfigure = plotdata.new_plotfigure(name='Pressure')
    plotfigure.show = True
    
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [-85,-55]
    #plotaxes.ylimits = [13,45]
    plotaxes.xlimits = [-75,-70]
    plotaxes.ylimits = [38,43]
    plotaxes.title = "Pressure Field"
    # plotaxes.afteraxes = gulf_after_axes
    plotaxes.scaled = True
    
    #pressure_limits = [surge_data.ambient_pressure / 100.0, 
    #                   2.0 * surge_data.ambient_pressure / 100.0]
    pressure_limits = [.999*surge_data.ambient_pressure / 100.0, 
                       1.001 * surge_data.ambient_pressure / 100.0]
    #pressure_limits = [-.000001*surge_data.ambient_pressure, 
    #                   .000001 * surge_data.ambient_pressure]
    surgeplot.add_pressure(plotaxes, bounds=pressure_limits)
    surgeplot.add_land(plotaxes)

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
