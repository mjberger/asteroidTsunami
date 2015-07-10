# encoding: utf-8
"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import datetime

import numpy as np

import clawpack.geoclaw.data as surge


#                           days   s/hour    hours/day            
days2seconds = lambda days: days * 60.0**2 * 24.0
seconds2days = lambda seconds: seconds / (60.0**2 * 24.0)

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data as clawdata

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = clawdata.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    #  this is whole region covered by bathy
    clawdata.lower[0] = -80.0      # west longitude
    clawdata.upper[0] = -55.0      # east longitude

    #clawdata.lower[1] = 13.0       # south latitude
    clawdata.lower[1] = 25.0       # south latitude
    clawdata.upper[1] = 45.0       # north latitude


    #clawdata.lower[0] = -75.0      # west longitude
    #clawdata.upper[0] = -69.0      # east longitude

    #clawdata.lower[1] = 37.0       # south latitude
    #clawdata.upper[1] = 43.0       # north latitude
 


    # Number of grid cells:
    degree_factor = 20
    clawdata.num_cells[0] = int(clawdata.upper[0] - clawdata.lower[0]) * degree_factor
    clawdata.num_cells[1] = int(clawdata.upper[1] - clawdata.lower[1]) * degree_factor

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3 + 1 + 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = True                # True to restart from prior results
    clawdata.restart_file = 'fort.chk01300'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    clawdata.tfinal = 1000000000  # seconds  # in case not set in output_1, used later for region times

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        #                 day     s/hour  hours/day
        
        #clawdata.tfinal = 1000  # seconds
        clawdata.tfinal = 14000  # seconds

        # Output occurrence per day, 24 = every hour, 4 = every 6 hours
        #recurrence = 70    
        #clawdata.num_output_times = int((clawdata.tfinal - clawdata.t0) 
        #                                    / recurrence)
        clawdata.num_output_times = 70

        clawdata.output_t0 = True  # output at initial (or restart) time?
        

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 100
        clawdata.total_steps = 2000
        clawdata.output_t0 = True

        
    clawdata.output_format = 'binary'      # 'ascii', 'binary'
    #clawdata.output_format = 'ascii'        # 'ascii' or 'netcdf' 

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'
    clawdata.output_aux_onlyonce = False    # output aux arrays only at t0



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 3



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True 

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = .10
    #clawdata.dt_initial = .50

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    # clawdata.cfl_desired = 0.75
    clawdata.cfl_desired = 0.1 
    clawdata.cfl_max = .2

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 2**16




    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity


    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

   # ---------------
    # Gauges:
    # ---------------
    gauges  = rundata.gaugedata.gauges 
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]

    gaugeno = 0
    #for d in [1570e3, 1590e3, 1610e3, 1630e3]:
        #gaugeno = gaugeno+1
        #x,y = latlong(d, theta_island, 40., Rearth)
        #gauges.append([gaugeno, x, y, 0., 1e10])

    # for manhattan test
    #gauges.append([1,  -74.008, 40.71, 0, 1e10])  # tip of manhattan
    #gauges.append([1,  -74.00, 40.71, 3000, 1e10])  # tip of manhattan using ref by 10

    gauges.append([1,  -74.01 , 40.705,    0, 1e10])  # tip of manhattan
    gauges.append([2,  -74.0, 40.5,      0,  1e10])
    gauges.append([3,  -74.05, 40.66,    0, 1e10])
    gauges.append([4,  -72.5, 40.0, 0,   1e10])

    #gauges.append([5,  -72.2, 39.8, 0,   1e10]) # old one near continental shelf
    gauges.append([5,  -64.8014, 32.4146, 0,   1e10])  # new one near Bermuda

    gauges.append([6,  -71.9, 39.2, 0,   1e10])
    gauges.append([7,  -71.5, 38.5, 0,   1e10])
    #gauges.append([8,  -70.5, 37.5, 0,   1e10])
    gauges.append([8,  -70.002, 37.9, 0,   1e10])  # move to be closer to blast center
    gauges.append([9,  -69.5, 36.5, 0,   1e10])
    gauges.append([10, -68.5, 35.5, 0,   1e10])

    gauges.append([11, -73.085, 40.7, 0,   1e10]) # long island 



    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 3

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [2700.1, 2900.100,3000.100]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 100

 


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata


    # max number of refinement levels:
    amrdata.amr_levels_max = 3

    # List of refinement ratios at each level (length at least mxnest-1)
    #amrdata.refinement_ratios_x = [10,4,4,6,16]
    #amrdata.refinement_ratios_y = [10,4,4,6,16]
    #amrdata.refinement_ratios_t = [10,4,4,6,16]
    amrdata.refinement_ratios_x = [12,8,4,6,16]
    amrdata.refinement_ratios_y = [12,8,4,6,16]
    amrdata.refinement_ratios_t = [12,8,4,6,16]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft','center','center',
                        'center', 'center']



    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.0  # Richardson tolerance

    amrdata.flag2refine = True
    # Note: in geoclaw the refinement tolerance is set as wave_tolerance below 
    # and flag2refine_tol is unused!
    #amrdata.flag2refine_tol = 0.5  # tolerance used in this routine


    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 1    


    # ---------------
    # Regions:
    # ---------------
    # == setregions.data values ==
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    #regions.append([2,7,3000,clawdata.tfinal,-74.5,-72.5,40.25,41.5])
    #nyc zoom
    #regions.append([3,7,9000,clawdata.tfinal,-74.5,-72.5,40.25,41.5])
    regions.append([2,2,0,clawdata.tfinal,-74.5,-72.5,40.25,41.5])

    #allow 2 levels on left of blast towards coastline
    #regions.append([1,2,0,clawdata.tfinal,-100.,-69.0,37.0,100.])
    #allow 2 levels all over
    regions.append([1,2,0,clawdata.tfinal,-100.,-50.0,20.0,100.])

    #allow 3 levels around pressure wave only
    regions.append([1,3,0,300,-100.,-69.0,37.0,100.])



    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py



    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print "*** Error, this rundata has no geodata attribute"
        raise AttributeError("Missing geodata attribute")
       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius =  6367500.0

    # == Forcing Options
    geo_data.coriolis_forcing = False   #True
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025 # Overridden below
    geo_data.friction_depth = 1e10

    # == Algorithm and Initial Conditions ==
    #geo_data.sea_level = 0.28  # Due to seasonal swelling of gulf
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 0.001

    # Refinement Criteria
    refine_data = rundata.refinement_data
    refine_data.wave_tolerance = .025 
    # refine_data.wave_tolerance = 0.5
    # refine_data.speed_tolerance = [0.25,0.5,1.0,2.0,3.0,4.0]
    # refine_data.speed_tolerance = [0.5,1.0,1.5,2.0,2.5,3.0]
    #refine_data.speed_tolerance = [1.0,2.0,3.0,4.0]
    refine_data.deep_depth = 1e6
    refine_data.max_level_deep = 3
    #refine_data.variable_dt_refinement_ratios = True
    refine_data.variable_dt_refinement_ratios = False

    # == settopo.data values ==
    topo_data = rundata.topo_data
    topo_data.topofiles = []
    # for topography, append lines of the form
    #   [topotype, minlevel, maxlevel, t1, t2, fname]
    # geodata.topofiles.append([3, 1, 3, rundata.clawdata.t0, 
    #                                    rundata.clawdata.tfinal, 
    #                                    '../bathy/atlantic_2min.tt3'])
    # topo_data.topofiles.append([3, 1, 3, rundata.clawdata.t0, 
    #                                    rundata.clawdata.tfinal, 
    #                                    '../bathy/atlantic_2min.tt3'])
    topo_data.topofiles.append([3, 1, 3, rundata.clawdata.t0, 
                                       rundata.clawdata.tfinal, 
                                       '../../bathy/atlantic_2min.tt3'])
    topo_data.topofiles.append([3, 1, 5, rundata.clawdata.t0, 
                                       rundata.clawdata.tfinal, 
                                       '../../bathy/newyork_3s.tt3'])

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == setfixedgrids.data values ==
    rundata.fixed_grid_data.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]
    
    # == fgmax.data values ==
    fgmax_files = rundata.fgmax_data.fgmax_files
    # for fixed grids append to this list names of any fgmax input files
    fgmax_files.append('fgmax_grid1.txt')
    #fgmax_files.append('fgmax_grid2.txt')
    rundata.fgmax_data.num_fgmax_val = 1  # Save depth only
    
    return rundata
    # end of function setgeo
    # ----------------------


def set_storm(rundata):

    data = rundata.stormdata

    # Physics parameters
    data.rho_air = 1.15
    data.ambient_pressure = 101.3e3 # Nominal atmos pressure

    # Source term controls - These are currently not respected
    data.wind_forcing = False
    data.drag_law = 1
    data.pressure_forcing = True
    
    # Source term algorithm parameters
    # data.wind_tolerance = 1e-4
    # data.pressure_tolerance = 1e-4 # Pressure source term tolerance

    # AMR parameters
    data.wind_refine = [20.0,40.0,60.0] # m/s
    #data.R_refine = [60.0e3,40e3,20e3]  # m
    data.R_refine = [0]  # m

    # Storm parameters
    data.storm_type = 0 # Type of storm
    #data.landfall = days2seconds(sandy_landfall.days) + sandy_landfall.seconds

    # Storm type 2 - Idealized storm track
    #data.storm_file = os.path.expandvars(os.path.join(os.getcwd(),'sandy.storm'))

    return data


def set_friction(rundata):

    data = rundata.frictiondata

    # Variable friction
    data.variable_friction = True

    # Region based friction
    # Entire domain
    data.friction_regions.append([rundata.clawdata.lower, 
                                  rundata.clawdata.upper,
                                  [np.infty,0.0,-np.infty],
                                  [0.050, 0.025]])

    return data


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()

    rundata.add_data(surge.SurgeData(),'stormdata')
    set_storm(rundata)
    rundata.add_data(surge.FrictionData(),'frictiondata')
    set_friction(rundata)

    rundata.write()
