!
subroutine set_pressure_field(maux,mbc,mx,my,xlow,ylow,dx,dy,time,aux,mptr)
!     ============================================
!
!     # set auxiliary array pressure field 
!     aux(pressure_index, i , j) = Pressure Field

    use amr_module, only: mcapa, xupper, yupper, xlower, ylower

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level
    use geoclaw_module, only: spherical_distance

    use storm_module, only: wind_index, pressure_index, set_storm_fields
    use storm_module, only: ambient_pressure

    use friction_module, only: friction_index, set_friction_field
    
    use topo_module
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: mbc,mx,my,maux,mptr
    real(kind=8), intent(in) :: xlow,ylow,dx,dy,time
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc) 
    
    ! Locals
    integer :: i,j,m,iloc,jloc
    real(kind=8) :: x,y,xm,ym,xp,yp,blastx_center,blasty_center
    real(kind=8) :: yc, xc,dist,dist_in_km,pressRatio, fallOff             
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
    real(kind=8) :: dx_in_radians, dy_in_radians,dsigma,maxRatio,currPress
    real (kind=8) :: sumPress,xuse,yuse,wtx,wty   ! to do simpsons rule to get cell avg of pressure
    real (kind=8) :: overPressure, computedOverPressure, maxOverPressure, maxPress,minOverPressure
    character(len=100) :: format_string
    real(kind=8),parameter :: pi= 3.14159265358979
  

    ! for two dimensional radially symmetric pressure wave. give blast center in lat/long units of domain
    ! this is in shallow part near NYC
    blasty_center = 17.45
    blastx_center = 117.38


    maxRatio        = 0.
    maxOverPressure = 0.
    minOverPressure = 10000000000.
    maxPress        = 0.0

    ! Set background pressure field , then overwrite
    aux(pressure_index, :, :) = ambient_pressure  !! units for geoclaw: amb is  101300 pascal (~101KPa)

    do j=1-mbc,my+mbc
       ym = ylow + (j - 1.d0) * dy
       y  = ylow + (j - 0.5d0) * dy
       yp = ylow + real(j,kind=8) * dy
       do i=1-mbc,mx+mbc
          xm = xlow + (i - 1.d0) * dx
          x  = xlow + (i - 0.5d0) * dx
          xp = xlow + real(i,kind=8) * dx

          ! prep for simpsons rule for conservative avg of pressure integral
          ! try to "see" it better on coarser grid. this uses 9 points in cell
          ! instead of midpoint eval.
          sumPress = 0.d0
          do jloc = 1, 3
             wty = 1.d0
             if (jloc .eq. 2) wty = 4.d0
             yuse = ym + (jloc-1)*.5d0*dy
             do iloc = 1, 3
                wtx = 1.d0
                if (iloc .eq. 2) wtx = 4.d0
                xuse = xm + (iloc-1)*.5d0*dx
                !  convert dist in lat-long angle to meters

                ! **** CHANGED to make 1D Test Problem *****
                !  for traveling wave test case, make into a 1d problem, only look at dist in x
                ! this is for one-dimensional pressure test (consant in y) 
                !dist_in_km = spherical_distance(xuse,yuse,blastx_center,yuse)/1000.
                
                ! this is for two-dimensional radially symmetric pressure wave
                dist_in_km = spherical_distance(xuse,yuse,blastx_center,blasty_center)/1000.

                overPressure = computedOverPressure(dist_in_km,time)/100. !  / 100 from percent to value
                !overPressure = computedOverPressure(dist_in_km,time)       !  this version returns dp

                maxOverPressure = max(maxOverPressure,abs(overPressure))
                minOverPressure = min(minOverPressure,abs(overPressure))
                pressRatio = 1.0 + overPressure
                !pressRatio = 1.0 + overPressure/ambient_pressure 

                sumPress = sumPress + wtx*wty*pressRatio
             end do
          end do

          sumPress = sumPress / 36.0   ! normalize for simpsons rule on rectangle
          maxRatio = max(maxRatio,sumPress)

          ! computed in terms of ratio to match with flowCart output done the line
          currPress                 = pressRatio * ambient_pressure 
          aux(pressure_index, i, j) = currPress
          maxPress                  = max(maxPress, currPress)
       enddo
    enddo

    format_string = "('time ',e12.5,' mptr ',i3,' max pressure ',e15.7,'  max abs. val. overPressure % ',e12.5), 'min ',e12.5"
    !write(22,format_string) time, mptr, maxPress, maxOverPressure, minOverPressure

end subroutine set_pressure_field

! ==========================================================================
! converted from MJA C program with parameters than model Chelyabinsk explosion:  520kt @ 29km alt
!
!        return overpressure (in percent of sea level STD ATM) as a function
!        of rad = radius from ground zero (km) and time in seconds
!        ...based on Friedlander wave, but scaled and tuned from actual runs.
!
double precision function  computedOverPressure(dist_in_km,time)

   implicit none

   ! Arguments
   real(kind=8), intent(in) :: dist_in_km, time

   ! local arguments - could be params, but for testing make them variables
   !double precision :: maxAmp = 6.d0    !/* ...max overpressre (%) = height of envelope */
   !double precision :: width  = 150.d0  !/* ...width(km) at which envelpe is 1/2 of max */
   !double precision :: thick  = 5.d0    ! /*...pulse width(km) where overpressure > 0   */
   !double precision :: speed = 0.3915d0 ! /*    ...(km/sec) speed of pulse on ground */
   double precision :: maxAmp, width, thick, speed
   double precision :: c, p_t, t, g, p, amplitude,tail

   ! Locals
   !real(kind=8),parameter :: ampl = 5., width = 5., tstar = .5;
    real(kind=8) :: ampl,  tstar  ! they could be params, but for debugging make them vars
    real(kind=8) :: blast_radius,a, airSpeed,rad

    maxAmp  = 800.d0 !from MJA sims  ! for Tunguska sized object from Scott Lawrence
    width   = 90.d0
    thick   = 12.d0 
    speed   = 0.3915d0

    t    =  time * speed
    rad = dist_in_km   ! different notation, before global edit

    ! ...pulse -- functional fit from various blast simulations 
    ! mja model computes %overpressure, so need to convert to units
    ! so 6% becomes .06*101300 , done above in calling program
    if ( dist_in_km <= t) then
       p_t  =  thick*2.d0
       c    =  width/2.35482d0
       g    =  maxAmp * exp(-15.d0*(rad/c)**2)  ! /* ...Gaussian envelope */
       tail = 100.d0/(1.d0+(.05*rad)**2)
       amplitude = g + tail
       computedOverPressure = amplitude*exp(-5.d0*(t-rad)/p_t)*(1.0d0 - 5.0d0*(t-rad)/p_t)
    else
       computedOverPressure = 0.d0
    endif

end function computedOverPressure
