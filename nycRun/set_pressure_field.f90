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
    real (kind=8) :: overPressure, computedOverPressure, maxOverPressure, maxPress
    character(len=100) :: format_string
    real(kind=8),parameter :: pi= 3.14159265358979
  
    ! this is for 1 dimensional wave
    !blastx_center = -17.999

    ! for two dimensional radially symmetric pressure wave. give blast center in lat/long units of domain
    ! this is in shallow part near NYC
    !blasty_center = 40.00
    !blastx_center = -73.0

    ! this is in deep part to right of shelf
    !blasty_center =  38.9
    !blastx_center = -71.9

    ! this is in mid_atlantic
    blasty_center = 35.00
    blastx_center = -68.0

    maxRatio        = 0.
    maxOverPressure = 0.
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

                !overPressure = computedOverPressure(dist_in_km,time)/100. !  / 100 from percent to value
                overPressure = computedOverPressure(dist_in_km,time)       !  this version returns dp
                maxOverPressure = max(maxOverPressure,abs(overPressure))
                pressRatio = 1.0 + overPressure/ambient_pressure 

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

    !format_string = "('time ',e12.5,' mptr ',i3,' max pressure ',e15.7,'  max abs. val. overPressure ',e12.5)"
    !write(*,format_string) time, mptr, maxPress, maxOverPressure

end subroutine set_pressure_field

! ==========================================================================

double precision function  computedOverPressure(dist_in_km,time)

!  # use model of overpressure (based on MJA calcs) as a function of x,y,time

   implicit none

   ! Arguments
   real(kind=8), intent(in) :: dist_in_km, time

   ! Locals
!   real(kind=8),parameter :: ampl = 5., width = 5., tstar = .5;
   real(kind=8) :: ampl, width, tstar  ! they could be params, but for debugging make them vars
   real(kind=8) :: blast_radius,a, airSpeed

  !!! dont know how to work variable time in yet   
  !!! blast_radius = 2*sqrt(time); ! note that it slows down as time increases

   ! if shock wave travels at mach 1.4, and sound speed at sea level is 340m/s
   ! then at time t it travels 340*time/1000 km.
   airSpeed  = 344.0  ! geoclaw is dimensional, meters per second, roughly mach 1

   ! 1.4 or 3 below refer to multiple of mach number for shock speed
   blast_radius = 1.4d0*airSpeed*time/1000.d0 
   !blast_radius = 3.0d0*airSpeed*time/1000.d0 

   ampl  = 9.d0  !mikes model number
   width = 5.d0
   tstar = .5d0


   ! later will make it depend on x or y direction and blend, to allow
   ! for different propagation speeds. For now radially ymmetric
!   if (abs(dist_in_km) <= blast_radius) then
   !    computedOverPressure = 2*a*exp(-.8*(blast_radius - dist_in_km)/tstar) *   &
   !               (.5d0 - 1.1d0*(blast_radius -dist_in_km)/tstar) 
   !   computedOverPressure = 1013.d0  !this is  1% overpressure, ambient_p = 101300

       computedOverPressure = 10130.d0*exp(-(dist_in_km-blast_radius)**2)  ! 10% overpressure, but a pulse

!   else
!       computedOverPressure = 0.0
!   endif

end function computedOverPressure
