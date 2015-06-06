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
    real(kind=8) :: yc, xc,dist,dist_in_m,pressRatio, fallOff             
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
    real(kind=8) :: dx_in_radians, dy_in_radians,dsigma,maxRatio,currPress
    real (kind=8) :: sumPress,xuse,yuse,wtx,wty   ! to do simpsons rule to get cell avg of pressure
    real (kind=8) :: overPressure, computedOverPressure, maxOverPressure, maxPress
    character(len=100) :: format_string
    real(kind=8),parameter :: pi= 3.14159265358979
    logical  :: simpson
  
    ! this is for 1 dimensional wave
    !blastx_center = -17.999
    simpson = .false.

    ! for two dimensional radially symmetric pressure wave. give blast center in cartesian coords center of domain
    blasty_center = 0.00
    blastx_center = -24000.0


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

          if (simpson) then
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
                  !dist_in_m = sqrt(xuse*xuse + yuse*yuse)
                  !dist_in_m = abs(xuse)    ! for 1d problem now

                  !overPressure = computedOverPressure(dist_in_m,time)/100. !  / 100 from percent to value
                  overPressure = computedOverPressure(blastx_center,xuse,time)       !  this version returns dp
                  maxOverPressure = max(maxOverPressure,abs(overPressure))
                  pressRatio = 1.0 + overPressure/ambient_pressure 
  
                  sumPress = sumPress + wtx*wty*pressRatio
               end do
            end do

            sumPress = sumPress / 36.0   ! normalize for simpsons rule on rectangle
            maxRatio = max(maxRatio,sumPress)
          else  ! use point eval
            xuse = x
            dist_in_m = abs(xuse)
            overPressure = computedOverPressure(blastx_center,xuse,time) 

            maxOverPressure = max(maxOverPressure,abs(overPressure))
            pressRatio = 1.0 + overPressure/ambient_pressure 
          endif

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

double precision function  computedOverPressure(blastx_center,xuse,time)

!  # use model of overpressure (based on MJA calcs) as a function of x,y,time

   implicit none

   ! Arguments
   real(kind=8), intent(in) :: xuse, time, blastx_center

   ! Locals
!   real(kind=8),parameter :: ampl = 5., width = 5., tstar = .5;
   real(kind=8) :: ampl, width, tstar  ! they could be params, but for debugging make them vars
   real(kind=8) :: airSpeed,blastloc,dist

  !!! dont know how to work variable time in yet   
  !!! blast_radius = 2*sqrt(time); ! note that it slows down as time increases

   ! if shock wave travels at mach 1.4, and sound speed at sea level is 340m/s
   ! then at time t it travels 340*time/1000 km.
   !airSpeed  = 344.0  ! geoclaw is dimensional, meters per second, roughly mach 1
   
   !for a 1 km ocean, sqrt(gh) ~ sqrt(10000) ~ 100 m/sec.
   ! choose larger speed s for moving pressure wave to test s>c.
   airSpeed = 50

    !  this decays nicely by itself so no need to test for radius
    !if (abs(dist_in_m) <= blast_radius) then
       ! 20% overpressure, in a pulse,centered at 0, approximately .5 km wide
       !computedOverPressure = 2.d0*10130.d0*exp(dist_in_m**2)  
    !else
    !   computedOverPressure = 0.0
    !endif

    computedOverPressure = 0.d0
    ! if less than one km, compute gaussian pressure pulse
    blastloc = blastx_center + airSpeed*time   ! current center of pulse
    dist = abs(blastloc-xuse)
    if (dist .lt. 1000) then
       computedOverPressure = 10130.d0*exp(-(dist/400.d0)**4)
    endif

end function computedOverPressure
