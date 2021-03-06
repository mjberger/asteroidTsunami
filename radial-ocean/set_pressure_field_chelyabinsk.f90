!
subroutine set_pressure_field(maux,mbc,mx,my,xlow,ylow,dx,dy,time,aux)
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
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlow,ylow,dx,dy,time
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m,iint,jint
    real(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
    real(kind=8) :: yc, xc,dist,dist_in_km,pressRatio,peakRatio, fallOff             
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
    real(kind=8) :: blastx_center, blasty_center,currPress
    real(kind=8) :: airSpeed, dx_in_radians, dy_in_radians,maxRatio
    real(kind=8) :: overPressure,computedOverPressure, maxOverPressure, maxPress
    integer :: iloc,jloc
    real (kind=8) :: sumPress,xuse,yuse,wtx,wty   ! to do simpsons rule to get cell avg of pressure
    character(len=100) :: format_string
  
    airSpeed  = 344.0  ! geoclaw is dimensional, meters per second, roughly mach 1
    blastx_center = 0.
    blasty_center = 40.
    maxRatio = 1.
    maxRatio = 0.
    maxOverPressure = 0.
    maxPress = 0.0

    ! Set pressure field 
    aux(pressure_index, :, :) = ambient_pressure  !! units for geoclaw: amb is  101300 pascal (~101KPa)

    do j=1-mbc,my+mbc
       ym = ylow + (j - 1.d0) * dy
       y = ylow + (j - 0.5d0) * dy
       yp = ylow + real(j,kind=8) * dy
       do i=1-mbc,mx+mbc
          xm = xlow + (i - 1.d0) * dx
          x = xlow + (i - 0.5d0) * dx
          xp = xlow + real(i,kind=8) * dx

          ! prep for simpsons rule for conservative avg of pressure integral
          ! try to "see" it better on coarser grid
          sumPress = 0.
          do jloc = 1, 3
             wty = 1.
             if (jloc .eq. 2) wty = 4.
             yuse = ym + (jloc-1)*.5d0*dy
             do iloc = 1, 3
                wtx = 1.
                if (iloc .eq. 2) wtx = 4.
                xuse = xm + (iloc-1)*.5d0*dx
                !  convert dist in lat-long angle to meters
                dist_in_km = spherical_distance(xuse,yuse,blastx_center,blasty_center)/1000.

                overPressure = computedOverPressure(dist_in_km,time)  !overpress was percentage
                maxOverPressure = max(maxOverPressure,abs(overPressure))
                pressRatio = 1.0 + overPressure 
 
               sumPress = sumPress + wtx*wty*pressRatio
             end do
          end do

          sumPress = sumPress / 36.0   ! normalize for simpsons rule on rectangle
          maxRatio = max(maxRatio,sumPress)

          currPress = pressRatio * ambient_pressure ! will add decay in time next
          aux(pressure_index, i, j) = currPress
          maxPress = max(maxPress,currPress)

       enddo
    enddo

    format_string = "('time ',e12.5,' max pressure ',e15.7,'  max abs. val. overPressure ',e12.5)"
    write(*,format_string) time, maxPress, maxOverPressure

end subroutine set_pressure_field

! ==========================================================================

double precision function  computedOverPressure(dist_in_km,time)

!  # use model of overpressure (based on MJA calcs) as a function of x,y,time
   implicit none

   ! Arguments
   real(kind=8), intent(in) :: dist_in_km, time

   ! Locals
   real(kind=8) :: blast_radius, airSoundSpeed, peakRatio

   airSoundSpeed = 340.
   peakRatio = 1.4

   ! if shock wave travels at mach 1.4, and sound speed at sea level is 340m/s
   ! then at time t it travels 340*time/1000 km.
   blast_radius = 1.4d0*airSoundSpeed*time/1000.d0 

  ! if (dist_in_km <= blast_radius) then   
      ! no decay in time here, just constant pressure, at the peak value
      ! remove +1 at end of formula since is overpressure, not total pressure
       computedOverPressure  = 1./(1.+5.*(dist_in_km/50.)**(2.5) ) *(peakRatio-1.) 
   !else
   !    computedOverPressure = 0.0
   !endif

end function computedOverPressure
