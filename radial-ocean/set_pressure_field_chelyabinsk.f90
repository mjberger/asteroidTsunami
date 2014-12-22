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
    real(kind=8) :: airSpeed, dx_in_radians, dy_in_radians,dsigma,pi,mindist,maxRatio
    integer :: iloc,jloc
    real (kind=8) :: sumPress,xuse,yuse,wtx,wty   ! to do simpsons rule to get cell avg of pressure
 
    airSpeed  = 344.0  ! geoclaw is dimensional, meters per second
    blastx_center = 0.
    blasty_center = 40.
    pi = 3.14159265357989
    mindist=1000000000.
    maxRatio = 1.

    ! Set pressure field  constnat in time for now
    aux(pressure_index, :, :) = ambient_pressure
!    if (time < 2000.d0) then
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
                   peakRatio = 1.036  !from popova, corr. to Chelyabinsk at 300kTNT, burst height 25km
                   pressRatio  = 1./(1.+5.*(dist_in_km/50.)**(2.5) ) *(peakRatio-1.) + 1.0
                   sumPress = sumPress + wtx*wty*pressRatio
                end do                   
                end do                   
                sumPress = sumPress / 36.0   ! normalize for simpsons rule on rectangle
                maxRatio = max(maxRatio,sumPress)

!!$                ! compute peak pressure at (x,y) then decay in time
!!$                if (dist>= airSpeed*time ) then  ! so refine within one cell on coarsest grid at time 0
!!$                   currPress = ambient_pressure
!!$                else

                 ! apply first as constant pressure once started see what it looks like
                  currPress = pressRatio * ambient_pressure ! will add decay in time next
!!$                endif

                aux(pressure_index, i, j) = currPress


             enddo
        enddo
!    end if
     !write(*,*)" max pressure ratio ",maxRatio
end subroutine set_pressure_field
