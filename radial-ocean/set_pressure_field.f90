subroutine set_pressure_field(maux,mbc,mx,my,xlow,ylow,dx,dy,time,aux)
!     ============================================
!
!     # set auxiliary array pressure field 
!     aux(pressure_index, i , j) = Pressure Field

    use amr_module, only: mcapa, xupper, yupper, xlower, ylower

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level

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
    real(kind=8) :: yc, xc,dist,pressRatio,peakRatio, fallOff             
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
     

    peakRatio = 1.04  !from popova report, corresponds to Chelyabinsk at 300kTNT,burst height 25km

    ! Set pressure field  constnat in time for now
    aux(pressure_index, :, :) = ambient_pressure
    if (time < 2000.d0) then
        do j=1-mbc,my+mbc
            ym = ylow + (j - 1.d0) * dy
            y = ylow + (j - 0.5d0) * dy
            yp = ylow + real(j,kind=8) * dy
            do i=1-mbc,mx+mbc
                xm = xlow + (i - 1.d0) * dx
                x = xlow + (i - 0.5d0) * dx
                xp = xlow + real(i,kind=8) * dx

                dist = sqrt(x**2 + (y-40.d0)**2) 
                if (dist < 1.d0) then
                    aux(pressure_index, i, j) = 2.d0 * aux(pressure_index, i, j) * (2000.d0 - time) / 2000.d0 * (1.d0 - dist)
                end if
             enddo
        enddo
    end if

end subroutine set_pressure_field
