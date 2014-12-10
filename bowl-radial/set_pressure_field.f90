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
    do j=1-mbc,my+mbc
        ym = ylow + (j - 1.d0) * dy
        y = ylow + (j - 0.5d0) * dy
        yp = ylow + real(j,kind=8) * dy
        do i=1-mbc,mx+mbc
            xm = xlow + (i - 1.d0) * dx
            x = xlow + (i - 0.5d0) * dx
            xp = xlow + real(i,kind=8) * dx

            !convert angles to kilometers, use approx
            ! 1 degree = 111 km, even for longitude
            

            if (coordinate_system == 2) then
                x = deg2rad * earth_radius**2 * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                y = ym * deg2rad
            endif
            

           ! use my fit to popova fig. S39 
            dist = sqrt(x**2 + y**2)  ! distance from center of domain at (0,0)
            fallOff =  1./(1. + 5.*dist**3)
            pressRatio = peakRatio * fallOff !1/(1. + 5.*dist**3)
            aux(pressure_index,i,j) = pressRatio * ambient_pressure ! Pressure field

         enddo
    enddo

end subroutine set_pressure_field
