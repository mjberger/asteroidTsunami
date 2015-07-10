c
c ----------------------------------------------------------------------
c
      subroutine makeXT(time,xstHov,ystHov,xendHov,yendHov,
     .                          nvar,naux,npts,ihovUnit)

      use amr_module
      use geoclaw_module, only: dry_tolerance

      implicit double precision (a-h,o-z)

!     dimension  hovLine(nvar+1,npts)  ! to include sealevel variable eta
      dimension  hovLine(2,npts),pt(npts), val(2)  ! output eta and press

      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)

c
c     this version write unstructured data in format x,t, variables
c     for use with tecplot.   Idea is that in the future will 
c     only output the "best" info but might not always be available
c     on all grids. have to see if tecplot can handle multiple data 
c     and take finest, or I have to.
c     but for now, called at multiples of coarse grid time step so
c     don't have to deal with this

c     make diagram using finest resolution for line, even if dont
c     have it everywhere and have interpolate
      dxfine = hxposs(mxnest)
      dyfine = hyposs(mxnest)
c     initialize all points on line to "unset" val, 
c     so will notice if there is a problem
      hovLine = -9999.

c     fill with interpolation using finest level spacing
c     fill from coarsest level up, so that will be
c     overwriting with finer level grids if they exist
c     grids are cell-centered.  hovLine starts at xstHov (node centered)
      level = 1
 10    if (level .gt. lfine) go to 90
         dx = hxposs(level)
         dy = hyposs(level)
         mptr = lstart(level)
 20      if (mptr .eq. 0) go to 30      
         nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         loc     = node(store1, mptr)
         locaux  = node(storeaux,mptr)
         mitot   = nx + 2*nghost
         mjtot   = ny + 2*nghost
         xlo = rnode(cornxlo,mptr)
         ylo = rnode(cornylo,mptr)
         xhi = rnode(cornxhi,mptr)
         yhi = rnode(cornyhi,mptr)
         xbegin = max(xlo,xstHov)
         xlast  = min(xhi,xendHov)
         ybegin = max(ylo,ystHov)
         ylast  = min(yhi,yendHov)
         dxhov = 0.d0  
         dyhov = 0.d0

         if (xbegin .le. xlast .and. ybegin .le. ylast) then ! transect intersects this grid
            ! starting index for insertion onto transect. Counting on
            ! no exact intersections here. Perturb hovLine off cartesian grid
            ! add 1 to get first point INSIDE grid on hovline
            iptx = int((xbegin-xstHov)/dxfine) + 1
            ipty = int((ybegin-ystHov)/dyfine) + 1

            if (xstHov .eq. xendHov) then  ! line is vertical
              ipt = ipty
              iptend = (ylast - ystHov)/dyfine
              dyhov = dyfine
              if (ystHov+(ipty-1)*dyfine .gt. ybegin .or.  
     .            ysthov+ipty*dyfine .lt. ybegin)  then
                 write(*,*)" trouble with hovline y index"
              endif
            else
              ipt = iptx
              iptend = (xlast - xstHov)/dxfine
              dxhov = dxfine
              if (xstHov+(iptx-1)*dxfine .gt. xbegin .or.  
     .            xsthov+iptx*dxfine .lt. xbegin)  then
                 write(*,*)" trouble with hovline x index"
              endif
            endif

           xhov = xstHov + ipt*dxhov   ! one of these dhov = 0 since
           yhov = ystHov + ipt*dyhov   ! only handle vert or horiz line

           do ii  = ipt, iptend
c            compute i,j index for source grid, interp val insert into hovLine
             ix = (xhov-(xlo+.5d0*dx))/dx + nghost  ! compute index into grid for interp
             jy = (yhov-(ylo+.5d0*dy))/dy + nghost  ! grid has ghost cells
             ix = max( min(ix,mitot-nghost), nghost+1)  ! natch, first case was right 
             jy = max( min(jy,mjtot-nghost), nghost+1)  ! on the border.

             call binterp(alloc(loc),alloc(locaux),nvar,naux,
     .                    mitot,mjtot,ix,jy,mptr,xhov,yhov,
     .                    xlo,ylo,dx,dy,val)
             hovLine(1,ii) = val(1)  ! eta
             hovLine(2,ii) = val(2)  ! press

             xhov = xhov + dxhov ! advance location on line
             yhov = yhov + dyhov
           end do
         endif 

         mptr = node(levelptr, mptr)
         go to 20  

 30      level = level + 1
         go to 10

 90     continue

c       output for tecplotting (using mja or matlab tesselation)
        if (ystHov .eq. yendHov) then
           do k = 1, npts
              pt(k) = xstHov+(k-1)*dxfine
           end do
        else
           do k = 1, npts  
              pt(k) = ystHov+(k-1)*dyfine
           end do
        endif

        write(ihovunit,100)(pt(k),time,hovLine(1,k),hovLine(2,k),
     .                      k=1,npts)
 100    format(4e15.7)

        return
        end
c
c ----------------------------------------------------------------------------------
c
       subroutine binterp(q,aux,nvar,naux,mitot,mjtot,
     .                    ix,jy,mptr,xhov,yhov,
     .                    xlo,ylo,dx,dy,val)

       use amr_module
       use storm_module, only: pressure_index, ambient_pressure

       implicit double precision (a-h,o-z)      

       dimension q(nvar,mitot,mjtot), aux(naux,mitot,mjtot)
       dimension val(2)

!     check that ix,jy is lower left corner enclosing xhov,yov
      xcorn   = xlo+(ix+.5-nghost)*dx
      ycorn   = ylo+(jy+.5-nghost)*dy

      eta00 = q(1,ix,   jy)   + aux(1,ix,  jy)
      eta01 = q(1,ix,   jy+1) + aux(1,ix,  jy+1)
      eta10 = q(1,ix+1, jy)   + aux(1,ix+1,jy)
      eta11 = q(1,ix+1, jy+1) + aux(1,ix+1,jy+1)

      press00 = aux(pressure_index,ix,  jy)
      press01 = aux(pressure_index,ix,  jy+1)
      press10 = aux(pressure_index,ix+1,jy)
      press11 = aux(pressure_index,ix+1,jy+1)

!     handle degenerate cases of linesensor in border around first 
!     and last cells by using pw constant
      if (xcorn .gt. xhov .or. ycorn .gt. yhov) then  
         eta   = eta11
         press = press11
      endif
      if (xcorn+dx .lt. xhov .or. ycorn+dy .lt. yhov) then
         eta  = eta00
         press = press00
      endif

      if (xcorn   .le. xhov .and. ycorn    .le. yhov .and.
     .    xcorn+dx.ge. xhov .and. ycorn+dy .ge.yhov) then ! best case
 
        alfax = (xhov-xcorn)/dx 
        alfay = (yhov-ycorn)/dy

        eta = (1.d0-alfax)*(1.d0-alfay)*eta00 + alfax*(1.d0-alfay)*eta10
     .         + alfay*(1.d0-alfax)*eta01     + alfax*alfay*eta11

        press = (1.d0-alfax)*(1.d0-alfay)*press00 + 
     .                 alfax*(1.d0-alfay)*press10 +
     .                 alfay*(1.d0-alfax)*press01 +
     .                        alfax*alfay*press11
      endif

      val(1) = eta
      val(2) = press/ambient_pressure - 1.d0  ! to get ordinary sized numbers for overpressure
              
      return
      end 
