c
c ----------------------------------------------------------------------
c
      subroutine makeXT(time,xstHov,ystHov,xendHov,yendHov,
     .                          nvar,naux,npts,ihovUnit)

      use amr_module
      use geoclaw_module, only: dry_tolerance

      implicit double precision (a-h,o-z)

!      dimension  hovLine(nvar+1,npts)  ! to include sealevel variable eta
      dimension  hovLine(2,npts),xpt(npts)  ! only output eta

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

c     make diagram using finest resolution
      dxfine = hxposs(mxnest)
      dyfine = hyposs(mxnest)
c     initialize line to "unset" val, 
      hovLine = -9999.

c     fill with interpolation using finest level spacing
c     overwritten with finer grids if exists
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
         if (xbegin .le. xlast .and. ybegin .le. ylast) then ! transect intersect this grid
            ! index for insertion onto transect
            iptx = (xbegin-xstHov)/dxfine + 1
            ipty = (ybegin-ystHov)/dyfine + 1
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
!           if (iptend .gt. npts) then
!              write(*,*)" still dont have correct index calc "
!           endif
           xhov = xbegin + .5d0 * dxhov
           yhov = ybegin + .5d0 * dyhov
           do ii  = ipt, iptend
c            compute i,j index for source grid, interp val insert into hovLine
             ix = (xhov-(xlo+.5d0*dx))/dx + nghost  ! compute index into grid for interp
             jy = (yhov-(ylo+.5d0*dy))/dy + nghost  ! grid has ghost cells
             ix = max( min(ix,mitot-nghost), nghost+1)  ! natch, first case is right 
             jy = max( min(jy,mjtot-nghost), nghost+1)  ! on the border.

             h = alloc(iadd(1,ix,jy))
             eta = h + alloc(iaddaux(1,ix,jy))

             etaInterp = binterp(alloc(loc),alloc(locaux),nvar,naux,
     .                          mitot,mjtot,ix,jy,mptr,xhov,yhov,
     .                          xlo,ylo,dx,dy,press)
!--             hovLine(1,ii) = h
!--             if (h > dry_tolerance) then
!--               hovLine(2,ii) = alloc(iadd(2,ix,jy))/h ! pw const for now
!--               hovLine(3,ii) = alloc(iadd(3,ix,jy))/h ! pw const for now
!--             else
!--                hovLine(2:3,ii) = 0.d0
!--             endif
!--             hovLine(4,ii) = eta

                hovLine(1,ii) = etaInterp
                hovLine(2,ii) = press

           xhov = xhov + dxhov ! advance location on line
           yhov = yhov + dyhov
           end do
         endif

         mptr = node(levelptr, mptr)
         go to 20  

 30      level = level + 1
         go to 10

 90     continue

c       output in tecplot format
        do ii = 1, npts
           xpt(ii) = xstHov+(ii-1)*dxFine
        end do

        write(ihovunit,100)(xpt(k),time,hovLine(1,k),hovLine(2,k),
     .                      k=1,npts)
 100    format(4e15.7)

        return
        end
c
c ----------------------------------------------------------------------------------
c
       double precision function binterp(q,aux,nvar,naux,mitot,mjtot,
     .                                   ix,jy,mptr,xhov,yhov,
     .                                   xlo,ylo,dx,dy,press)

       use amr_module
       use storm_module, only: pressure_index

       implicit double precision (a-h,o-z)      

       dimension q(nvar,mitot,mjtot), aux(naux,mitot,mjtot)

!     check that ix,jy is lower left corner enclosing xhov,yov
      xcorn = xlo+(ix-.5-nghost)*dx
      ycorn = ylo+(jy-.5-nghost)*dy

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
      if (xcorn .gt. xhov)  val = eta11
      if (ycorn .gt. yhov)  val = eta11
      if (xcorn+dx .lt. xhov) val = eta00
      if (ycorn+dy .lt. yhov) val = eta00

      if (xcorn   .le. xhov .and. ycorn    .le. yhov .and.
     .    xcorn+dx.ge. xhov .and. ycorn+dy .ge.yhov) then ! best case
 
        alfax = (xhov-xcorn)/dx 
        alfay = (yhov-ycorn)/dy

        val = (1.d0-alfax)*(1.d0-alfay)*eta00 + alfax*(1.d0-alfay)*eta10
     .         + alfay*(1.d0-alfax)*eta10      + alfax*alfay*eta11

        press = (1.d0-alfax)*(1.d0-alfay)*press00 + alfax*(1.d0-alfay)*press10
     .         + alfay*(1.d0-alfax)*press10      + alfax*alfay*press11
      endif

      binterp = val
          
      return
      end function
