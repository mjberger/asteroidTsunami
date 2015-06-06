c
c ----------------------------------------------------------------------
c
      subroutine makeHovmoller(time,xstHov,ystHov,xendHov,yendHov,
     .                         nvar,naux,npts,ihovUnit)

      use amr_module
      use geoclaw_module, only: dry_tolerance

      implicit double precision (a-h,o-z)

!      dimension  hovLine(nvar+1,npts)  ! to include sealevel variable eta
      dimension  hovLine(npts)  ! only output eta

      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)


c     open for appending, need complete time history
c     ihovUnit = 44
c     if (time .eq. t0) then
c        write(*,*)" opening file and rewinding for hovmoller data"
c        open(unit=ihovUnit,file='hovSlice.data',status='unknown',
c    .         position='rewind',form='formatted')
c     else
c        write(*,*)" opening file for appending "
c        open(unit=ihovUnit,file='hovSlice.data',status='unknown',
c    .         position='append',form='formatted')
c     endif


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
     .                          xlo,ylo,dx,dy)
!--             hovLine(1,ii) = h
!--             if (h > dry_tolerance) then
!--               hovLine(2,ii) = alloc(iadd(2,ix,jy))/h ! pw const for now
!--               hovLine(3,ii) = alloc(iadd(3,ix,jy))/h ! pw const for now
!--             else
!--                hovLine(2:3,ii) = 0.d0
!--             endif
!--             hovLine(4,ii) = eta
                hovLine(ii) = etaInterp

           xhov = xhov + dxhov ! advance location on line
           yhov = yhov + dyhov
           end do
         endif

         mptr = node(levelptr, mptr)
         go to 20  

 30      level = level + 1
         go to 10

 90     continue
        if (npts > 4000) then
           write(*,*)" increase format statement:hov longer than 1 line"
        endif
!        write(ihovunit,100) (hovLine(4,ii),ii=1,npts)
        write(ihovunit,100) (hovLine(ii),ii=1,npts)
 100    format(4000e15.7)

        write(*,*) "hovmoller diagram with ",npts," points"
c       close(ihovunit)

        return
        end
