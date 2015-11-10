c
c  read output file only until reach grid at level 2, then close file
c  output again with only level 1 grids. need to adjust grid count in
c  fort.t file, but not in fort.q (determined experimentally)
c  this is to short circuit very slow python plotting program
c

       program valinAll

       implicit double precision (a-h, o-z)
       real(kind=8), allocatable :: qeta(:)
       integer frameNumSt,frameNumend,frameNum, fTmp
       character*10  fname1,fname2,fname3
       character*21  fname4

       frameNumSt = 32   ! starting frame number
       frameNumEnd = 116

       do fTmp = frameNumSt, frameNumEnd

       fname1 = 'fort.t'
       fname2 = 'fort.q'
       fname3 = 'fort.b'
       fname4 = 'level1Only/fort.b'

       
c      make file names for this frame
       frameNum = fTmp
       do ipos = 10,7,-1
          idigit = mod(frameNum,10)
          fname1(ipos:ipos) = char(ichar('0')  + idigit)
          fname2(ipos:ipos) = char(ichar('0')  + idigit)
          fname3(ipos:ipos) = char(ichar('0')  + idigit)
          fname4(11+ipos:11+ipos) = char(ichar('0')  + idigit)
          frameNum = frameNum / 10
       end do

       open(unit=1,file=fname1,status='old',form='formatted')
       open(unit=2,file=fname2,status='old',form='formatted')
       open(unit=3,file=fname3,status='old',access='stream')
       open(unit=4,file=fname4,status='unknown',access='stream')

       read(1,*) time
       read(1,*) meqn
       read(1,*) ngridNum
       read(1,*) naux
       read(1,*) ndim
       read(1,*) nghost


       ngrids = 0

       do 

          read(2,100,err=101) numGrid
 100      format(i5)
 101      continue   ! error due to ***** from too small field for Numgrid
          read(2,100) level
          read(2,100) mx
          read(2,100) my
          read(2,102) xlow
 102      format(e18.8)
          read(2,102) ylow
          read(2,102) dx
          read(2,102) dy
          read(2,*) 


          if (level == 1) then
             ngrids = ngrids + 1
             mitot = mx+2*nghost
             mjtot = my+2*nghost

             allocate(qeta(4*mitot*mjtot))
             read(3) qeta
             write(4)qeta
             deallocate(qeta)
          else
             close(4)
             rewind(1)
             write(1,102) time
             write(1,103) meqn
             write(1,103) ngrids
             write(1,103) naux
             write(1,103) ndim
             write(1,103) nghost
 103         format(i6)
             exit
          endif


      end do
      end do

      stop    
      end

