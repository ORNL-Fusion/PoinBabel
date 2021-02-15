      subroutine r8evbicub_GCvars(xget,yget,x,nx,y,ny,ilinx,iliny,
     >     fBR,fBPHI,fBZ,fER,fEPHI,fEZ,fgradBR,fgradBPHI,
     >     fgradBZ,fcurlBR,fcurlBPHI,fcurlBZ,inf2,ict,fvalBR,fvalBPHI,
     >     fvalBZ,fvalER,fvalEPHI,fvalEZ,fvalgradBR,fvalgradBPHI,
     >     fvalgradBZ,fvalcurlBR,fvalcurlBPHI,fvalcurlBZ,ier)
C
C  evaluate a 2d cubic Spline interpolant on a rectilinear
C  grid -- this is C2 in both directions.
C
C  this subroutine calls two subroutines:
C     herm2xy  -- find cell containing (xget,yget)
C     fvbicub  -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inf2
!============
      integer nx,ny                     ! grid sizes
      REAL*8 xget,yget                    ! target of this interpolation
      REAL*8 x(nx)                        ! ordered x grid
      REAL*8 y(ny)                        ! ordered y grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced
C
      REAL*8 fBR(0:3,inf2,ny)   ! function data
      REAL*8 fBPHI(0:3,inf2,ny)   ! function data
      REAL*8 fBZ(0:3,inf2,ny)   ! function data
      REAL*8 fER(0:3,inf2,ny)   ! function data
      REAL*8 fEPHI(0:3,inf2,ny)   ! function data
      REAL*8 fEZ(0:3,inf2,ny)   ! function data
      REAL*8 fgradBR(0:3,inf2,ny)   ! function data
      REAL*8 fgradBPHI(0:3,inf2,ny)   ! function data
      REAL*8 fgradBZ(0:3,inf2,ny)   ! function data
      REAL*8 fcurlBR(0:3,inf2,ny)   ! function data
      REAL*8 fcurlBPHI(0:3,inf2,ny) ! function data
      REAL*8 fcurlBZ(0:3,inf2,ny) ! function data
      
C
C       f 2nd dimension inf2 must be .ge. nx
C       contents of f:
C
C  f(0,i,j) = f @ x(i),y(j)
C  f(1,i,j) = d2f/dx2 @ x(i),y(j)
C  f(2,i,j) = d2f/dy2 @ x(i),y(j)
C  f(3,i,j) = d4f/dx2dy2 @ x(i),y(j)
C
C      (these are spline coefficients selected for continuous 2-
C      diffentiability, see mkbicub[w].for)
C
      integer ict(6)                    ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return d2f/dx2  (0, don't)
C  ict(5)=1 -- return d2f/dy2  (0, don't)
C  ict(6)=1 -- return d2f/dxdy (0, don't)
c                   the number of non zero values ict(1:6)
c                   determines the number of outputs...
c
c  new dmc December 2005 -- access to higher derivatives (even if not
c  continuous-- but can only go up to 3rd derivatives on any one coordinate.
c     if ict(1)=3 -- want 3rd derivatives
c          ict(2)=1 for d3f/dx3
c          ict(3)=1 for d3f/dx2dy
c          ict(4)=1 for d3f/dxdy2
c          ict(5)=1 for d3f/dy3
c               number of non-zero values ict(2:5) gives no. of outputs
c     if ict(1)=4 -- want 4th derivatives
c          ict(2)=1 for d4f/dx3dy
c          ict(3)=1 for d4f/dx2dy2
c          ict(4)=1 for d4f/dxdy3
c               number of non-zero values ict(2:4) gives no. of outputs
c     if ict(1)=5 -- want 5th derivatives
c          ict(2)=1 for d5f/dx3dy2
c          ict(3)=1 for d5f/dx2dy3
c               number of non-zero values ict(2:3) gives no. of outputs
c     if ict(1)=6 -- want 6th derivatives
c          d6f/dx3dy3 -- one value is returned.
C
C output arguments:
C =================
C
      REAL*8 fvalBR(*)          ! output data
      REAL*8 fvalBPHI(*)          ! output data
      REAL*8 fvalBZ(*)          ! output data
      REAL*8 fvalER(*)          ! output data
      REAL*8 fvalEPHI(*)          ! output data
      REAL*8 fvalEZ(*)          ! output data
      REAL*8 fvalgradBR(*)          ! output data
      REAL*8 fvalgradBPHI(*)          ! output data
      REAL*8 fvalgradBZ(*)          ! output data
      REAL*8 fvalcurlBR(*)          ! output data
      REAL*8 fvalcurlBPHI(*)          ! output data
      REAL*8 fvalcurlBZ(*)          ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the fourth output (depends on ict(...) spec)
C  fval(5) receives the fourth output (depends on ict(...) spec)
C  fval(6) receives the fourth output (depends on ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,0,0,1]
C   on output fval= [f,df/dx,df/dy,d2f/dxdy], elements 5 & 6 not referenced.
C
C    on input ict = [1,0,0,0,0,0]
C   on output fval= [f] ... elements 2 -- 6 never referenced.
C
C    on input ict = [0,0,0,1,1,0]
C   on output fval= [d2f/dx2,d2f/dy2] ... elements 3 -- 6 never referenced.
C
C    on input ict = [0,0,1,0,0,0]
C   on output fval= [df/dy] ... elements 2 -- 6 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i(1),j(1)                       ! cell indices
C
C  normalized displacement from (x(i),y(j)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C
      REAL*8 xparam(1),yparam(1)
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      REAL*8 hx(1),hy(1)
      REAL*8 hxi(1),hyi(1)
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C
C  ** the interface is very similar to herm2ev.for; can use herm2xy **
C---------------------------------------------------------------------
C
      i(1)=0
      j(1)=0
      call r8herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i(1),j(1),xparam(1),yparam(1),hx(1),hxi(1),hy(1),hyi(1),ier)
      if(ier.ne.0) return
c
      call r8fvbicub(ict,1,1,
     >     fvalBR,i,j,xparam,yparam,hx,hxi,hy,hyi,fBR,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalBPHI,i,j,xparam,yparam,hx,hxi,hy,hyi,fBPHI,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalBZ,i,j,xparam,yparam,hx,hxi,hy,hyi,fBZ,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalER,i,j,xparam,yparam,hx,hxi,hy,hyi,fER,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalEPHI,i,j,xparam,yparam,hx,hxi,hy,hyi,fEPHI,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalEZ,i,j,xparam,yparam,hx,hxi,hy,hyi,fEZ,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalgradBR,i,j,xparam,yparam,hx,hxi,hy,hyi,fgradBR,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalgradBPHI,i,j,xparam,yparam,hx,hxi,hy,hyi,fgradBPHI,inf2,
     >     ny)
      call r8fvbicub(ict,1,1,
     >     fvalgradBZ,i,j,xparam,yparam,hx,hxi,hy,hyi,fgradBZ,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalcurlBR,i,j,xparam,yparam,hx,hxi,hy,hyi,fcurlBR,inf2,ny)
      call r8fvbicub(ict,1,1,
     >     fvalcurlBPHI,i,j,xparam,yparam,hx,hxi,hy,hyi,fcurlBPHI,inf2,
     >     ny)
      call r8fvbicub(ict,1,1,
     >     fvalcurlBZ,i,j,xparam,yparam,hx,hxi,hy,hyi,fcurlBZ,inf2,ny)
C
      return
      end
