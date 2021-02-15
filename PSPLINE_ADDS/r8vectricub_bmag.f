      subroutine r8vectricub_bmag(ict,ivec,xvec,yvec,zvec,ivd,
     >     fvalBR,fvalBPHI,fvalBZ,
     >     nx,xpkg,ny,ypkg,nz,zpkg,
     >     fsplBR,fsplBPHI,fsplBZ,
     >     inf2,inf3,iwarn,ier)
c
c  vectorized spline evaluation routine -- 2d *compact* spline
c  1.  call vectorized zone lookup routine
c  2.  call vectorized spline evaluation routine
c
c--------------------------
c  input:
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iwarn1,iwarn2,iwarn3
!============
      integer ict(6)                    ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
c        ict(4)=1 for d2f/dx2 (don't evaluate if ict(4)=0)
c        ict(5)=1 for d2f/dy2 (don't evaluate if ict(5)=0)
c        ict(6)=1 for d2f/dxdy (don't evaluate if ict(6)=0)
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
c
c  list of (x,y) pairs:
c
      REAL*8 xvec(ivec)                   ! x-locations at which to evaluate
      REAL*8 yvec(ivec)         ! y-locations at which to evaluate
      REAL*8 Zvec(ivec)                   ! y-locations at which to evaluate
c
      integer ivd                       ! 1st dimension of output array
c
c    ivd -- 1st dimension of fval, .ge.ivec
c
c output:
      REAL*8 fvalBR(ivd,*)      ! output array
      REAL*8 fvalBPHI(ivd,*)    ! output array
      REAL*8 fvalBZ(ivd,*)      ! output array
c
c  fval(1:ivec,1) -- values as per 1st non-zero ict(...) element
c  fval(1:ivec,2) -- values as per 2nd non-zero ict(...) element
c   --etc--
c
c input:
      integer nx,ny,nz                     ! dimension of spline grids
      REAL*8 xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      REAL*8 ypkg(ny,4)         ! y grid "package" (cf genxpkg)
      REAL*8 Zpkg(ny,4)                   ! y grid "package" (cf genxpkg)
      integer inf2,inf3                      ! fspl 3rd array dimension, .ge.nx
      REAL*8 fsplBR(0:3,inf2,ny)  ! (compact) spline coefficients
      REAL*8 fsplBPHI(0:3,inf2,ny)! (compact) spline coefficients
      REAL*8 fsplBZ(0:3,inf2,ny)  ! (compact) spline coefficients      
c     
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
c
      integer ix(ivec)                  ! zone indices {j}
      REAL*8 dxn(ivec)                    ! normalized displacements w/in zones
      REAL*8 hx(ivec)                     ! h(j) vector
      REAL*8 hxi(ivec)                    ! 1/h(j) vector
c
      integer iy(ivec)                  ! zone indices {j}
      REAL*8 dyn(ivec)                    ! normalized displacements w/in zones
      REAL*8 hy(ivec)                     ! h(j) vector
      REAL*8 hyi(ivec)          ! 1/h(j) vector
c
      integer iz(ivec)                  ! zone indices {j}
      REAL*8 dzn(ivec)                    ! normalized displacements w/in zones
      REAL*8 hz(ivec)                     ! h(j) vector
      REAL*8 hzi(ivec)                    ! 1/h(j) vector
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
c
      if(nx.lt.2) then
         write(6,*) ' ?vectricub:  nx.lt.2:  nx = ',nx
         ier=1
      endif
c
      if(ny.lt.2) then
         write(6,*) ' ?vectricub:  ny.lt.2:  ny = ',ny
         ier=1
      endif
c
      if(nz.lt.2) then
         write(6,*) ' ?vectricub:  nz.lt.2:  nz = ',nz
         ier=1
      endif
c
      
      if(ivec.le.0) then
         write(6,*) ' ?vecbicub:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vecbicub:  output vector dimension less than input ',
     >      'vector dimension.'
         write(6,*) ' ivec=',ivec,' ivd=',ivd
         ier=1
      endif
c
      if(ier.ne.0) return
c
c  vectorized lookup
c
      ix=0
      iy=0
      iz=0
      call r8xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
      call r8xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
      call r8xlookup(ivec,zvec,nz,zpkg,2,iz,dzn,hz,hzi,iwarn3)
      iwarn=iwarn1+iwarn2+iwarn3
c
c  vectorized evaluation
c
      call r8fvtricub(ict,ivec,ivd,fvalBR,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplBR,inf2,inf3,nz)
      call r8fvtricub(ict,ivec,ivd,fvalBPHI,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplBPHI,inf2,inf3,nz)
      call r8fvtricub(ict,ivec,ivd,fvalBZ,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplBZ,inf2,inf3,nz)

      return
      end
