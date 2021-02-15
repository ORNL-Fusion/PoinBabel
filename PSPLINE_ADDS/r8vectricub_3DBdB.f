      subroutine r8vectricub_3DBdB(ict,ivec,xvec,yvec,zvec,ivd,
     >     fvalBR,fvalBPHI,fvalBZ,fvaldBRdR,fvaldBPHIdR,fvaldBZdR,
     >     fvaldBRdPHI,fvaldBPHIdPHI,fvaldBZdPHI,
     >     fvaldBRdZ,fvaldBPHIdZ,fvaldBZdZ,
     >     fvalER,fvalEPHI,fvalEZ,
     >     nx,xpkg,ny,ypkg,nz,zpkg,fsplBR,fsplBPHI,fsplBZ,
     >     fspldBRdR,fspldBPHIdR,fspldBZdR,
     >     fspldBRdPHI,fspldBPHIdPHI,fspldBZdPHI,
     >     fspldBRdZ,fspldBPHIdZ,fspldBZdZ,
     >     fsplER,fsplEPHI,fsplEZ,inf4,inf5,iwarn,ier)
c
c  vectorized spline evaluation routine -- 3d *compact* spline
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
! idecl:  explicitize implicit REAL declarations:
      REAL*8 stat
!============
      integer ict(10)           ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
c        ict(4)=1 for df/dz  (don't evaluate if ict(4)=0)
c        ict(5)=1 for d2f/dx2 (don't evaluate if ict(5)=0)
c        ict(6)=1 for d2f/dy2 (don't evaluate if ict(6)=0)
c        ict(7)=1 for d2f/dz2 (don't evaluate if ict(7)=0)
c        ict(8)=1 for d2f/dxdy (don't evaluate if ict(8)=0)
c        ict(9)=1 for d2f/dxdz (don't evaluate if ict(9)=0)
c        ict(10)=1 for d2f/dydz (don't evaluate if ict(10)=0)
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
c
c  list of (x,y,z) triples:
c
      REAL*8 xvec(ivec)                   ! x-locations at which to evaluate
      REAL*8 yvec(ivec)                   ! y-locations at which to evaluate
      REAL*8 zvec(ivec)                   ! z-locations at which to evaluate
c
      integer ivd                       ! 1st dimension of output array
c
c    ivd -- 1st dimension of fval, .ge.ivec
c
c output:
      REAL*8 fvalBR(ivd,*)      ! output array
      REAL*8 fvalBPHI(ivd,*)    ! output array
      REAL*8 fvalBZ(ivd,*)      ! output array
      REAL*8 fvaldBRdR(ivd,*)      ! output array
      REAL*8 fvaldBPHIdR(ivd,*)    ! output array
      REAL*8 fvaldBZdR(ivd,*)   ! output array
      REAL*8 fvaldBRdPHI(ivd,*)   ! output array
      REAL*8 fvaldBPHIdPHI(ivd,*)    ! output array
      REAL*8 fvaldBZdPHI(ivd,*)   ! output array
      REAL*8 fvaldBRdZ(ivd,*)   ! output array
      REAL*8 fvaldBPHIdZ(ivd,*)    ! output array
      REAL*8 fvaldBZdZ(ivd,*)      ! output array
      REAL*8 fvalER(ivd,*)      ! output array
      REAL*8 fvalEPHI(ivd,*)    ! output array
      REAL*8 fvalEZ(ivd,*)      ! output array

c
c  fval(1:ivec,1) -- values as per 1st non-zero ict(...) element
c  fval(1:ivec,2) -- values as per 2nd non-zero ict(...) element
c   --etc--
c
c input:
      integer nx,ny,nz                  ! dimension of spline grids
      REAL*8 xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      REAL*8 ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
      REAL*8 zpkg(nz,4)                   ! z grid "package" (cf genxpkg)
      integer inf4                      ! fspl 4th array dimension, .ge.nx
      integer inf5                      ! fspl 5th array dimension, .ge.ny
      REAL*8 fsplBR(8,inf4,inf5,nz)  ! (compact) spline coefficients
      REAL*8 fsplBPHI(8,inf4,inf5,nz)! (compact) spline coefficients
      REAL*8 fsplBZ(8,inf4,inf5,nz) ! (compact) spline coefficients
      REAL*8 fspldBRdR(8,inf4,inf5,nz) ! (compact) spline coefficients
      REAL*8 fspldBPHIdR(8,inf4,inf5,nz)! (compact) spline coefficients
      REAL*8 fspldBZdR(8,inf4,inf5,nz) ! (compact) spline coefficients
      REAL*8 fspldBRdPHI(8,inf4,inf5,nz) ! (compact) spline coefficients
      REAL*8 fspldBPHIdPHI(8,inf4,inf5,nz)! (compact) spline coefficients
      REAL*8 fspldBZdPHI(8,inf4,inf5,nz) ! (compact) spline coefficients
      REAL*8 fspldBRdZ(8,inf4,inf5,nz) ! (compact) spline coefficients
      REAL*8 fspldBPHIdZ(8,inf4,inf5,nz)! (compact) spline coefficients
      REAL*8 fspldBZdZ(8,inf4,inf5,nz)  ! (compact) spline coefficients
      REAL*8 fsplER(8,inf4,inf5,nz)  ! (compact) spline coefficients
      REAL*8 fsplEPHI(8,inf4,inf5,nz)! (compact) spline coefficients
      REAL*8 fsplEZ(8,inf4,inf5,nz) ! (compact) spline coefficients
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
c
      integer, dimension(:), allocatable :: ix,iy,iz
      REAL*8, dimension(:), allocatable :: dxn,dyn,dzn
      REAL*8, dimension(:), allocatable :: hx,hxi,hy,hyi,hz,hzi
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
         write(6,*) ' ?vectricub:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vectricub:  output vector dimension less than input ',
     >      'vector dimension.'
         write(6,*) ' ivec=',ivec,' ivd=',ivd
         ier=1
      endif
c
      if(ier.ne.0) return
c
      allocate(ix(ivec), iy(ivec), iz(ivec),
     >   dxn(ivec), dyn(ivec), dzn(ivec),
     >   hx(ivec),  hy(ivec),  hz(ivec),
     >   hxi(ivec), hyi(ivec), hzi(ivec), stat=ier)
c
      if(ier.ne.0) then
         write(6,*)
     >      ' ?vectricub: memory allocation failure.'
         ier=99
      endif
c
      if(ier.ne.0) return
c
c  vectorized lookup
c
      ix=0; iy=0; iz=0
      call r8xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
      call r8xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
      call r8xlookup(ivec,zvec,nz,zpkg,2,iz,dzn,hz,hzi,iwarn3)
      iwarn=max(iwarn1,iwarn2,iwarn3)
c
c  vectorized evaluation
c
      call r8fvtricub(ict,ivec,ivd,fvalBR,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplBR,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvalBPHI,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplBPHI,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvalBZ,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplBZ,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBRdR,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBRdR,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBPHIdR,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBPHIdR,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBZdR,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBZdR,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBRdPHI,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBRdPHI,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBPHIdPHI,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBPHIdPHI,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBZdPHI,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBZdPHI,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBRdZ,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBRdZ,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBPHIdZ,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBPHIdZ,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvaldBZdZ,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fspldBZdZ,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvalER,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplER,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvalEPHI,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplEPHI,inf4,inf5,nz)
      call r8fvtricub(ict,ivec,ivd,fvalEZ,ix,iy,iz,dxn,dyn,dzn,
     >     hx,hxi,hy,hyi,hz,hzi,fsplEZ,inf4,inf5,nz)

c
      deallocate(ix,iy,iz,dxn,dyn,dzn,hx,hy,hz,hxi,hyi,hzi)
c
      return
      end
