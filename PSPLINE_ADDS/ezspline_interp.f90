!/////
! R8 !
!/////

!!!
!!! 1-d
!!!

subroutine EZspline_interp1_r8(spline_o, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  real(ezspline_r8) p1  ! the location where the interpolation is sought
  real(ezspline_r8) f   ! the interpolation
  integer, intent(out) :: ier

  integer ifail
  integer, parameter :: ict(3)=(/1, 0, 0/)
  real(ezspline_r8) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isLinear == 1) then

     call r8pc1ev(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call r8evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, ansr, ifail)

  else

     call r8herm1ev(p1, &
          &         spline_o%x1(1), spline_o%n1, &
          &         spline_o%ilin1,&
          &         spline_o%fspl(1,1),  &
          &         ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_r8


subroutine EZspline_interp1_array_r8(spline_o, k, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k) ! location arrays
  real(ezspline_r8), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: i, ifail
  integer, parameter :: ict(3)=(/1,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isLinear == 1) then

     call r8vecpc1(ict, k, p1, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then

     call r8vecspline(ict, k, p1, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else

     call r8vecherm1(ict, k, p1, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_array_r8

!!!
!!! 2-d
!!!

subroutine EZspline_interp2_r8(spline_o, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  real(ezspline_r8) p1, p2  ! the location where the interpolation is sought
  real(ezspline_r8) f          ! the interpolation
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6)=(/1, 0, 0, 0, 0, 0 /)
  real(ezspline_r8) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3),  &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call r8pc2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call r8evbicub(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else

     call r8herm2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_r8

subroutine EZspline_interp2_GCvars_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ, spline_ogradBR, &
     spline_ogradBPHI, spline_ogradBZ, spline_ocurlbR, spline_ocurlbPHI, &
     spline_ocurlbZ, p1, p2, fBR, fBPHI, fBZ, &
     fER, fEPHI, fEZ, &
     fgradBR, fgradBPHI, fgradBZ, fcurlbR, fcurlbPHI, fcurlbZ, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline2_r8) spline_oER,spline_oEPHI,spline_oEZ
  type(EZspline2_r8) spline_ogradBR,spline_ogradBPHI,spline_ogradBZ
  type(EZspline2_r8) spline_ocurlbR,spline_ocurlbPHI,spline_ocurlbZ
  real(ezspline_r8), intent(in) :: p1, p2
  real(ezspline_r8), intent(out):: fBR, fBPHI, fBZ
  real(ezspline_r8), intent(out):: fER, fEPHI, fEZ
  real(ezspline_r8), intent(out):: fgradBR, fgradBPHI, fgradBZ
  real(ezspline_r8), intent(out):: fcurlbR, fcurlbPHI, fcurlbZ
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6)=(/1, 0, 0, 0, 0, 0 /)
  real(ezspline_r8) ansrBR(1), ansrBPHI(1), ansrBZ(1)
  real(ezspline_r8) ansrER(1), ansrEPHI(1), ansrEZ(1)
  real(ezspline_r8) ansrgradBR(1), ansrgradBPHI(1), ansrgradBZ(1)
  real(ezspline_r8) ansrcurlBR(1), ansrcurlBPHI(1), ansrcurlBZ(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier =94
     return
  endif


  call r8evbicub_GCvars(p1, p2,  &
       &   spline_oBR%x1(1), spline_oBR%n1, &
       &   spline_oBR%x2(1), spline_oBR%n2, &
       &   spline_oBR%ilin1, spline_oBR%ilin2, &
       &   spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
       &   spline_oBZ%fspl(1,1,1), &
       &   spline_oER%fspl(1,1,1), spline_oEPHI%fspl(1,1,1), &
       &   spline_oEZ%fspl(1,1,1), &
       &   spline_ogradBR%fspl(1,1,1), spline_ogradBPHI%fspl(1,1,1), &
       &   spline_ogradBZ%fspl(1,1,1), &
       &   spline_ocurlBR%fspl(1,1,1), spline_ocurlBPHI%fspl(1,1,1), &
       &   spline_ocurlBZ%fspl(1,1,1), &
       &   spline_oBR%n1, &
       &   ict, ansrBR, ansrBPHI, ansrBZ,  ansrER, ansrEPHI, ansrEZ, &
       &   ansrgradBR, ansrgradBPHI, ansrgradBZ, ansrcurlBR, &
       &   ansrcurlBPHI, ansrcurlBZ, ifail)

  fBR=ansrBR(1)
  fBPHI=ansrBPHI(1)
  fBZ=ansrBZ(1)
  fER=ansrER(1)
  fEPHI=ansrEPHI(1)
  fEZ=ansrEZ(1)
  fgradBR=ansrgradBR(1)
  fgradBPHI=ansrgradBPHI(1)
  fgradBZ=ansrgradBZ(1)
  fcurlBR=ansrcurlBR(1)
  fcurlBPHI=ansrcurlBPHI(1)
  fcurlBZ=ansrcurlBZ(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_GCvars_r8

subroutine EZspline_interp2_array_r8(spline_o, k1, k2, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r8) :: p1(k1), p2(k2) ! location arrays
  real(ezspline_r8) :: f(k1,k2)  ! interpolated function array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8gridintrp2d( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8gridpc2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then

     call r8gridbicub( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call r8gridherm2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_array_r8

subroutine EZspline_interp2_cloud_r8(spline_o, k, p1, p2, f, ier)
  ! list of coordinate doublets
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k) ! location arrays
  real(ezspline_r8), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif


  if (spline_o%isHybrid == 1) then

     call r8vecintrp2d(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8vecpc2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call r8vecbicub(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call r8vecherm2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_cloud_r8

subroutine EZspline_interp2_2DBdB_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_odBRdR, spline_odBPHIdR, spline_odBZdR, &
     spline_odBRdPHI, spline_odBPHIdPHI, spline_odBZdPHI, &
     spline_odBRdZ, spline_odBPHIdZ, spline_odBZdZ, &
     spline_oER, spline_oEPHI, spline_oEZ, k, p1, p2, fBR, fBPHI, fBZ, &
     fdBRdR, fdBPHIdR, fdBZdR,fdBRdPHI, fdBPHIdPHI, fdBZdPHI, &
     fdBRdZ, fdBPHIdZ, fdBZdZ,fER, fEPHI, fEZ, ier)
  use EZspline_obj
    implicit none
  type(EZspline2_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline2_r8) spline_odBRdR,spline_odBPHIdR,spline_odBZdR
  type(EZspline2_r8) spline_odBRdPHI,spline_odBPHIdPHI,spline_odBZdPHI
  type(EZspline2_r8) spline_odBRdZ,spline_odBPHIdZ,spline_odBZdZ
  type(EZspline2_r8) spline_oER,spline_oEPHI,spline_oEZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  real(ezspline_r8), intent(out):: fdBRdR(k), fdBPHIdR(k), fdBZdR(k)
  real(ezspline_r8), intent(out):: fdBRdPHI(k), fdBPHIdPHI(k), fdBZdPHI(k)
  real(ezspline_r8), intent(out):: fdBRdZ(k), fdBPHIdZ(k), fdBZdZ(k)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vecbicub_2DBdB(ict, k, p1, p2, k, fBR, fBPHI, fBZ, &
       fdBRdR, fdBPHIdR, fdBZdR, &
       fdBRdPHI, fdBPHIdPHI, fdBZdPHI, &
       fdBRdZ, fdBPHIdZ, fdBZdZ, &
       fER, fEPHI,fEZ, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
       & spline_oBZ%fspl(1,1,1), &
       & spline_odBRdR%fspl(1,1,1), spline_odBPHIdR%fspl(1,1,1), &
       & spline_odBZdR%fspl(1,1,1), &
       & spline_odBRdPHI%fspl(1,1,1), spline_odBPHIdPHI%fspl(1,1,1), &
       & spline_odBZdPHI%fspl(1,1,1), &
       & spline_odBRdZ%fspl(1,1,1), spline_odBPHIdZ%fspl(1,1,1), &
       & spline_odBZdZ%fspl(1,1,1), &
       spline_oER%fspl(1,1,1), &
       & spline_oEPHI%fspl(1,1,1), spline_oEZ%fspl(1,1,1), &   
       & spline_oBR%n1, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_2DBdB_cloud_r8

subroutine EZspline_interp2_2DBdB1_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ,spline_oER, spline_oEPHI, spline_oEZ, spline_oPSIp, &
     k, p1, p2, fBR, fBPHI, fBZ, &
     fER, fEPHI, fEZ, fPSIp, ier)
  use EZspline_obj
    implicit none
  type(EZspline2_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline2_r8) spline_oER,spline_oEPHI,spline_oEZ
  type(EZspline2_r8) spline_oPSIp
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out):: fBR(k,3), fBPHI(k,3), fBZ(k,3)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  real(ezspline_r8), intent(out):: fPSIp(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ictB(6) = (/1,1,1,0,0,0/)
  integer, parameter :: ictE(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vecbicub_2DBdB1(ictB, ictE, k, p1, p2, k, fBR, fBPHI, fBZ, &
       & fER, fEPHI, fEZ, fPSIp, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
       & spline_oBZ%fspl(1,1,1), &
       & spline_oER%fspl(1,1,1), &
       & spline_oEPHI%fspl(1,1,1), spline_oEZ%fspl(1,1,1), &
       & spline_oPSIp%fspl(1,1,1), &
       & spline_oBR%n1, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_2DBdB1_cloud_r8

subroutine EZspline_interp2_GCvars_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ, spline_ogradBR, &
     spline_ogradBPHI, spline_ogradBZ, spline_ocurlbR, spline_ocurlbPHI, &
     spline_ocurlbZ, k, p1, p2, fBR, fBPHI, fBZ, &
     fER, fEPHI, fEZ, &
     fgradBR, fgradBPHI, fgradBZ, fcurlbR, fcurlbPHI, fcurlbZ, ier)
  ! list of coordinate doublets
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline2_r8) spline_oER,spline_oEPHI,spline_oEZ
  type(EZspline2_r8) spline_ogradBR,spline_ogradBPHI,spline_ogradBZ
  type(EZspline2_r8) spline_ocurlbR,spline_ocurlbPHI,spline_ocurlbZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k) ! location arrays
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  real(ezspline_r8), intent(out):: fgradBR(k), fgradBPHI(k), fgradBZ(k)
  real(ezspline_r8), intent(out):: fcurlbR(k), fcurlbPHI(k), fcurlbZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vecbicub_GCvars(ict, k, p1, p2, k, fBR, fBPHI, fBZ, fER, fEPHI, &
       fEZ, fgradBR, fgradBPHI, fgradBZ, fcurlbR, fcurlbPHI, fcurlbZ, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
       & spline_oBZ%fspl(1,1,1), spline_oER%fspl(1,1,1), &
       & spline_oEPHI%fspl(1,1,1), spline_oEZ%fspl(1,1,1), &
       & spline_ogradBR%fspl(1,1,1), spline_ogradBPHI%fspl(1,1,1), &
       & spline_ogradBZ%fspl(1,1,1), spline_ocurlbR%fspl(1,1,1), &
       & spline_ocurlbPHI%fspl(1,1,1), spline_ocurlbZ%fspl(1,1,1), &
       & spline_oBR%n1, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_GCvars_cloud_r8

subroutine EZspline_interp2_FOvars_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ,k, p1, p2, fBR, &
     fBPHI, fBZ, fER, fEPHI, fEZ, ier)
  use EZspline_obj
  type(EZspline2_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline2_r8) spline_oER,spline_oEPHI,spline_oEZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vecbicub_FOvars(ict, k, p1, p2, k, fBR, fBPHI, fBZ, fER, fEPHI, &
       & fEZ, spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
       & spline_oBZ%fspl(1,1,1), spline_oER%fspl(1,1,1), &
       & spline_oEPHI%fspl(1,1,1), spline_oEZ%fspl(1,1,1), &
       & spline_oBR%n1, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_FOvars_cloud_r8

subroutine EZspline_interp2_FOmars_cloud_r8(spline_oA, spline_oBRRe, &
     spline_oBPHIRe, spline_oBZRe, spline_oBRIm, spline_oBPHIIm, &
     spline_oBZIm, k, p1, p2, fA, fBRRe, fBPHIRe, fBZRe, &
     fBRIm, fBPHIIm, fBZIm, ier)
  use EZspline_obj
  type(EZspline2_r8) spline_oA
  type(EZspline2_r8) spline_oBRRe,spline_oBPHIRe,spline_oBZRe
  type(EZspline2_r8) spline_oBRIm,spline_oBPHIIm,spline_oBZIm
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out):: fA(k,3)
  real(ezspline_r8), intent(out):: fBRRe(k), fBPHIRe(k), fBZRe(k)
  real(ezspline_r8), intent(out):: fBRIm(k), fBPHIIm(k), fBZIm(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ictA(6) = (/1,1,1,0,0,0/)
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oA) .or. spline_oA%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vecbicub_FOmars(ictA, ict, k, p1, p2, k, fA, &
       & fBRRe, fBPHIRe, fBZRe, fBRIm, fBPHIIm, fBZIm, &
       & spline_oA%n1, spline_oA%x1pkg(1,1), &
       & spline_oA%n2, spline_oA%x2pkg(1,1), &
       & spline_oA%fspl(1,1,1), &
       & spline_oBRRe%fspl(1,1,1), spline_oBPHIRe%fspl(1,1,1), &
       & spline_oBZRe%fspl(1,1,1), spline_oBRIm%fspl(1,1,1), &
       & spline_oBPHIIm%fspl(1,1,1), spline_oBZIm%fspl(1,1,1), &
       & spline_oA%n1, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_FOmars_cloud_r8

subroutine EZspline_interp2_collision_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ,&
     spline_one,spline_oTe,spline_oZeff,k, p1, p2, fBR, &
     fBPHI, fBZ, fER, fEPHI, fEZ, &
     fne, fTe, fZeff, ier)
  use EZspline_obj
  type(EZspline2_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline2_r8) spline_oER,spline_oEPHI,spline_oEZ
  type(EZspline2_r8) spline_one,spline_oTe,spline_oZeff
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  real(ezspline_r8), intent(out):: fne(k), fTe(k), fZeff(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vecbicub_collisions(ict, k, p1, p2, k, fBR, fBPHI, fBZ, fER, fEPHI, &
       & fEZ, fne, fTe, fZeff, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
       & spline_oBZ%fspl(1,1,1), spline_oER%fspl(1,1,1), &
       & spline_oEPHI%fspl(1,1,1), spline_oEZ%fspl(1,1,1), &
       & spline_one%fspl(1,1,1), spline_oTe%fspl(1,1,1), &
       & spline_oZeff%fspl(1,1,1), &
       & spline_oBR%n1, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_collision_cloud_r8

subroutine EZspline_interp2_bmag_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, k, p1, p2, fBR, fBPHI, fBZ, ier)
  use EZspline_obj
  type(EZspline2_r8) spline_oBR,spline_oBPHI,spline_oBZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vecbicub_bmag(ict, k, p1, p2, k, fBR, fBPHI, fBZ, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
       & spline_oBZ%fspl(1,1,1), spline_oBR%n1, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_bmag_cloud_r8

!!!
!!! 3-d
!!!

subroutine EZspline_interp3_r8(spline_o, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  real(ezspline_r8) p1, p2, p3 ! the location where the interpolation is sought
  real(ezspline_r8) f          ! the interpolation

  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  real(ezspline_r8) ansr(1)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8evintrp3d(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), size(spline_o%fspl,4), &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call r8pc3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call r8evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else

     call r8herm3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_r8

subroutine EZspline_interp3_array_r8(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer :: k1, k2, k3
  real(ezspline_r8) :: p1(k1), p2(k2), p3(k3)  ! location arrays
  real(ezspline_r8) :: f(k1,k2,k3)  ! interpolant array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8gridintrp3d( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8gridpc3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call r8gridtricub( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else

     call r8gridherm3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_array_r8

subroutine EZspline_interp3_cloud_r8(spline_o, k, p1, p2, p3, f, ier)
  ! list of coordinate triplets
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)  ! location arrays
  real(ezspline_r8), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8vecintrp3d(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8vecpc3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if (spline_o%isHermite == 0) then
     !
     call r8vectricub(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call r8vecherm3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_cloud_r8


subroutine EZspline_interp3_GCvars_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ, spline_ogradBR, &
     spline_ogradBPHI, spline_ogradBZ, spline_ocurlbR, spline_ocurlbPHI, &
     spline_ocurlbZ, k, p1, p2, p3, fBR, fBPHI, fBZ, &
     fER, fEPHI, fEZ, &
     fgradBR, fgradBPHI, fgradBZ, fcurlbR, fcurlbPHI, fcurlbZ, ier)
  ! list of coordinate doublets
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline3_r8) spline_oER,spline_oEPHI,spline_oEZ
  type(EZspline3_r8) spline_ogradBR,spline_ogradBPHI,spline_ogradBZ
  type(EZspline3_r8) spline_ocurlbR,spline_ocurlbPHI,spline_ocurlbZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k) ! location arrays
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  real(ezspline_r8), intent(out):: fgradBR(k), fgradBPHI(k), fgradBZ(k)
  real(ezspline_r8), intent(out):: fcurlbR(k), fcurlbPHI(k), fcurlbZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(10) = (/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vectricub_GCvars(ict, k, p1, p2, p3, k, fBR, fBPHI, fBZ, fER, fEPHI, &
       fEZ, fgradBR, fgradBPHI, fgradBZ, fcurlbR, fcurlbPHI, fcurlbZ, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%n3, spline_oBR%x3pkg(1,1), &
       & spline_oBR%fspl(1,1,1,1), spline_oBPHI%fspl(1,1,1,1), &
       & spline_oBZ%fspl(1,1,1,1), spline_oER%fspl(1,1,1,1), &
       & spline_oEPHI%fspl(1,1,1,1), spline_oEZ%fspl(1,1,1,1), &
       & spline_ogradBR%fspl(1,1,1,1), spline_ogradBPHI%fspl(1,1,1,1), &
       & spline_ogradBZ%fspl(1,1,1,1), spline_ocurlbR%fspl(1,1,1,1), &
       & spline_ocurlbPHI%fspl(1,1,1,1), spline_ocurlbZ%fspl(1,1,1,1), &
       & spline_oBR%n1, spline_oBR%n2, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_GCvars_cloud_r8

subroutine EZspline_interp3_3DBdB_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_odBRdR, spline_odBPHIdR, spline_odBZdR, &
     spline_odBRdPHI, spline_odBPHIdPHI, spline_odBZdPHI, &
     spline_odBRdZ, spline_odBPHIdZ, spline_odBZdZ, &
     spline_oER, spline_oEPHI, spline_oEZ, k, p1, p2, p3, fBR, fBPHI, fBZ, &
     fdBRdR, fdBPHIdR, fdBZdR,fdBRdPHI, fdBPHIdPHI, fdBZdPHI, &
     fdBRdZ, fdBPHIdZ, fdBZdZ,fER, fEPHI, fEZ, ier)
  use EZspline_obj
    implicit none
  type(EZspline3_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline3_r8) spline_odBRdR,spline_odBPHIdR,spline_odBZdR
  type(EZspline3_r8) spline_odBRdPHI,spline_odBPHIdPHI,spline_odBZdPHI
  type(EZspline3_r8) spline_odBRdZ,spline_odBPHIdZ,spline_odBZdZ
  type(EZspline3_r8) spline_oER,spline_oEPHI,spline_oEZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  real(ezspline_r8), intent(out):: fdBRdR(k), fdBPHIdR(k), fdBZdR(k)
  real(ezspline_r8), intent(out):: fdBRdPHI(k), fdBPHIdPHI(k), fdBZdPHI(k)
  real(ezspline_r8), intent(out):: fdBRdZ(k), fdBPHIdZ(k), fdBZdZ(k)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(10) = (/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vectricub_3DBdB(ict, k, p1, p2, p3, k, fBR, fBPHI, fBZ, &
       fdBRdR, fdBPHIdR, fdBZdR, &
       fdBRdPHI, fdBPHIdPHI, fdBZdPHI, &
       fdBRdZ, fdBPHIdZ, fdBZdZ, &
       fER, fEPHI,fEZ, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%n3, spline_oBR%x3pkg(1,1), &
       & spline_oBR%fspl(1,1,1,1), spline_oBPHI%fspl(1,1,1,1), &
       & spline_oBZ%fspl(1,1,1,1), &
       & spline_odBRdR%fspl(1,1,1,1), spline_odBPHIdR%fspl(1,1,1,1), &
       & spline_odBZdR%fspl(1,1,1,1), &
       & spline_odBRdPHI%fspl(1,1,1,1), spline_odBPHIdPHI%fspl(1,1,1,1), &
       & spline_odBZdPHI%fspl(1,1,1,1), &
       & spline_odBRdZ%fspl(1,1,1,1), spline_odBPHIdZ%fspl(1,1,1,1), &
       & spline_odBZdZ%fspl(1,1,1,1), &
       spline_oER%fspl(1,1,1,1), &
       & spline_oEPHI%fspl(1,1,1,1), spline_oEZ%fspl(1,1,1,1), &   
       & spline_oBR%n1, spline_oBR%n2, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_3DBdB_cloud_r8

subroutine EZspline_interp3_3DBdB1_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ, spline_oPSIp, &
     k, p1, p2, p3, fBR, fBPHI, fBZ, &
     fER, fEPHI, fEZ, fPSIp, ier)
  use EZspline_obj
    implicit none
  type(EZspline3_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline3_r8) spline_oER,spline_oEPHI,spline_oEZ
  type(EZspline3_r8) spline_oPSIp
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r8), intent(out):: fBR(k,4), fBPHI(k,4), fBZ(k,4)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  real(ezspline_r8), intent(out):: fPSIp(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ictB(10) = (/1,1,1,1,0,0,0,0,0,0/)
  integer, parameter :: ictE(10) = (/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vectricub_3DBdB1(ictB, ictE, k, p1, p2, p3, k, fBR, fBPHI, fBZ, &
       & fER, fEPHI, fEZ, fPSIp, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%n3, spline_oBR%x3pkg(1,1), &
       & spline_oBR%fspl(1,1,1,1), spline_oBPHI%fspl(1,1,1,1), &
       & spline_oBZ%fspl(1,1,1,1), &
       & spline_oER%fspl(1,1,1,1), &
       & spline_oEPHI%fspl(1,1,1,1), spline_oEZ%fspl(1,1,1,1), &
       & spline_oPSIp%fspl(1,1,1,1), &
       & spline_oBR%n1, spline_oBR%n2, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_3DBdB1_cloud_r8

subroutine EZspline_interp3_FOvars_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ, &
     k, p1, p2, p3, fBR, fBPHI, fBZ, fER, fEPHI, fEZ, ier)
  use EZspline_obj
  type(EZspline3_r8) spline_oBR,spline_oBPHI,spline_oBZ
  type(EZspline3_r8) spline_oER,spline_oEPHI,spline_oEZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  real(ezspline_r8), intent(out):: fER(k), fEPHI(k), fEZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vectricub_FOvars(ict, k, p1, p2, p3, k, &
       & fBR, fBPHI, fBZ, fER, fEPHI, fEZ,&
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%n3, spline_oBR%x3pkg(1,1), &
       & spline_oBR%fspl(1,1,1,1), spline_oBPHI%fspl(1,1,1,1), &
       & spline_oBZ%fspl(1,1,1,1), spline_oER%fspl(1,1,1,1), &
       & spline_oEPHI%fspl(1,1,1,1), spline_oEZ%fspl(1,1,1,1), &
       & spline_oBR%n1, spline_oBR%n2, iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_FOvars_cloud_r8

!/////
! R4 !
!/////

!!!
!!! 1-d
!!!

subroutine EZspline_interp1_r4(spline_o, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  real(ezspline_r4) p1  ! the location where the interpolation is sought
  real(ezspline_r4) f   ! the interpolation
  integer, intent(out) :: ier

  integer ifail
  integer, parameter :: ict(3)=(/1, 0, 0/)
  real(ezspline_r4) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isLinear == 1) then

     call pc1ev(p1, &
          &         spline_o%x1(1), spline_o%n1, &
          &         spline_o%ilin1,&
          &         spline_o%fspl(1,1),  &
          &         ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, ansr, ifail)

  else

     call herm1ev(p1, &
          &         spline_o%x1(1), spline_o%n1, &
          &         spline_o%ilin1,&
          &         spline_o%fspl(1,1),  &
          &         ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_r4


subroutine EZspline_interp1_array_r4(spline_o, k, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k) ! location arrays
  real(ezspline_r4), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: i, ifail
  integer, parameter :: ict(3)=(/1,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isLinear == 1) then

     call vecpc1(ict, k, p1, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn,ifail)

  else if (spline_o%isHermite == 0) then

     call vecspline(ict, k, p1, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else

     call vecherm1(ict, k, p1, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_array_r4

!!!
!!! 2-d
!!!

subroutine EZspline_interp2_r4(spline_o, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  real(ezspline_r4) p1, p2  ! the location where the interpolation is sought
  real(ezspline_r4) f          ! the interpolation
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6)=(/1, 0, 0, 0, 0, 0 /)
  real(ezspline_r4) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call pc2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then
     call evbicub(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else

     call herm2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_r4

subroutine EZspline_interp2_array_r4(spline_o, k1, k2, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r4) :: p1(k1), p2(k2) ! location arrays
  real(ezspline_r4) :: f(k1,k2)  ! interpolated function array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call gridintrp2d( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call gridpc2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then

     call gridbicub( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call gridherm2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97


end subroutine EZspline_interp2_array_r4

subroutine EZspline_interp2_cloud_r4(spline_o, k, p1, p2, f, ier)
  ! list of coordinate doublets
  use EZspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k), p2(k) ! location arrays
  real(ezspline_r4), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif


  if (spline_o%isHybrid == 1) then

     call vecintrp2d(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call vecpc2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call vecbicub(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call vecherm2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_cloud_r4


!!!
!!! 3-d
!!!

subroutine EZspline_interp3_r4(spline_o, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  real(ezspline_r4) p1, p2, p3 ! the location where the interpolation is sought
  real(ezspline_r4) f          ! the interpolation

  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  real(ezspline_r4) ansr(1)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call evintrp3d(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), size(spline_o%fspl,4), &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call pc3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else

     call herm3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_r4

subroutine EZspline_interp3_array_r4(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer :: k1, k2, k3
  real(ezspline_r4) :: p1(k1), p2(k2), p3(k3)  ! location arrays
  real(ezspline_r4) :: f(k1,k2,k3)  ! interpolant array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call gridintrp3d( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call gridpc3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call gridtricub( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else

     call gridherm3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_array_r4

subroutine EZspline_interp3_cloud_r4(spline_o, k, p1, p2, p3, f, ier)
  ! list of coordinate triplets
  use EZspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)  ! location arrays
  real(ezspline_r4), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call vecintrp3d(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call vecpc3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if (spline_o%isHermite == 0) then
     !
     call vectricub(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call vecherm3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_cloud_r4

subroutine EZspline_interp3_bmag_cloud_r8(spline_oBR, spline_oBPHI, &
     spline_oBZ, k, p1, p2, p3, fBR, fBPHI, fBZ, ier)
  use EZspline_obj
  type(EZspline3_r8) spline_oBR,spline_oBPHI,spline_oBZ
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r8), intent(out):: fBR(k), fBPHI(k), fBZ(k)
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_oBR) .or. spline_oBR%isReady /= 1) then
     ier = 94
     return
  endif


  call r8vectricub_bmag(ict, k, p1, p2, p3, k, fBR, fBPHI, fBZ, &
       & spline_oBR%n1, spline_oBR%x1pkg(1,1), &
       & spline_oBR%n2, spline_oBR%x2pkg(1,1), &
       & spline_oBR%n3, spline_oBR%x3pkg(1,1), &
       & spline_oBR%fspl(1,1,1,1), spline_oBPHI%fspl(1,1,1,1), &
       & spline_oBZ%fspl(1,1,1,1), spline_oBR%n1, &
       & spline_oBR%n2,iwarn, ifail)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_bmag_cloud_r8
