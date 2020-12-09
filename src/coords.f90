module coords
  !! @note Module containing subroutines to calculate the position of
  !! the simulated particles in toroidal and cylindrical coordinates. @endnote
  use types
  use constants

  IMPLICIT NONE

  PUBLIC :: cyl_check_if_confined_p
  
CONTAINS


  subroutine cyl_check_if_confined_p(pchunk,a,R0,Xcyl_R,Xcyl_Z,flag)
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp),DIMENSION(pchunk),  INTENT(IN)      :: Xcyl_R
    REAL(rp),DIMENSION(pchunk),  INTENT(IN)      :: Xcyl_Z
    INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag
    REAL(rp),  INTENT(IN)                            :: a,R0
    !! Distance to plasma edge as measured from the magnetic axis.
    INTEGER                                  :: cc
    
    !$OMP SIMD
!    !$OMP& aligned(Xcyl_R,Xcyl_Z,flag)
    do cc=1_idef,pchunk
       if (sqrt((Xcyl_R(cc)-R0)**2+Xcyl_Z(cc)**2) .gt. a) flag(cc)=0_is
    end do
    !$OMP END SIMD

  end subroutine cyl_check_if_confined_p
  
end module coords
