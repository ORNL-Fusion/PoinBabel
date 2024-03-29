#ifdef FIO
!*******************************************************************************
!  @file korc_fio_interface.f90
!  @brief Interface for the fio interpolation library.
!*******************************************************************************

MODULE PB_fio
  USE, INTRINSIC :: iso_c_binding
  USE PB_types
  USE PB_input
  USE PB_HDF5
  USE mpi

  IMPLICIT NONE

  INTEGER (C_INT), PARAMETER :: FIO_SUCCESS        = 0
  INTEGER (C_INT), PARAMETER :: FIO_OUT_OF_BOUNDS  = 10002
  INTEGER (C_INT), PARAMETER :: FIO_NO_DATA        = 10006

  INTEGER (C_INT), PARAMETER :: FIO_M3DC1_SOURCE   = 3
  INTEGER (C_INT), PARAMETER :: FIO_NIMROD_SOURCE   = 5

  INTEGER (C_INT), PARAMETER :: FIO_TIMESLICE      = 1

  INTEGER (C_INT), PARAMETER :: FIO_SPECIES        = 3

  INTEGER (C_INT), PARAMETER :: FIO_ELECTRON       = 1
  INTEGER (C_INT), PARAMETER :: FIO_MAIN_ION       = -1

  INTEGER (C_INT), PARAMETER :: FIO_DENSITY        = 102
  INTEGER (C_INT), PARAMETER :: FIO_TEMPERATURE    = 103

  INTEGER (C_INT), PARAMETER :: FIO_ELECTRIC_FIELD = 1001
  INTEGER (C_INT), PARAMETER :: FIO_MAGNETIC_FIELD = 1003
  INTEGER (C_INT), PARAMETER :: FIO_VECTOR_POTENTIAL = 1002

  INTEGER (C_INT), PARAMETER :: FIO_TIME = 7001

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_add_field(icfield, ifield, op, fac)      &
          BIND(C, NAME='fio_add_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: icfield
       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
       INTEGER (C_INT), VALUE, INTENT(IN) :: op
       REAL (C_DOUBLE), VALUE, INTENT(IN) :: fac
     END FUNCTION fio_add_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_allocate_search_hint(isrc, hint)         &
          BIND(C, NAME='fio_allocate_search_hint')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
       TYPE (C_PTR), INTENT(OUT)          :: hint
     END FUNCTION fio_allocate_search_hint
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_field(ifield)                      &
          BIND(C, NAME='fio_close_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
     END FUNCTION fio_close_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_series(iseries)                    &
          BIND(C, NAME='fio_close_series')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
     END FUNCTION fio_close_series
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_source(isrc)                       &
          BIND(C, NAME='fio_close_source')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
     END FUNCTION fio_close_source
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_create_compound_field(ifield)            &
          BIND(C, NAME='fio_create_compound_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), INTENT(IN) :: ifield
     END FUNCTION fio_create_compound_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_deallocate_search_hint(isrc, hint)       &
          BIND(C, NAME='fio_deallocate_search_hint')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
       TYPE (C_PTR), INTENT(INOUT)        :: hint
     END FUNCTION fio_deallocate_search_hint
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_eval_field(ifield, x, v, hint)         &
          BIND(C, NAME='fio_eval_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
       REAL (C_DOUBLE), INTENT(IN)        :: x
       REAL (C_DOUBLE), INTENT(OUT)       :: v
       TYPE (C_PTR), VALUE, INTENT(IN)    :: hint
     END FUNCTION fio_eval_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_eval_field_deriv(ifield, x, v, hint)   &
          BIND(C, NAME='fio_eval_field_deriv')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
       REAL (C_DOUBLE), INTENT(IN)        :: x
       REAL (C_DOUBLE), INTENT(OUT)       :: v
       TYPE (C_PTR), VALUE, INTENT(IN)    :: hint
     END FUNCTION fio_eval_field_deriv
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_eval_series(iseries, x, v)             &
          BIND(C, NAME='fio_eval_series')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
       REAL (C_DOUBLE), INTENT(IN)        :: x
       REAL (C_DOUBLE), INTENT(OUT)       :: v
     END FUNCTION fio_eval_series
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_options(isrc)                        &
          BIND(C, NAME='fio_get_options')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
     END FUNCTION fio_get_options
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_available_fields(isrc, n, f)         &
          BIND(C, NAME='fio_get_available_fields')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN)         :: isrc
       INTEGER (C_INT), INTENT(OUT)               :: n
       INTEGER (C_INT), DIMENSION(:), INTENT(OUT) :: f
     END FUNCTION fio_get_available_fields
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_open_source(itype, filename, handle)     &
          BIND(C, NAME='fio_open_source')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN)        :: itype
       CHARACTER (kind=C_CHAR,len=1), INTENT(IN) :: filename
       INTEGER (C_INT), INTENT(OUT)              :: handle
     END FUNCTION fio_open_source
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_field(isrc, type, handle)            &
          BIND(C, NAME='fio_get_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(in) :: isrc
       INTEGER (C_INT), VALUE, INTENT(IN) :: type
       INTEGER (C_INT), INTENT(INOUT)     :: handle
     END FUNCTION fio_get_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_set_int_option(iopt, v)                  &
          BIND(C, NAME='fio_set_int_option')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(in) :: iopt
       INTEGER (C_INT), VALUE, INTENT(in) :: v
     END FUNCTION fio_set_int_option
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_real_field_parameter(ifield, t, p) &
          BIND(C, NAME='fio_get_real_field_parameter')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       integer(c_int), intent(in), value :: ifield, t
       real(c_double) :: p
     END FUNCTION fio_get_real_field_parameter
  END INTERFACE

CONTAINS

  SUBROUTINE initialize_m3d_c1(params, F, spp,init)

    TYPE(KORC_PARAMS), INTENT(INOUT)           :: params
    TYPE(FIELDS), INTENT(INOUT)                :: F
    TYPE(SPECIES), INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN)  :: init
    
    INTEGER                                    :: ii
    INTEGER                                    :: pp
    INTEGER                                    :: status
    INTEGER                                    :: isrc
    real(c_double)  ::  time0,time1
    INTEGER (C_INT)                         :: FIO_tmp
    TYPE(C_PTR) :: hint_tmp
    real(rp), DIMENSION(3) :: x
    REAL(rp), DIMENSION(3)         :: Btmp   

    status = fio_open_source(FIO_M3DC1_SOURCE,           &
         TRIM(params%magnetic_field_filename)            &
         // C_NULL_CHAR, F%isrc)       


    isrc=F%isrc
    
    status = fio_get_options(isrc)
       
    status = fio_set_int_option(FIO_TIMESLICE, params%time_slice)
    
    status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%FIO_B)    

    hint_tmp=c_null_ptr
    x(1)=spp%Xtrace(1)
    x(2)=spp%Xtrace(2)
    x(3)=spp%Xtrace(3)

    status = fio_eval_field(F%FIO_B, x(1),                      &
         Btmp(1),hint_tmp)

    F%Bo = -Btmp(2)
    
    do pp = 1, spp%ppp
       status = fio_allocate_search_hint(isrc, spp%vars%hint(pp))
    end do

    spp%vars%cart = .false.

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,*) 'Calculate B',F%FIO_B
    end if
       
  END SUBROUTINE initialize_m3d_c1
  
  SUBROUTINE initialize_nimrod(params, F, spp,init)

    TYPE(KORC_PARAMS), INTENT(INOUT)           :: params
    TYPE(FIELDS), INTENT(INOUT)                :: F
    TYPE(SPECIES), INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN)  :: init
    
    INTEGER                                    :: ii
    INTEGER                                    :: pp
    INTEGER                                    :: status
    INTEGER                                    :: isrc
    real(c_double)  ::  time0,time1
    INTEGER (C_INT)                         :: FIO_tmp
    TYPE(C_PTR) :: hint_tmp
    real(rp), DIMENSION(3) :: x
    REAL(rp), DIMENSION(3)         :: Btmp
    CHARACTER(150) :: filename

    
    status = fio_open_source(FIO_NIMROD_SOURCE,           &
         TRIM(params%magnetic_field_filename)            &
         // C_NULL_CHAR, F%isrc)       


    isrc=F%isrc
    
    status = fio_get_options(isrc)       
    
    status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%FIO_B)    

    hint_tmp=c_null_ptr
    x(1)=spp%Xtrace(1)
    x(2)=spp%Xtrace(2)
    x(3)=spp%Xtrace(3)

    status = fio_eval_field(F%FIO_B, x(1),Btmp(1),hint_tmp)

    F%Bo = -Btmp(2)
    
    do pp = 1,spp%ppp
       status = fio_allocate_search_hint(isrc, spp%vars%hint(pp))
    end do

    spp%vars%cart = .false.

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,*) 'Calculate B',F%FIO_B
    end if
       
  END SUBROUTINE initialize_nimrod
  

  SUBROUTINE finalize_fio(params, F)
    TYPE(KORC_PARAMS), INTENT(IN)           :: params
    TYPE(FIELDS), INTENT(IN)                :: F
    INTEGER                                    :: status
    INTEGER                                    :: ii

    status = fio_close_field(F%FIO_B)       

    status=fio_close_source(F%isrc)
    
  end SUBROUTINE FINALIZE_FIO


END MODULE PB_fio
#endif
