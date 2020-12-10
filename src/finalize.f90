module PB_finalize
  !! @note Module containing subroutines to terminate parallel
  !! communications and free memory.
  use PB_types
  use PB_fields
  use PB_hpc

  IMPLICIT NONE

  PUBLIC :: finalize_communications,&
       deallocate_variables

CONTAINS

  
  subroutine finalize_communications(params)
    !! @note Interface to function that finalizes MPI communications.
    !! See [[korc_hpc]].
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    !! Core KORC simulation parameters.  

    call finalize_mpi(params)
  end subroutine finalize_communications


  subroutine deallocate_variables(params,F,spp)
    !! @note Subroutine to free allocatable simulation variables.    
    TYPE(KORC_PARAMS), INTENT(INOUT) 			:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT) 			:: F
    !! An instance of KORC's derived type FIELDS containing all the
    !! information about the fields used in the simulation. See
    !! [[korc_types]] and [[korc_fields]].
    TYPE(SPECIES), INTENT(INOUT) :: spp
    !! An instance of KORC's derived type SPECIES containing all the
    !! information of different electron species. See [[korc_types]].

    
    DEALLOCATE(spp%vars%Y)
    DEALLOCATE(spp%vars%B)
    DEALLOCATE(spp%vars%PSI_P)
    DEALLOCATE(spp%vars%flagCon)
    DEALLOCATE(spp%vars%Y0)
    DEALLOCATE(spp%vars%k1)
    DEALLOCATE(spp%vars%k2)
    DEALLOCATE(spp%vars%k3)
    DEALLOCATE(spp%vars%k4)
    DEALLOCATE(spp%vars%k5)
    DEALLOCATE(spp%vars%k6)      

    call DEALLOCATE_FIELDS_ARRAYS(F)
  end subroutine deallocate_variables

end module PB_finalize
