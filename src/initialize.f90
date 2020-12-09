module initialize
  !! @note Module with subroutines to load simulation parameters 
  !! and to define the time step in the simulation.@endnote
  use types
  use constants
  use hpc
  use PB_HDF5
  use field
  use spatial_distribution
  use coords
  use input

  IMPLICIT NONE

  PRIVATE :: set_paths,&
       load_params
  PUBLIC :: initialize_parameters,&
       initialize_particles

CONTAINS

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! ** SUBROUTINES FOR INITIALIZING KORC PARAMETERS ** !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !


  subroutine load_params(params)
    TYPE (KORC_PARAMS), INTENT(INOUT) 	:: params
    INTEGER 				:: imax
    !! Auxiliary variable used to parse the output_list
    INTEGER 				:: imin
    !! Auxiliary variable used to parse the output_list
    INTEGER 				:: ii
    !! Iterator used to parse the output_list
    INTEGER 				:: jj
    !! Iterator used to parse the output_list
    INTEGER 				:: num_outputs
    !! Auxiliary variable used to parse the output_list
    INTEGER, DIMENSION(2) 		:: indices



    params%num_punctures = num_punctures
    params%dx = dx

    params%field_model = TRIM(field_model)
    params%magnetic_field_filename = TRIM(magnetic_field_filename)
    params%time_slice = time_slice
    params%rmax = rmax
    params%rmin = rmin
    params%zmax = zmax
    params%zmin = zmin

    if (HDF5_error_handling) then
       params%HDF5_error_handling = 1_idef
    else
       params%HDF5_error_handling = 0_idef
    end if

    params%pchunk=pchunk

    params%phi_section=phi_section


    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'(/,"* * * * * SIMULATION PARAMETERS * * * * *")')
       write(output_unit_write,'("Number of punctures: ",I4)') params%num_punctures
       write(output_unit_write,'("Step length: ",E17.10)') params%dx
       write(output_unit_write,*) 'Magnetic field model: ',TRIM(params%field_model)
       write(output_unit_write,*) 'Magnetic field file: ',TRIM(params%magnetic_field_filename)
    end if
    write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * *",/)')

  end subroutine load_params

  subroutine initialize_parameters(params)
    !! @note Interface for calling initialization subroutines @endnote
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    INTEGER 							:: mpierr
    !! MPI error status.

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    call read_namelist(params,params%path_to_inputs,.true.,params%path_to_outputs)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    call load_params(params)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  end subroutine initialize_parameters

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * * * SUBROUTINES FOR INITIALIZING PARTICLES * * * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  subroutine initialize_particles(params,F,spp)
    !! @note Subroutine that loads the information of the initial condition 
    !! of the different particle species. This subroutine calls
    !! the subroutine that generates the initial energy and pitch angle 
    !! distribution functions. @endnote
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN) 					:: F
    !! An instance of KORC's derived type FIELDS containing all the information 
    !! about the fields used in the simulation. See [[korc_types]]
    !!and [[korc_fields]].
    TYPE(SPECIES), INTENT(OUT) 	:: spp

    INTEGER 				                       	:: ii
    !! Iterator of spp structure.
    INTEGER 							:: mpierr
    !! MPI error status.
    !REAL(rp), DIMENSION(:), ALLOCATABLE :: Xtrace

    ! Allocate array containing variables of particles for each species

    spp%ppp=ppp
    spp%spatial_distrib = TRIM(spatial_distrib)

    spp%Xtrace = Xtrace

    ALLOCATE( spp%vars%punct(params%num_punctures,spp%ppp,2))
    ALLOCATE( spp%vars%Y(spp%ppp,3) )
    ALLOCATE( spp%vars%B(spp%ppp,3) )
    ALLOCATE( spp%vars%PSI_P(spp%ppp) )
    ALLOCATE( spp%vars%flagCon(spp%ppp) )
#ifdef FIO
    ALLOCATE( spp%vars%hint(spp%ppp))
#endif

    !     write(output_unit_write,'("0 size of PSI_P: ",I16)') size(spp%vars%PSI_P)
    spp%vars%punct = 0.0_rp
    spp%vars%Y = 0.0_rp    
    spp%vars%B = 0.0_rp
    spp%vars%PSI_P = 0.0_rp
    spp%vars%flagCon = 1_is

    ALLOCATE( spp%vars%Y0(spp%ppp,3) )
    ALLOCATE( spp%vars%k1(spp%ppp,3) )
    ALLOCATE( spp%vars%k2(spp%ppp,3) )
    ALLOCATE( spp%vars%k3(spp%ppp,3) )
    ALLOCATE( spp%vars%k4(spp%ppp,3) )
    ALLOCATE( spp%vars%k5(spp%ppp,3) )
    ALLOCATE( spp%vars%k6(spp%ppp,3) )

    spp%vars%Y0 = 0.0_rp
    spp%vars%k1 = 0.0_rp
    spp%vars%k2 = 0.0_rp
    spp%vars%k3 = 0.0_rp
    spp%vars%k4 = 0.0_rp
    spp%vars%k5 = 0.0_rp
    spp%vars%k6 = 0.0_rp

  end subroutine initialize_particles


  subroutine set_up_particles_ic(params,F,spp)
    !! @note Subroutine with calls to subroutines to load particles' 
    !! information if it is a restarting simulation, or to initialize the
    !! spatial and velocity distribution of each species if it is a new  
    !! simulation. @endnote
    TYPE(KORC_PARAMS), INTENT(INOUT) 				:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT) 					:: F
    !! An instance of KORC's derived type FIELDS containing all 
    !! the information about the fields used in the simulation. 
    !! See [[korc_types]] and [[korc_fields]].
    TYPE(SPECIES), INTENT(INOUT)       :: spp
    !! An instance of KORC's derived type SPECIES containing all 
    !! the information of different electron species. See [[korc_types]].
    !! An instance of the KORC derived type PROFILES.
    INTEGER                                                    :: ii
    !! Species iterator.

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("* * * * INITIALIZING SPATIAL DISTRIBUTION * * * *")')
    end if
    call intitial_spatial_distribution(params,spp,F)
    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
    end if


  end subroutine set_up_particles_ic

end module initialize
