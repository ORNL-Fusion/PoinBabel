program main
  !! @note  Main function of KORC. @endnote
  !! The main program contains the calls to the main functions and subroutines. 
  !! Also, it contains the variables that control
  !! the behavior of the core of KORC and all other external/optional modules.

  use PB_types
  use PB_fields
  use PB_hpc
  use PB_HDF5
  use PB_ppusher
  use PB_interp
  use PB_initialize
  use PB_finalize
  use PB_input
#ifdef FIO
  use PB_fio
#endif
  
  implicit none

  TYPE(KORC_PARAMS) :: params
  !! Contains the parameters that control the core of KORC: 
  !! time steping, output list, etc.
  TYPE(SPECIES) :: spp
  !! Contains the initial parameters of each species, which 
  !! can be different electrons with different
  !! distribution functions.
  TYPE(FIELDS) :: F
  !! F: Contains the parameters of the analytical magnetic 
  !! and electric fields, or in the case of using 
  !! external fields it contains the data used in the interpolations. 
  !!See [[korc_fields(module)]] for details.
  INTEGER(ip) :: it 
  !! Time iteration
  INTEGER 				:: mpierr
    
  call initialize_communications(params)
  !!<h2>Order of KORC operations</h2>
  !!
  !!<h3>Communication and Timing</h3>
  !! <h4>1\. Parallel Communications</h4>
  !!
  !! Subroutine [[initialize_communications]] in [[korc_hpc]] that 
  !! initializes MPI and OpenMP communications.
  
  call timing(params)
  !! <h4>2\. Timers</h4>
  !!
  !! Subroutine [[timing_KORC]] in [[korc_hpc]] that times the 
  !! execution of any parallel sections of KORC.
  
  ! * * * INITIALIZATION STAGE * * *!

  call initialize_HDF5()
  !!<h3>Initialization</h3>
  !!
  !! <h4>1\. HDF5</h4>
  !!
  !! Subroutine [[initialize_HDF5]] in [[korc_HDF5]] that initializes
  !! HDF5 library. 
  
  call initialize_parameters(params)
  !! <h4>2\. Initialize korc parameters</h4>
  !!
  !! Subroutine [[initialize_korc_parameters]] in [[korc_initialize]] that 
  !! initializes paths and KORC parameters through [[load_korc_params]]
  !! on MPI processes.

  call initialize_fields(params,F)
  !! <h4>3\. Initialize fields</h4>
  !!
  !! Subroutine [[initialize_fields]] in [[korc_fields]] that initializes 
  !! parameters of the EM fields, either analytically or from an external HDF5
  !! file. Reads in &amp;analytical_fields_params and 
  !! &amp;externalPlasmaModel namelists from input file.
  
  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if

  call initialize_particles(params,F,spp) ! Initialize particles
  !! <h4>5\. Initialize Particle Velocity Phase Space</h4>
  !! 
  !! Subroutine [[initialize_particles]] in [[korc_initialize]] that 
  !! initializes particle parameters from &amplasma_species namelist, 
  !! allocates arrays for individual particles, including location, velocity, 
  !! local EM fields and plasma profiles, etc., and 
  !! calls [[initial_energy_pitch_dist]] to assign particles' energy and pitch
  !! angle according to the chosen distribution.

!  write(output_unit_write,'("init eta: ",E17.10)') spp(1)%vars%eta

#ifdef FIO
  if (TRIM(params%field_model) .eq. 'M3D_C1') then

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * INITIALIZING M3D-C1 INTERFACE * * * *"
     endif
     
     call initialize_m3d_c1(params, F, spp,.true.)

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * * * * * * * * * * * * * * * * * * * *"
     endif

  elseif (TRIM(params%field_model) .eq. 'NIMROD') THEN

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * INITIALIZING NIMROD INTERFACE * * * *"
     endif
     
     call initialize_nimrod(params, F, spp,.true.)

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * * * * * * * * * * * * * * * * * * * *"
     endif
     
  endif
#endif  

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if

#ifdef PSPLINE
  call initialize_fields_interpolant(params,F)
#endif
  !! <h4>15\. Initialize Fields Interpolant</h4>
  !!
  !! Subroutine [[initialize_fields_interpolant]] in [[korc_interp]] calls
  !! EZspline
  !! subroutines EZspline_init for memory allocation and boundary condition
  !! setup
  !! and EZspline_setup to compute the necessary cubic coefficients needed
  !! for subsequent
  !! field interpolations. The magnetic field can be defined in terms of an
  !! axisymmetric
  !! scalar flux function, axisymmetric field, or 3D field, while the
  !! electric field
  !! can be defined as an axisymmetric or 3D field.
  
  if (params%mpi_params%rank .EQ. 0) then
     write(output_unit_write,'("* * * * INITIALIZING INITIAL CONDITIONS * * * *",/)')
  end if
  call set_up_particles_ic(params,F,spp)
  
  if (params%mpi_params%rank .EQ. 0) then
     write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
  end if
  
  ! * * * INITIALIZATION STAGE * * *
  
  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  ! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !
  
  call save_simulation_parameters(params,spp,F)       

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  ! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !


  call timing(params)

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if

  if (params%field_model.eq.'ANALYTICAL') then     
     call adv_eqn_top(params,F,spp)
  else if (params%field_model.eq.'PSPLINE') then
#ifdef PSPLINE
     if (F%Bflux) then
        call adv_interp_psi_top(params,F,spp)
     else if (F%Bfield.and.F%axisymmetric_fields) then
        call adv_interp_2DB_top(params,F,spp)
     else if (F%Bfield.and.(.not.F%axisymmetric_fields)) then
        call adv_interp_3DB_top(params,F,spp)
     endif
#else
     write(6,*) 'Turn on PSPLINE!'
#endif     
  else if ((params%field_model.eq.'NIMROD').or. &
       (params%field_model.eq.'M3D_C1')) then
#ifdef FIO
     call adv_fio_top(params,F,spp)
#else
     write(6,*) 'Turn on FIO!'
#endif
  endif

  call save_simulation_outputs(params,spp,F)
  
  call timing(params)

  ! * * * FINALIZING SIMULATION * * * 
  call finalize_HDF5()

#ifdef PSPLINE
  call finalize_interpolants(params)
#endif
  
#ifdef FIO
  if (TRIM(params%field_model) .eq. 'M3D_C1'.or. &
      TRIM(params%field_model) .eq. 'NIMROD') then
     call finalize_FIO(params,F)
  end if
#endif
  
  ! DEALLOCATION OF VARIABLES
  call deallocate_variables(params,F,spp)

  
  call finalize_communications(params)
  ! * * * FINALIZING SIMULATION * * *

  if (params%mpi_params%rank .EQ. 0) then
     write(output_unit_write,'("PoinBabel ran successfully!")')
     close(output_unit_write)
  end if
  
end program main

