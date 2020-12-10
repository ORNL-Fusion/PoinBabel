module PB_input
  !! @note Module with subroutines to read in all namelists in supplied
  !! input file and broadcast to all mpi processes.
  USE PB_types
  USE PB_hpc
  
  IMPLICIT NONE

  !! Default values for all inputs
  !! -----------------------------------------------
  !! input_parameters
  !! -----------------------------------------------
  INTEGER :: num_punctures = 1.E3
    ! Total aimed simulation time in seconds   
    ! Run 10 mu s If transients exist put 5 extra mu s.
  REAL(rp) :: dx = 1.E-2
    ! Time step as fraction of relativistic gyro-period
  CHARACTER(30) :: field_model = 'ANALYTICAL'
  CHARACTER(150) :: magnetic_field_filename = 'C1.h5'
  !  magnetic_field_filename = 'JFIT_D3D_164409_1405ms.h5'
  INTEGER :: time_slice = 000
  REAL(rp) :: rmax =  1.60
  REAL(rp) :: rmin =  0.15
  REAL(rp) :: zmax =  1.65
  REAL(rp) :: zmin = -1.65
  LOGICAL :: HDF5_error_handling = .TRUE.
  INTEGER :: pchunk = 1
  REAL(rp) :: phi_section

  !! -----------------------------------------------
  !! plasma_species
  !! As these inputs are vectors with dimension given by the number of species
  !! indicate default values for num_species=1 below, after the input_parameter
  !! namelist is read
  !! -----------------------------------------------
  INTEGER :: ppp = 1
  CHARACTER(30) :: spatial_distrib
    !! String describing the type of initial spatial distribution for
    !! each electron species.
    ! Options are: 'UNIFORM', 'DISK', 'TORUS', 'EXPONENTIAL-TORUS',
    ! 'GAUSSIAN-TORUS', 'ELLIPTIC-TORUS', 'EXPONENTIAL-ELLIPTIC-TORUS',
    ! 'GAUSSIAN-ELLIPTICAL-TORUS', '2D-GAUSSIAN-ELLIPTIC-TORUS-MH',
    ! 'AVALANCHE-4D','TRACER','SPONG-3D','HOLLMANN-3D','HOLLMANN-3D-PSI',
    ! 'MH_psi'
  REAL(rp), DIMENSION(3)  :: Xtrace = (/1.5_rp,0._rp,0._rp/)
    ! Initial position of tracer particle for debugging with
    ! spatial_distribution='TRACER'
  CHARACTER(30) :: position_filename = 'positions.dat'

  !! -----------------------------------------------
  !! analytical_fields_params
  !! -----------------------------------------------
  CHARACTER(30) :: current_direction = 'ANTI-PARALLEL' 
    ! 'PARALLEL' or 'ANTI-PARALLEL'
  REAL(rp) :: Bo = 2.2 
    ! In Teslas. ITER: 5.42 DIII-D: 2.19
  REAL(rp) :: minor_radius = 0.7 
    ! Minor radius in meters. ITER: 1.5 DIII-D: 0.5
  REAL(rp) :: major_radius = 1.7 
    ! Major radius in meters. ITER: 6.5 DIII-D: 1.6
  REAL(rp) :: qa = 5 
    ! Safety factor at the separatrix (r=a)
  REAL(rp) :: qo = 1.5 
    ! Safety factor at the separatrix (r=a)


  !! -----------------------------------------------
  !! externalPlasmaModel
  !! -----------------------------------------------
  LOGICAL :: Bfield = .FALSE.
  LOGICAL :: axisymmetric_fields = .FALSE.
  LOGICAL :: Bflux = .FALSE.
  REAL(rp) :: psip_conv = 1

CONTAINS

  subroutine read_namelist(params,infile,echo_in,outdir)

    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    CHARACTER(*), INTENT(IN) :: infile,outdir
    LOGICAL, INTENT(IN) :: echo_in

    INTEGER :: read_stat,nc
    INTEGER :: number_of_namelists=0,il,inst
    INTEGER, DIMENSION(20) :: namel_order=0
    CHARACTER(20) :: tempfile
    CHARACTER(128) :: ctmp
    CHARACTER(128) :: outfile
    LOGICAL :: reading
    INTEGER :: mpierr
    INTEGER :: tmp
    
    !! Namelist declarations
    NAMELIST /input_parameters/ field_model,magnetic_field_filename, &
         num_punctures,dx,HDF5_error_handling,time_slice,rmax, &
         rmin,zmax,zmin,pchunk,phi_section
    NAMELIST /plasma_species/ ppp,spatial_distrib,Xtrace,position_filename
    NAMELIST /analytical_fields_params/ Bo,minor_radius,major_radius,&
         qa,qo,current_direction
    NAMELIST /externalPlasmaModel/ Bfield, Bflux,axisymmetric_fields,psip_conv

!!-----------------------------------------------------------------------
!!     open input file.
!!     Remove comments from input file and put into temporary file.
!!-----------------------------------------------------------------------
    tempfile='tempinput.korc'
    if (params%mpi_params%rank.eq.0) then
       CALL rmcoment(infile,tempfile)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)   
    OPEN(UNIT=default_unit_open,FILE=tempfile,STATUS='OLD',POSITION='REWIND')
!!-----------------------------------------------------------------------
!!    check namelist file for namelist order and number.
!!-----------------------------------------------------------------------
    DO
       READ(UNIT=default_unit_open,FMT='(a)',IOSTAT=read_stat) ctmp 
       IF (read_stat/=0) EXIT
       nc=LEN_TRIM(ctmp)
       IF (nc<1) CYCLE
       ctmp=ADJUSTL(ctmp)
       reading=.false.
       IF (ctmp(1:1)=='&') THEN
          number_of_namelists=number_of_namelists+1
!!-----------------------------------------------------------------------
!!         trim all but the namelist name.
!!-----------------------------------------------------------------------
          DO il=2,nc+1
             IF (ctmp(il:il)/=' ') THEN
                IF (.NOT.reading) inst=il
                reading=.true.
                CYCLE
             ENDIF
             IF (ctmp(il:il)==' '.AND.reading) THEN
                ctmp=ctmp(inst:il-1)
                EXIT
             ENDIF
          ENDDO
          BACKSPACE(default_unit_open)
!!-----------------------------------------------------------------------
!!         select and read namelist.
!!-----------------------------------------------------------------------
          SELECT CASE(TRIM(ctmp))
          CASE('input_parameters')
             READ(UNIT=default_unit_open,NML=input_parameters,IOSTAT=read_stat)
          CASE('plasma_species')
             READ(UNIT=default_unit_open,NML=plasma_species,IOSTAT=read_stat)
          CASE('analytical_fields_params')
             READ(UNIT=default_unit_open,NML=analytical_fields_params,IOSTAT=read_stat)
          CASE('externalPlasmaModel')
             READ(UNIT=default_unit_open,NML=externalPlasmaModel,IOSTAT=read_stat)
          CASE DEFAULT
             write(output_unit_write,*) (TRIM(ctmp)//' is an unrecognized &
                  &namelist.')
             call PB_abort(13)
          END SELECT
          IF (read_stat/=0) then
             write(output_unit_write,*) ('Error reading namelist '//TRIM(ctmp)//'.')             
             call PB_abort(13)
          end if
       ENDIF
    ENDDO

!!-----------------------------------------------------------------------
!!     close input file.
!!       Delete it since it is the temporary file
!!-----------------------------------------------------------------------
    if (params%mpi_params%rank.ne.0) then
       CLOSE(default_unit_open)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    if (params%mpi_params%rank.eq.0) then
       CLOSE(default_unit_open,STATUS='DELETE')
    end if
!!-----------------------------------------------------------------------
!!     echo the input parameters to the output file.
!!-----------------------------------------------------------------------


    IF (echo_in) THEN
       if (params%mpi_params%rank .EQ. 0) then

          WRITE(output_unit_write,'(a,/)') 'VALUE OF ALL INPUTS:'
          WRITE(UNIT=output_unit_write,NML=input_parameters)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=plasma_species)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=analytical_fields_params)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=externalPlasmaModel)
          WRITE(output_unit_write,'(/)')
             
       end if
    end if

!!---------------------------------------------------------
!!     some tests
!!---------------------------------------------------------


    !write(6,*) TRIM(magnetic_field_filename),len(TRIM(magnetic_field_filename))

    tmp=len(TRIM(magnetic_field_filename))
    if (magnetic_field_filename(tmp-2:tmp).ne.'.h5'.and. &
         magnetic_field_filename(tmp-5:tmp-5).ne.'.') then
       if(params%mpi_params%rank.eq.0) then
          write(6,*) &
               'Check that enough characters are allocated for&
               & magnetic field filename!'
       end if
       call PB_abort(13)
    end if 
      
    end subroutine read_namelist

    SUBROUTINE rmcoment(fileold,filenew)
      
    CHARACTER(*), INTENT(IN) :: fileold,filenew
    CHARACTER(128) :: line
    INTEGER, PARAMETER :: nold=55,nnew=56
    INTEGER cmax, ios
    LOGICAL :: file_exist
!!-----------------------------------------------------------------------
!!     Open files, but make sure the old one exists first.
!!-----------------------------------------------------------------------
    INQUIRE(FILE=fileold,EXIST=file_exist)
    IF(.NOT. file_exist) THEN
       PRINT *,'The file "',fileold,'" could not be found.'
       STOP
    ENDIF

    OPEN(UNIT=default_unit_open,FILE=fileold,status="OLD",form='formatted')
    OPEN(UNIT=default_unit_write,FILE=filenew,status='REPLACE')

!!-----------------------------------------------------------------------
!!     Strip comments.     Note: line lengths limited to 127 characters
!!-----------------------------------------------------------------------
    DO
       READ(UNIT=default_unit_open,FMT='(a)',IOSTAT=ios) line
       IF (ios /= 0) EXIT
       cmax=1
       DO WHILE(line(cmax:cmax).NE.'!' .AND. cmax .LE. 127)
          cmax=cmax+1
       ENDDO
       IF(cmax .GT. 1) WRITE(default_unit_write,'(a)') line(1:cmax-1)
    ENDDO

!!-----------------------------------------------------------------------
!!     Close files and exit
!!-----------------------------------------------------------------------
    CLOSE(default_unit_open)
    CLOSE(default_unit_write)

  END SUBROUTINE rmcoment
  
end module PB_input
