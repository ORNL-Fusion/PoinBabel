module types
  !! @note Module containing the definition of KORC derived types and
  !! KORC variables, the building blocks of the code. @endnote
#ifdef FIO
  USE, INTRINSIC :: iso_c_binding
#endif
  implicit none

  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Real and integer precisions * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  INTEGER, PUBLIC, PARAMETER 	:: is = KIND(INT(1,1)) 
  !! Definition of 1 Byte (8 bits) Fortran KORC integer type.
  INTEGER, PUBLIC, PARAMETER 	:: ip = KIND(INT(1,8)) 
  !! Definition of 8 Bytes (64 bits) Fortran KORC integer type.
  INTEGER, PUBLIC, PARAMETER 	:: idef = KIND(1) 
  !! Definition of the default KORC integer type on the system where
  !! KORC is compiled.
  INTEGER, PUBLIC, PARAMETER 	:: rdef = KIND(1.0) 
  !! Definition of the default KORC real type on the system where
  !! KORC is compiled.
#ifdef DOUBLE_PRECISION
  INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(0.d0) 
  !! Definition of the KORC double precision real type.
#elif SINGLE_PRECISION
  INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(1.0) 
  !! Definition of the KORC single precision real type.
#endif
  REAL(rp), PUBLIC, PARAMETER :: korc_zero = 1.0E-15 
  !! Definition of the zero in KORC.

  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Real and integer precisions * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  INTEGER, PUBLIC, PARAMETER 	:: MAX_STRING_LENGTH = 1000 
  !! Default length of a KORC_STRING variable.
  INTEGER, PUBLIC, PARAMETER 	:: default_unit_open = 101 
  !! Default file unit for opening and reading from an external text file.
  INTEGER, PUBLIC, PARAMETER 	:: default_unit_write = 201 
  !! Default file unit for opening and writing to external an external text file.
  INTEGER, PUBLIC, PARAMETER 	:: output_unit_write = 202 
  !! Default file unit for opening and writing to external an external text file.

  !> @note KORC string type of length MAX_STRING_LENGTH. @endnote
  TYPE, PUBLIC :: KORC_STRING 
     CHARACTER(MAX_STRING_LENGTH) :: str
  END TYPE KORC_STRING

  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Basic korc array structures * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  TYPE, PUBLIC :: V_FIELD_3D
     !! @note KORC 3-D vector field type @endnote
     !! This KORC type represents a 3-D vector field varible in
     !! cylindrical coordinates. For example, this could be the 3-D magnetic
     !! field, which can be written as $$\mathbf{B}(R,\phi,Z) = B_R(R,\phi,Z)
     !! \hat{R} + B_\phi(R,\phi,Z) \hat{\phi} + B_Z(R,\phi,Z) \hat{Z}.$$
     !! All the members (components) of the V_FIELD_3D type follow the
     !! following index convention:
     !! (\(R\) index,\(\phi\) index, \(Z\) index)

     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: R 
     !! \(R\) component of the vector field variable.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: PHI 
     !! \(\phi\) component of the vector field variable.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Z 
     !! \(Z\) component of the vector field variable.
  END TYPE V_FIELD_3D

  TYPE, PUBLIC :: V_FIELD_2D
     !! @note KORC 2-D vector field type @endnote
     !! This KORC type represents a 2-D vector field varible in cylindrical
     !! coordinates. For example, this could be the magnetic
     !! field in an axisymmetric plasma, which can be written as
     !! $$\mathbf{B}(R,Z) = B_R(R,Z) \hat{R} + B_\phi(R,Z) \hat{\phi} + B_Z(R,Z)
     !! \hat{Z}.$$
     !! All the members (components) of the V_FIELD_2D type follow the
     !! following index convention:
     !! (\(R\) index,\(Z\) index).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: R 
     !! \(R \) component of the vector field variable.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: PHI 
     !! \(\phi \) component of the vector field variable.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z 
     !! \(Z \) component of the vector field variable.
  END TYPE V_FIELD_2D

  TYPE, PUBLIC :: V_FIELD_1D
     !! @note KORC 1-D vector field type @endnote
     !! This KORC type represents a 1-D vector field varible in cylindrical
     !! coordinates. For example, this could be the magnetic
     !! field in an axisymmetric plasma, which can be written as
     !! $$\mathbf{B}(r) = B_R(r) \hat{R} + B_\phi(r) \hat{\phi} + B_Z(r)
     !! \hat{Z}.$$
     !! All the members (components) of the V_FIELD_1D type follow the
     !! following index convention:
     !! (\(r\) index).
     REAL(rp), DIMENSION(:), ALLOCATABLE :: R 
     !! \(R \) component of the vector field variable.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: PHI 
     !! \(\phi \) component of the vector field variable.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: Z 
     !! \(Z \) component of the vector field variable.
  END TYPE V_FIELD_1D
  
  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Basic korc array structures * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  TYPE, PRIVATE :: KORC_MPI
     !! @note KORC derived type to keep relevant MPI parameters. @endnote
     INTEGER :: nmpi 
     !! Number of MPI processes.
     INTEGER :: rank 
     !! Rank in WORLD COMMON communicator.
     INTEGER :: rank_topo 
     !! Rank in mpi_topo communicator
     INTEGER :: mpi_topo 
     !! MPI communicator for a certain topology.
  END TYPE KORC_MPI


  TYPE, PUBLIC :: KORC_PARAMS
     !! @note Core KORC parameters. @endnote
     !!  This KORC derived type contains the variables that control KORC's
     !! core behavior.

     CHARACTER(MAX_STRING_LENGTH) :: path_to_inputs 
     !! Absolute path to KORC's input file.
     CHARACTER(MAX_STRING_LENGTH) :: path_to_outputs 
     !! Absolute path to the outputs' folder.
     INTEGER 			:: num_omp_threads 
     !! Number of open MP threads per MPI process used in the simulation.
     INTEGER 			:: num_punctures
     !! Total simulation time in seconds.
     REAL(rp) 			:: dx 
     !! Time step in the simulation as a fraction of the relativistic electron
     !! gyro-period \(\tau_e = 2\pi\gamma m_e/eB_0\).
     REAL(rp) 			:: time = 0.0_rp 
     !! Current physical time in the simulation.
     INTEGER(ip) 			:: ito = 0_ip 
     !! Initial time iteration in the simulation, this is different from zero
     !! in case is a restarting simulation.
     INTEGER(ip) 			:: it = 0_ip 
     !! Current time iteration in the simulation, this is different from zero
     !! in case is a restarting simulation.   
     INTEGER(ip) 			:: t_steps
     INTEGER(ip) 			:: t_it_SC=1_ip 
     CHARACTER(MAX_STRING_LENGTH) :: field_model
     CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_filename 
     INTEGER 			:: HDF5_error_handling 
     !! Flag to indicate whether we allow HDF5 to show warnings 
     !! during runtime (HDF5_error_handling=1) or not (HDF5_error_handling=0)
     TYPE(KORC_MPI) 		:: mpi_params 
     !! An instance of the KORC_MPI derived type.
     INTEGER                      :: time_slice !< M3D-C1 time slice to use.
     REAL(rp)                     :: rmax !< Maximum r for M3D-C1 fields.
     REAL(rp)                     :: rmin !< Minimum r for M3D-C1 fields.
     REAL(rp)                     :: zmax !< Maximum z for M3D-C1 fields.
     REAL(rp)                     :: zmin !< Minimum z for M3D-C1 fields.
     INTEGER  :: pchunk
     !! number of particles per vectorized chunk
  END TYPE KORC_PARAMS


  TYPE, PUBLIC :: PARTICLES
     !! @note Derived type containing all the electrons' variables
     !!in the simulation. @endnote
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE    :: punct
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Y 
     !! Coordinates of the electrons' position in cylindrical or toroidal
     !! coordinates.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Y0 
     !! Placeholder coordinates of the electrons' position in cylindrical
     !! coordinates for GC orbit model.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: B 
     !! Cartesian components of the magnetic field felt by each electron.
     !! dim(B) = dim(X).
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PSI_P 
     !! RHS of equations of motion for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k1
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k2
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k3
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k4
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k5
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k6
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     INTEGER(is), DIMENSION(:), ALLOCATABLE 	:: flagCon
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: AUX 
     !! An auxiliary scalar variable for each electron.
#ifdef FIO
     TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: hint
     !! Hint for M3D_C1 interpolation.
#endif
     LOGICAL                                :: cart
  END TYPE PARTICLES


  TYPE, PUBLIC :: SPECIES
     !! @note Derived type containing the initial parameters of each electron
     !! ensemble followed in a KORC simulation. @endnote
     INTEGER :: ppp
     TYPE(PARTICLES) 			:: vars 
     !! An instance of the PARTICLES derived type.
     CHARACTER(MAX_STRING_LENGTH) 	:: spatial_distrib
     !! String describing the type of initial spatial distribution for
     !! each electron species.
     REAL(rp),DIMENSION(3) :: Xtrace
     !! Initial position in Cartesian coordinates for tracer particle
  END TYPE SPECIES


  TYPE, PRIVATE :: A_FIELD
     !! @note Derived type having all the parameters of the analytical
     !! magnetic field included in KORC. @endnote
     !! The analytical magnetic field is given by:
     !! $$\mathbf{B}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}}
     !! \left[ B_0 \hat{e}_\zeta  + B_\vartheta(r) \hat{e}_\vartheta \right],$$
     !! where \(\eta = r/R_0\) is the aspect ratio, the constant \(B_0\)
     !! denotes the magnitude of the toroidal magnetic field, 
     !! and \(B_\vartheta(r) = \eta B_0/q(r)\) is the poloidal magnetic
     !! field with safety factor 
     !! \(q(r) = q_0\left( 1 + \frac{r^2}{\lambda^2} \right)\). The
     !! constant \(q_0\) is the safety factor at 
     !! the magnetic axis and the constant \(\lambda\) is obtained from
     !! the values of \(q_0\) and \(q(r)\) 
     !! at the plasma edge \(r=r_{edge}\).

     REAL(rp) 			:: Bo 
     !! Magnitude of the toroidal magnetic field \(B_0\).
     REAL(rp) 			:: a 
     !! Plasma edge \(r_{edge}\) as measured from the magnetic axis.
     REAL(rp) 			:: Ro 
     !! Radial position of the magnetic axis \(R_0\)
     REAL(rp) 			:: qa 
     !! Safety factor at the plasma edge.
     REAL(rp) 			:: qo 
     !! Safety factor at the magnetic axis \(q_0\).
     REAL(rp) 			:: lambda 
     !! \(\lambda\) parameter of \(q(r)\).
     REAL(rp) 			:: Bpo 
     !! @deprecated Parameter not used anymore. @todo Delete parameter.
     REAL(rp) 			:: Bp_sign 
     !! Sign of \(B_\vartheta(r)\). This depends on current_direction,
     !! Bp_sign=1 for 
     !! current_direction='PARALLEL', and Bp_sign=-1 for
     !! current_direction='ANTI-PARALLEL'.
     CHARACTER(MAX_STRING_LENGTH) :: current_direction 
     !! Direction of plasma current: PARALLEL or ANTI-PARALLEL to the
     !! toroidal magnetic field.
  END TYPE A_FIELD

  TYPE, PRIVATE :: MESH
     !! Derived type with the cylindrical coordinates of the grid nodes 
     !! at which the pre-computed plasma profiles and fields are known.

     REAL(rp), DIMENSION(:), ALLOCATABLE :: R 
     !! Radial grid.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: PHI 
     !! Azimuthal grid.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: Z 
     !! Z grid.
  END TYPE MESH


  TYPE, PUBLIC :: FIELDS
     !! @note Derived type with all the variables and data of analytical 
     !! and pre-computed electric and magnetic fields. @endnote

     TYPE(A_FIELD) 	 			:: AB 
     !! An instance of the KORC derived data type A_FIELD.
     TYPE(V_FIELD_3D) 				:: E_3D 
     !! KORC 3-D vector field of the pre-computed electric field.
     TYPE(V_FIELD_3D) 				:: B_3D
     TYPE(V_FIELD_3D) 				:: dBdR_3D
     TYPE(V_FIELD_3D) 				:: dBdPHI_3D
     TYPE(V_FIELD_3D) 				:: dBdZ_3D 
     !! KORC 3-D vector field of the pre-computed magnetic field.
     TYPE(V_FIELD_2D) 				:: E_2D 
     !! KORC 2-D vector field of the pre-computed electric field.
     TYPE(V_FIELD_2D) 				:: B_2D
     TYPE(V_FIELD_2D) 				:: B1Re_2D
     TYPE(V_FIELD_2D) 				:: B1Im_2D
     TYPE(V_FIELD_2D) 				:: dBdR_2D
     TYPE(V_FIELD_2D) 				:: dBdPHI_2D
     TYPE(V_FIELD_2D) 				:: dBdZ_2D 
     !! KORC 3-D vector field of the pre-computed magnetic field.
     TYPE(V_FIELD_2D) 				:: gradB_2D
     TYPE(V_FIELD_3D) 				:: gradB_3D 
     !! KORC 3-D vector field of the gradient of the magnitude of the
     !! pre-computed magnetic field.
     TYPE(V_FIELD_2D) 				:: curlb_2D
     TYPE(V_FIELD_3D) 				:: curlb_3D 
     !! KORC 3-D vector field of the curl of the unit vector in the
     !! direction of the pre-computed magnetic field.
     TYPE(V_FIELD_1D) 				:: E_SC_1D
     TYPE(V_FIELD_1D) 				:: J0_SC_1D
     TYPE(V_FIELD_1D) 				:: J1_SC_1D
     TYPE(V_FIELD_1D) 				:: J2_SC_1D
     TYPE(V_FIELD_1D) 				:: J3_SC_1D
     TYPE(V_FIELD_1D) 				:: A1_SC_1D
     TYPE(V_FIELD_1D) 				:: A2_SC_1D
     TYPE(V_FIELD_1D) 				:: A3_SC_1D
     
     REAL(rp), DIMENSION(:), ALLOCATABLE :: r_1D
     REAL(rp), DIMENSION(:), ALLOCATABLE :: PSIP_1D
     REAL(rp), DIMENSION(:), ALLOCATABLE :: dMagPsiSqdPsiP
     REAL(rp), DIMENSION(:), ALLOCATABLE :: ddMagPsiSqdPsiPSq
     TYPE(MESH) 		 		:: X 
     !! An instance of the KORC derived type MESH.
     CHARACTER(MAX_STRING_LENGTH) :: E_model
     !! Name for dynamical, analytic, electric field model to be added to
     REAL(rp)  :: E_dyn
     REAL(rp)  :: E_pulse
     REAL(rp)  :: E_width
     REAL(rp)  :: PSIP_min
     REAL(rp)  :: PSIp_lim,PSIp_0
     REAL(rp)  :: AMP
     REAL(rp)  :: MARS_AMP_Scale
     !! interpolated E field
     INTEGER 			:: res_double
     INTEGER, DIMENSION(3) 			:: dims 
     !! Dimensions of the KORC vector field. dims=(number of grid 
     !! nodes along \(R\), number of grid nodes along \(\phi\), 
     !! number of grid nodes along \(Z\)).
     INTEGER 			:: dim_1D
     INTEGER 			:: subcycle_E_SC
     REAL(rp)  :: dt_E_SC,Ip_exp,Ip0
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: PSIp
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: PSIp_FS
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE 	:: PSIp3D 
     !! 2-D array for storing the data of the poloidal magnetic flux.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: FLAG2D 
     !! 2-D array defining the simulation domain where pre-computed data exist.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE      :: FLAG3D 
     !! 3-D array defining the simulation domain where pre-computed data exist.
     REAL(rp) 					:: Eo 
     !! Characteristic electric field.
     REAL(rp) 					:: Bo 
     !! Characteristic magnetic field.
     REAL(rp) 					:: Ro 
     !! Radial position of the magnetic axis.
     REAL(rp) 					:: Zo 
     !! \(Z\) position of the magnetic axis.
     LOGICAL 					:: Bfield 
     !! Flag to indicate whether a pre-computed magnetic field will be
     !! used (Bfield=T) or not (Bfield=F).
     LOGICAL 					:: dBfield 
     !! Flag to indicate whether a pre-computed magnetic field will be
     !! used (Bfield=T) or not (Bfield=F).
     LOGICAL 					:: Bflux
     LOGICAL 					:: Bflux3D
     LOGICAL 					:: B1field
     !! Flag to indicate whether a pre-computed poloidal magnetic flux will
     !! be used (Bflux=T) or not (Bflux=F).
     LOGICAL 					:: Efield 
     !! Flag to indicate whether a pre-computed electric field will be used
     !! (Efield=T) or not (Efield=F).
     LOGICAL 					:: Bfield_in_file 
     !! Flag to indicate if a pre-computed magnetic field is in the input file.
     LOGICAL 					:: dBfield_in_file
     LOGICAL 					:: B1field_in_file 
     !! Flag to indicate if a pre-computed magnetic field is in the input file.
     LOGICAL 					:: Bflux_in_file 
     !! Flag to indicate if a pre-computed poloidal magnetic flux is in the
     !! input file.
     LOGICAL 					:: Efield_in_file 
     !! Flag to indicate if a pre-computed electric field is in the input file.
     LOGICAL 					:: axisymmetric_fields 
     !! Flag to indicate if the pre-computed fields are axisymmetric.
     LOGICAL 					:: Dim2x1t
     LOGICAL 					:: E_2x1t,ReInterp_2x1t
     REAL(rp)  :: t0_2x1t
     INTEGER  :: ind0_2x1t,ind_2x1t
#ifdef FIO
     INTEGER  :: isrc
     INTEGER (C_INT)                         :: FIO_B
     !! An M3D-C1 magnetic field.
     INTEGER (C_INT)                         :: FIO_E
     !! An M3D-C1 Electric field.
     INTEGER (C_INT)                         :: FIO_A
     !! An M3D-C1 vector potential.
     real(c_double)  ::  time0,time1
#endif
     REAL(rp)  :: psip_conv
     
  END TYPE FIELDS

end module types
