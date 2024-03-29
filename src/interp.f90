module PB_interp
  !! @note Module containing functions and subroutines for performing
  !! interpolations using the PSPLINE library. @endnote
  !! For a detailed documentation of the PSPLINE library we refer the
  !! user to "https://w3.pppl.gov/ntcc/PSPLINE/".
  use PB_types
  use PB_hpc
  use PB_constants

#ifdef PSPLINE
  use EZspline_obj	! psplines module
  use EZspline          ! psplines module
#endif

#ifdef FIO
  use PB_fio
#endif

  !$ use OMP_LIB

  IMPLICIT NONE

#ifdef PSPLINE

  TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for
     !! cylindrical components of vector fields
     !! \(\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi +
     !! F_Z\hat{e}_Z\). Real precision of 8 bytes. @endnote
     TYPE(EZspline3)    :: A     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline3)    :: R
     !! Interpolant of \(F_R(R,\phi,Z)\).
     TYPE(EZspline3)    :: PHI
     !! Interpolant of \(F_\phi(R,\phi,Z)\).
     TYPE(EZspline3)    :: Z
     !! Interpolant of \(F_Z(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.

     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both
     !! ends of the \(\phi\) direction.

     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_3D_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of vector fields \(\mathbf{F}(R,Z) =
     !! F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline2)    :: R
     !! Interpolant of \(F_R(R,Z)\).
     TYPE(EZspline2)    :: PHI
     !! Interpolant of \(F_\phi(R,Z)\).
     TYPE(EZspline2)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2D_FIELDS_INTERPOLANT



  TYPE, PRIVATE :: KORC_INTERPOLANT_DOMAIN
     !! @note Derived type containing 2-D and 3-D arrays with the information of
     !! the spatial domain where the fields and profiles are known.
     !! This info is used for detecting when a particle is lost, and therefore not
     !! followed anymore. @endnote
     INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE      :: FLAG1D
     !! 2-D array with info of the spatial domain where the axisymmetric fields
     !! and plasma profiles are known.
     INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE      :: FLAG2D
     !! 2-D array with info of the spatial domain where the axisymmetric fields
     !! and plasma profiles are known.
     INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE    :: FLAG3D
     !! 3-D array with info of the spatial domain where the 3-D fields and plasma
     !! profiles are known.

     REAL(rp)                                          :: Ro
     !! Smaller radial position of the fields and profiles domain.
     REAL(rp)                                          :: Zo
     !! Smaller vertical position of the fields and profiles domain
     REAL(rp)                                          :: To

     REAL(rp)                                          :: Drm
     REAL(rp)                                          :: DPSIP
     REAL(rp)                                          :: DR
     !! Separation between grid points along the radial direction.
     REAL(rp)                                          :: DPHI  !
     ! Separation between grid points along the azimuthal direction.
     REAL(rp)                                          :: DT  !
     ! Separation between grid points along the azimuthal direction.
     REAL(rp)                                          :: DZ
     !! Separation between grid points along the vertical direction.
  END TYPE KORC_INTERPOLANT_DOMAIN

  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: b1Refield_2d
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: b1Imfield_2d
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_3d
  !! An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE         :: fields_domain
  !! An instance of KORC_INTERPOLANT_DOMAIN used for interpolating fields.
  INTEGER                                        :: ezerr
  !! Error status during PSPLINE interpolations.

#endif


CONTAINS

#ifdef PSPLINE
  subroutine initialize_fields_interpolant(params,F)
    !! @note Subroutine that initializes fields interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of fields in the KORC-dervied-type variable F.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)       :: F
    !! An instance of KORC's derived type FIELDS containing all the information
    !! about the fields used in the simulation.
    !! See [[korc_types]] and [[korc_fields]].
    integer :: ii,jj

    if ((params%field_model .EQ. 'PSPLINE').or. &
         (params%field_model .EQ. 'MARS')) then

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * INITIALIZING FIELDS INTERPOLANT * * * *")')
       end if

       ! * * * * * * * * MAGNETIC FIELD * * * * * * * * !
       if (F%Bflux) then

             write(output_unit_write,*) '2D poloidal flux function'

             if (EZspline_allocated(bfield_2d%A)) &
             call Ezspline_free(bfield_2d%A, ezerr)

             bfield_2d%NR = F%dims(1)
             bfield_2d%NZ = F%dims(3)

             ! Initializing poloidal flux interpolant
             call EZspline_init(bfield_2d%A,bfield_2d%NR,bfield_2d%NZ, &
             bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             bfield_2d%A%x1 = F%X%R
             bfield_2d%A%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z
             !write(output_unit_write,'("PSIp",E17.10)') F%PSIp3D(:,F%ind_2x1t,:)

             if (F%Dim2x1t) THEN

                call EZspline_setup(bfield_2d%A, F%PSIp3D(:,F%ind_2x1t,:), ezerr, .TRUE.)
                call EZspline_error(ezerr)

                !write(output_unit_write,'("bfield_2d%A: ",E17.10)') bfield_2d%A%fspl(1,:,:)

                if (.not.ALLOCATED(fields_domain%FLAG2D)) &
                ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))

                fields_domain%FLAG2D = F%FLAG3D(:,F%ind_2x1t,:)

             else

                call EZspline_setup(bfield_2d%A, F%PSIp, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                !write(output_unit_write,'("bfield_2d%A: ",E17.10)') bfield_2d%A%fspl(1,:,:)

                if (.not.ALLOCATED(fields_domain%FLAG2D)) &
                ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))

                fields_domain%FLAG2D = F%FLAG2D

             end if

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))

       end if

       if(F%B1field) then

             write(output_unit_write,*) '2D n=1 MARS magnetic fields'

             b1Refield_2d%NR = F%dims(1)
             b1Refield_2d%NZ = F%dims(3)

             ! Initializing BR1Re interpolant
             call EZspline_init(b1Refield_2d%R,b1Refield_2d%NR, &
                  b1Refield_2d%NZ,b1Refield_2d%BCSR,b1Refield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2d%R%x1 = F%X%R
             b1Refield_2d%R%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2d%R, F%B1Re_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Re_2D%R(:,200))
             !write(6,*) 'B1Re_2D%R',F%B1Re_2D%R(:,200)*params%cpp%Bo

             ! Initializing BPHI1Re interpolant
             call EZspline_init(b1Refield_2d%PHI,b1Refield_2d%NR, &
                  b1Refield_2d%NZ,b1Refield_2d%BCSR,b1Refield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2d%PHI%x1 = F%X%R
             b1Refield_2d%PHI%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2d%PHI, F%B1Re_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing BZ1_Re interpolant
             call EZspline_init(b1Refield_2d%Z,b1Refield_2d%NR, &
                  b1Refield_2d%NZ,b1Refield_2d%BCSR,b1Refield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2d%Z%x1 = F%X%R
             b1Refield_2d%Z%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2d%Z, F%B1Re_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             b1Imfield_2d%NR = F%dims(1)
             b1Imfield_2d%NZ = F%dims(3)

             ! Initializing BR1RIm interpolant
             call EZspline_init(b1Imfield_2d%R,b1Imfield_2d%NR, &
                  b1Imfield_2d%NZ,b1Imfield_2d%BCSR,b1Imfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2d%R%x1 = F%X%R
             b1Imfield_2d%R%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2d%R, F%B1Im_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing BPHI1Im interpolant
             call EZspline_init(b1Imfield_2d%PHI,b1Imfield_2d%NR, &
                  b1Imfield_2d%NZ,b1Imfield_2d%BCSR,b1Imfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2d%PHI%x1 = F%X%R
             b1Imfield_2d%PHI%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2d%PHI, F%B1Im_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)


             ! Initializing BZ1_Im interpolant
             call EZspline_init(b1Imfield_2d%Z,b1Imfield_2d%NR, &
                  b1Imfield_2d%NZ,b1Imfield_2d%BCSR,b1Imfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2d%Z%x1 = F%X%R
             b1Imfield_2d%Z%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2d%Z, F%B1Im_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)


       end if

       if (F%Bfield) then
          if (F%axisymmetric_fields) then

             write(output_unit_write,*) '2D magnetic field'

             bfield_2d%NR = F%dims(1)
             bfield_2d%NZ = F%dims(3)

             ! Initializing R component
             call EZspline_init(bfield_2d%R,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%R%x1 = F%X%R
             bfield_2d%R%x2 = F%X%Z

             call EZspline_setup(bfield_2d%R, F%B_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component
             call EZspline_init(bfield_2d%PHI,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%PHI%x1 = F%X%R
             bfield_2d%PHI%x2 = F%X%Z

             call EZspline_setup(bfield_2d%PHI, F%B_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Z component
             call EZspline_init(bfield_2d%Z,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%Z%x1 = F%X%R
             bfield_2d%Z%x2 = F%X%Z

             call EZspline_setup(bfield_2d%Z, F%B_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             if (.not.ALLOCATED(fields_domain%FLAG2D)) &
                  ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))

             fields_domain%FLAG2D = F%FLAG2D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
          else

             write(output_unit_write,*) '3D magnetic field'


             if (stel_sym) then
                bfield_3d%BCSPHI = (/ 0, 0 /)
             endif

             bfield_3d%NR = F%dims(1)
             bfield_3d%NPHI = F%dims(2)
             bfield_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(bfield_3d%R, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%R%x1 = F%X%R
             bfield_3d%R%x2 = F%X%PHI
             bfield_3d%R%x3 = F%X%Z

             call EZspline_setup(bfield_3d%R, F%B_3D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component of interpolant
             call EZspline_init(bfield_3d%PHI, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%PHI%x1 = F%X%R
             bfield_3d%PHI%x2 = F%X%PHI
             bfield_3d%PHI%x3 = F%X%Z

             call EZspline_setup(bfield_3d%PHI, F%B_3D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(output_unit_write,*) bfield_3d%PHI%x2


             ! Initializing Z component of interpolant
             call EZspline_init(bfield_3d%Z, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%Z%x1 = F%X%R
             bfield_3d%Z%x2 = F%X%PHI
             bfield_3d%Z%x3 = F%X%Z

             call EZspline_setup(bfield_3d%Z, F%B_3D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ALLOCATE(fields_domain%FLAG3D(bfield_3d%NR,bfield_3d%NPHI, &
                  bfield_3d%NZ))
             fields_domain%FLAG3D = F%FLAG3D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             if (.not.F%stel_sym) then
                fields_domain%DPHI = 2.0_rp*C_PI/(bfield_3d%NPHI-1)
             else
                fields_domain%DPHI = C_PI/F%nsymm/bfield_3d%NPHI
             endif
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
          end if
       end if


       fields_domain%Ro = F%X%R(1)
       fields_domain%Zo = F%X%Z(1)

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * * * INTERPOLANT INITIALIZED * * * * * *",/)')
       end if
    else if (params%field_model(1:10) .EQ. 'ANALYTICAL') then
       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * USING ANALYTICAL MAGNETIC FIELD * * * *",/)')
       end if
    end if
  end subroutine initialize_fields_interpolant

  subroutine check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag)

    use omp_lib

    !! @note Subrotuine that checks if particles in the simulation are within
    !! the spatial domain where interpolants and fields are known. @endnote
    !! External fields and interpolants can have different spatial domains where
    !! they are defined. Therefore, it is necessary to
    !! check if a given particle has left these spatial domains to
    !! stop following it, otherwise this will cause an error in the simulation.
    INTEGER, INTENT(IN)  :: pchunk
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp), DIMENSION(pchunk),  INTENT(IN)      :: Y_R,Y_PHI,Y_Z
    INTEGER(is), DIMENSION(pchunk), INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed in the
    !! simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Species iterator.
    INTEGER             :: thread_num

    thread_num = OMP_GET_THREAD_NUM()

    if (ALLOCATED(fields_domain%FLAG3D)) then

       !$OMP SIMD
       !       !$OMP&  aligned(IR,IPHI,IZ)
       do pp=1_idef,pchunk

          IR = INT(FLOOR((Y_R(pp)  - fields_domain%Ro + &
               0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IPHI = INT(FLOOR((Y_PHI(pp)  + 0.5_rp*fields_domain%DPHI)/ &
               fields_domain%DPHI) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(fields_domain%Zo) + &
               0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

          !write(6,*) thread_num,'flag:',shape(fields_domain%FLAG3D)

          !if (((IR.GT.bfield_3d%NR).OR.(IZ.GT.bfield_3d%NZ))) then



          !write(6,*) thread_num,'R,PHI,Z',Y_R,Y_PHI,Y_Z
          !write(6,*) thread_num,'query:',IR,IPHI,IZ
          !write(6,*) thread_num,'total:',bfield_3d%NR,bfield_3d%NPHI,bfield_3d%NZ

          !end if


          if (((IR.GT.bfield_3d%NR).OR.(IZ.GT.bfield_3d%NZ)).OR. &
               ((IR.LT.1).OR.(IZ.LT.1)).OR. &
               (fields_domain%FLAG3D(IR,IPHI,IZ).NE.1_is)) then
             flag(pp) = 0_is

          end if
       end do
       !$OMP END SIMD

    else if (ALLOCATED(fields_domain%FLAG2D)) then
       !$OMP SIMD
       !       !$OMP& aligned(IR,IZ)
       do pp=1_idef,pchunk

          IR = INT(FLOOR((Y_R(pp)  - fields_domain%Ro + &
               0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(fields_domain%Zo) + &
               0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)



          if ((fields_domain%FLAG2D(IR,IZ).NE.1_is).OR. &
               ((IR.GT.bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ))) then
             flag(pp) = 0_is

          end if
       end do
       !$OMP END SIMD
       !       write(output_unit_write,'("Shit''s not fucked.")')
    end if
  end subroutine check_if_in_fields_domain_p


  subroutine interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)
    INTEGER, INTENT(IN)  :: pchunk
    TYPE(FIELDS), INTENT(IN)                               :: F
    REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_R,B_PHI,B_Z
    INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flagCon
    !  INTEGER(ip) :: ezerr

    call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flagCon)

    call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z, &
         pchunk,Y_R,Y_Z,B_R,B_PHI,B_Z,ezerr)
    call EZspline_error(ezerr)


  end subroutine interp_2DB_p

  subroutine interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)
    INTEGER, INTENT(IN)  :: pchunk
    TYPE(FIELDS), INTENT(IN)                               :: F
    REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(pchunk)   :: YPHI,YZ
    REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_R,B_PHI,B_Z
    INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flagCon
    !  INTEGER(ip) :: ezerr
    INTEGER(ip)                                            :: pp
    REAL(rp) :: nsymm

    if (F%stel_sym) then

       nsymm=F%nsymm

       YPHI=Y_PHI
       YZ=Y_Z

       do pp=1_idef,pchunk
          YPHI(pp)=modulo(YPHI(pp),2*C_PI/nsymm)
          if (YPHI(pp).ge.C_PI/nsymm) then
             YPHI(pp)=2*C_PI/nsymm-YPHI(pp)
             YZ(pp)=-YZ(pp)
          end if
       end do

    else
       YPHI=Y_PHI
       YZ=Y_Z
    endif

    call check_if_in_fields_domain_p(pchunk,F,Y_R,YPHI,YZ,flagCon)

    call EZspline_interp(bfield_3d%R,bfield_3d%PHI,bfield_3d%Z, &
         pchunk,Y_R,YPHI,YZ,B_R,B_PHI,B_Z,ezerr)
    call EZspline_error(ezerr)

    if (F%stel_sym) then

       YPHI=Y_PHI

       do pp=1_idef,pchunk
          YPHI(pp)=modulo(YPHI(pp),2*C_PI/nsymm)
          if (YPHI(pp).ge.C_PI/nsymm) then
             B_R(pp)=-B_R(pp)
          end if
       end do

    endif

    !write(6,*) 'R,PHI,Z:0',Y_R,Y_PHI,Y_Z
    !write(6,*) 'R,PHI,Z:1',Y_R,YPHI,YZ
    !write(6,*) 'BR,BPHI,BZ',B_R,B_PHI,B_Z

  end subroutine interp_3DB_p


  subroutine calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,flagCon)
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
    TYPE(FIELDS), INTENT(IN)                               :: F
    REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
    INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flagCon
    REAL(rp), DIMENSION(pchunk,2)  :: A
    REAL(rp)  :: psip_conv
    INTEGER :: cc

    psip_conv=F%psip_conv

    call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flagCon)

    ! FR = (dA/dZ)/R
    call EZspline_gradient(bfield_2d%A, pchunk, Y_R, Y_Z, &
         A, ezerr)
    call EZspline_error(ezerr)

    !write(output_unit_write,'("dPSIp/dR: ",E17.10)') A(:,1)
    !write(output_unit_write,'("dPSIp/dZ: ",E17.10)') A(:,2)
    !write(output_unit_write,'("Y_R: ",E17.10)') Y_R

    !$OMP SIMD
    do cc=1_idef,pchunk
       B_R(cc) = psip_conv*A(cc,2)/Y_R(cc)

       ! FPHI = Fo*Ro/R
       B_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)

       ! FR = -(dA/dR)/R

       !     write(output_unit_write,'("R*B_Z: ",E17.10)') B_Z(1)

       B_Z(cc) = -psip_conv*A(cc,1)/Y_R(cc)
    end do
    !$OMP END SIMD

  end subroutine calculate_magnetic_field_p

  subroutine interp_3Dmars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
       flag_cache)
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    INTEGER, INTENT(IN)  :: pchunk
    TYPE(FIELDS), INTENT(IN)                               :: F
    REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(pchunk)   :: B0_R,B0_PHI,B0_Z
    REAL(rp),DIMENSION(pchunk)   :: B1_R,B1_PHI,B1_Z
    REAL(rp),DIMENSION(pchunk)   :: B1Re_R,B1Re_PHI,B1Re_Z
    REAL(rp),DIMENSION(pchunk)   :: B1Im_R,B1Im_PHI,B1Im_Z
    REAL(rp),DIMENSION(pchunk)   :: PSIp
    REAL(rp),DIMENSION(pchunk)   :: cP,sP
    REAL(rp), DIMENSION(pchunk,3)  :: A
    !  INTEGER(ip) :: ezerr
    INTEGER                                      :: cc
    !! Particle chunk iterator.
    INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
    REAL(rp) :: psip_conv
    REAL(rp) :: amp

    psip_conv=F%psip_conv
    amp=F%AMP

    call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

    call EZspline_interp(bfield_2d%A,b1Refield_2d%R,b1Refield_2d%PHI, &
         b1Refield_2d%Z,b1Imfield_2d%R,b1Imfield_2d%PHI,b1Imfield_2d%Z, &
         pchunk,Y_R,Y_Z,A,B1Re_R,B1Re_PHI,B1Re_Z,B1Im_R,B1Im_PHI,B1Im_Z,ezerr)
    call EZspline_error(ezerr)

    !$OMP SIMD
    do cc=1_idef,pchunk
       PSIp(cc)=A(cc,1)

       B0_R(cc) = psip_conv*A(cc,3)/Y_R(cc)
       B0_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)
       B0_Z(cc) = -psip_conv*A(cc,2)/Y_R(cc)

       cP(cc)=cos(Y_PHI(cc))
       sP(cc)=sin(Y_PHI(cc))

       B1_R(cc) = amp*(B1Re_R(cc)*cP(cc)-B1Im_R(cc)*sP(cc))
       B1_PHI(cc) = amp*(B1Re_PHI(cc)*cP(cc)-B1Im_PHI(cc)*sP(cc))
       B1_Z(cc) = amp*(B1Re_Z(cc)*cP(cc)-B1Im_Z(cc)*sP(cc))

       B_R(cc) = B0_R(cc)+B1_R(cc)
       B_PHI(cc) = B0_PHI(cc)+B1_PHI(cc)
       B_Z(cc) = B0_Z(cc)+B1_Z(cc)

    end do
    !$OMP END SIMD

    !write(6,*) '(R,PHI,Z)',Y_R,Y_PHI,Y_Z
    !write(6,*) 'psi',PSIp*params%cpp%Bo*params%cpp%length**2
    !write(6,*) 'dpsidR',A(:,2)*params%cpp%Bo*params%cpp%length
    !write(6,*) 'dpsidZ',A(:,3)*params%cpp%Bo*params%cpp%length
    !write(6,*) 'B0',B0_R,B0_PHI,B0_Z
    !write(6,*) 'AMP',amp
    !write(6,*) 'B1Re',B1Re_R,B1Re_PHI,B1Re_Z
    !write(6,*) 'B1Im',B1Im_R,B1Im_PHI,B1Im_Z
    !write(6,*) 'B1',B1_R,B1_PHI,B1_Z
    !write(6,*) 'B',B_R,B_PHI,B_Z


  end subroutine interp_3Dmars_p

  !> @brief Subroutine that frees memory allocated for PSPLINE interpolants.
  !!
  !! @param[in] params Core KORC simulation parameters.
  subroutine finalize_interpolants(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params

    if (params%field_model(1:8) .EQ. 'PSPLINE') then
       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * FINALIZING FIELD INTERPOLANT * * * *")')
       end if

       if (EZspline_allocated(bfield_3d%R)) call Ezspline_free(bfield_3d%R, ezerr)
       if (EZspline_allocated(bfield_3d%PHI)) &
            call Ezspline_free(bfield_3d%PHI,ezerr)
       if (EZspline_allocated(bfield_3d%Z)) call Ezspline_free(bfield_3d%Z, ezerr)
       if (EZspline_allocated(bfield_2d%A)) call Ezspline_free(bfield_2d%A, ezerr)
       if (EZspline_allocated(bfield_2d%R)) call Ezspline_free(bfield_2d%R, ezerr)
       if (EZspline_allocated(bfield_2d%PHI)) &
            call Ezspline_free(bfield_2d%PHI,ezerr)
       if (EZspline_allocated(bfield_2d%Z)) call Ezspline_free(bfield_2d%Z, ezerr)



       if (ALLOCATED(fields_domain%FLAG2D)) DEALLOCATE(fields_domain%FLAG2D)
       if (ALLOCATED(fields_domain%FLAG3D)) DEALLOCATE(fields_domain%FLAG3D)

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * FIELD INTERPOLANT FINALIZED * * * *")')
       end if
    end if
  end subroutine finalize_interpolants
#endif

#ifdef FIO
  subroutine get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,flag,hint)
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: B_R,B_PHI,B_Z
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp,pchunk
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Btmp

    pchunk=params%pchunk

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)
          x(2) = Y_PHI(pp)
          x(3) = Y_Z(pp)

          !             prtcls%hint(pp)=c_null_ptr

          status = fio_eval_field(F%FIO_B, x(1),                      &
               Btmp(1),hint(pp))

          if (status .eq. FIO_SUCCESS) then
             B_R(pp)=Btmp(1)
             B_PHI(pp)=Btmp(2)
             B_Z(pp)=Btmp(3)
          else if (status .eq. FIO_NO_DATA) then
             flag(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             flag(pp) = 0_is
          end if

       end if
    end do

  end subroutine get_fio_magnetic_fields_p
#endif


end module PB_interp
