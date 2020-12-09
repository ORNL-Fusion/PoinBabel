module field
  !! @note Module containing subroutines to initialize externally
  !! generated fields, and to calculate the electric and magnetic
  !! fields when using an analytical model. @endnote
  use types
  use hpc
  use interp
  use PB_HDF5
  USE PB_input

  IMPLICIT NONE

  PUBLIC :: analytical_fields_p,&
       initialize_fields,&
       load_field_data_from_hdf5,&
       load_dim_data_from_hdf5,&
       ALLOCATE_2D_FIELDS_ARRAYS,&
       ALLOCATE_3D_FIELDS_ARRAYS,&
       DEALLOCATE_FIELDS_ARRAYS
  PRIVATE :: ALLOCATE_V_FIELD_2D,&
       ALLOCATE_V_FIELD_3D

CONTAINS

  subroutine analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)
    INTEGER, INTENT(IN)  :: pchunk
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp),  INTENT(IN),DIMENSION(pchunk)      :: Y_R,Y_PHI,Y_Z
    REAL(rp),  INTENT(OUT),DIMENSION(pchunk)     :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(pchunk)                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp),DIMENSION(pchunk)                              :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp),DIMENSION(pchunk)                                :: q
    !! Safety profile \(q(r)\).
    REAL(rp),DIMENSION(pchunk)                             :: cT,sT,rm
    REAL(rp)                            :: Bo,Ro,Zo,lam,qo
    INTEGER                                      :: cc
    !! Particle chunk iterator.

    Bo=F%AB%Bo
    Ro=F%Ro
    Zo=F%Zo
    lam=F%AB%lambda
    qo=F%AB%qo
    
    
    !$OMP SIMD
    !    !$OMP& aligned(cT,sT,cZ,sZ,eta,q,Bp,Bzeta,B_X,B_Y,B_Z, &
    !    !$OMP& Ezeta,E_X,E_Y,E_Z,T_T,T_Z,T_R)
    do cc=1_idef,pchunk

       rm(cc)=sqrt((Y_R(cc)-Ro)**2+(Y_Z(cc)-Zo)**2)
       
       cT(cc)=(Y_R(cc)-Ro)/rm(cc)
       sT(cc)=(Y_Z(cc)-Zo)/rm(cc)


       q(cc) = qo*(1.0_rp + (rm(cc)*rm(cc)/(lam*lam)))
       Bp(cc) = rm(cc)*Ro*Bo/(Y_R(cc)*Y_R(cc)*q(cc))
       Bzeta(cc) = Ro*Bo/Y_R(cc)


       B_R(cc) = - Bp(cc)*sT(cc)
       B_PHI(cc) = -Bzeta(cc)
       B_Z(cc) = Bp(cc)*cT(cc)

    end do
    !$OMP END SIMD

  end subroutine analytical_fields_p

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * * *  SUBROUTINES FOR INITIALIZING FIELDS   * * * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  subroutine initialize_fields(params,F)
    
    !! @note Subroutine that initializes the analytical or externally
    !! calculated electric and magnetic fields. @endnote
    !! In this subroutine we load the parameters of the electric and
    !! magnetic fields from the namelists 'analytical_fields_params' and
    !! 'externalPlasmaModel' in the input file.
    TYPE(KORC_PARAMS), INTENT(INOUT)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(OUT)      :: F
    INTEGER                        :: ii
    !! Iterators for creating mesh for GC model with analytic fields
    INTEGER                        :: kk
    !! Iterators for creating mesh for GC model with analytic fields
    !INTEGER                        :: nR
    !! Number of mesh points in R for grid in GC model of analytical field
    !INTEGER                        :: nZ,nPHI
    !! Number of mesh points in Z for grid in GC model of analytical field
    real(rp)                       :: rm
    !! Minor radius at each position in the grid for
    !! GC model of analytical field
    real(rp)                       :: qr
    !! Safety factor at each position in the grid for
    !! GC model of analytical field
    real(rp)                       :: theta
    !! Poloidal angle at each position in the grid for
    !! GC model of analytical field
    logical :: test
    !integer :: res_double
    real(rp) :: RMAX,RMIN,ZMAX,ZMIN


#ifdef FIO
    F%FIO_B = -1
    F%FIO_E = -1
#endif
    
    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'(/,"* * * * * * * * INITIALIZING FIELDS * * * * * * * *")')
    end if

    if (params%field_model(1:10).eq.'ANALYTICAL') then

       F%AB%Bo = Bo
       F%AB%a = minor_radius
       F%AB%Ro = major_radius
       F%Ro = major_radius
       F%Zo = 0.0_rp
       F%AB%qa = qa
       F%AB%qo = qo
       F%AB%lambda = F%AB%a/SQRT(qa/qo - 1.0_rp)
       F%AB%Bpo = F%AB%lambda*F%AB%Bo/(F%AB%qo*F%AB%Ro)
       F%AB%current_direction = TRIM(current_direction)
       SELECT CASE (TRIM(F%AB%current_direction))
       CASE('PARALLEL')
          F%AB%Bp_sign = 1.0_rp
       CASE('ANTI-PARALLEL')
          F%AB%Bp_sign = -1.0_rp
       CASE DEFAULT
       END SELECT
       F%Bo = F%AB%Bo      
       
       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("ANALYTIC")')
          write(output_unit_write,'("Magnetic field: ",E17.10)') F%Bo
          write(output_unit_write,'("Electric field: ",E17.10)') F%Eo

       end if

       
!    CASE('PSPLINE')
    else if (params%field_model(1:7).eq.'PSPLINE') then

       F%Bfield = Bfield
       F%Bflux = Bflux
       F%axisymmetric_fields = axisymmetric_fields
       F%psip_conv=psip_conv

       
       call load_dim_data_from_hdf5(params,F)
       !sets F%dims for 2D or 3D data

       
       call which_fields_in_file(params,F%Bfield_in_file,F%Bflux_in_file)

       
       if (F%Bflux.AND..NOT.F%Bflux_in_file) then
          write(output_unit_write,'("ERROR: Magnetic flux to be used but no data in file!")')
          call PB_ABORT(18)
       end if

       if (F%Bfield.AND..NOT.F%Bfield_in_file) then
          write(output_unit_write,'("ERROR: Magnetic field to be used but no data in file!")')
          call PB_ABORT(18)
       end if
       

       if (F%axisymmetric_fields) then
          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield,F%Bflux)
       else
          call ALLOCATE_3D_FIELDS_ARRAYS(params,F,F%Bfield)          
       end if
       !allocates 2D or 3D data arrays (fields and spatial)
       
       
       call load_field_data_from_hdf5(params,F)
                   
       if (params%mpi_params%rank .EQ. 0) then

          write(output_unit_write,'("PSPLINE")')
          write(output_unit_write,'("Magnetic field: ",E17.10)') F%Bo
          write(output_unit_write,'("Electric field: ",E17.10)') F%Eo

       end if
       
      
       if (params%mpi_params%rank.EQ.0) then

          if (F%axisymmetric_fields) then
             if (F%Bflux) then
                write(output_unit_write,'("PSIp(r=0)",E17.10)') F%PSIp(F%dims(1)/2,F%dims(3)/2)
                write(output_unit_write,'("BPHI(r=0)",E17.10)') F%Bo
             else
                write(output_unit_write,'("BR(r=0)",E17.10)') F%B_2D%R(F%dims(1)/2,F%dims(3)/2)
                write(output_unit_write,'("BPHI(r=0)",E17.10)') &
                     F%B_2D%PHI(F%dims(1)/2,F%dims(3)/2)
                write(output_unit_write,'("BZ(r=0)",E17.10)') F%B_2D%Z(F%dims(1)/2,F%dims(3)/2)
             end if
          end if


       end if

    end if

    if (params%mpi_params%rank.eq.0) then
       write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * * *",/)')
    end if
       
  end subroutine initialize_fields

  
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Subroutines for getting the fields data from HDF5 files
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  !> @brief Subroutine that loads the size of the arrays having the electric and magnetic field data.
  !! @details All the information of externally calculated fields must be given in a rectangular, equally spaced mesh in the \((R,\phi,Z)\) space of cylindrical coordinates.
  !! If the fields are axisymmetric, then the fields must be in a rectangular mesh on the \(RZ\)-plane.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] F An instance of the KORC derived type FIELDS.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @dims Array containing the size of the mesh with the data of the electric and magnetic fields. dims(1) = dimension along the \(R\) coordinate,
  !! dims(2) = dimension along the \(\phi\) coordinate, and dims(3) = dimension along the \(Z\) coordinate.
  !! @param h5error HDF5 error status.
  !! @param rdamum Temporary variable keeping the read data.
  subroutine load_dim_data_from_hdf5(params,F)
    TYPE(KORC_PARAMS), INTENT(IN)                  :: params
    TYPE(FIELDS), INTENT(INOUT)                    :: F
    CHARACTER(MAX_STRING_LENGTH)                   :: filename
    CHARACTER(MAX_STRING_LENGTH)                   :: gname
    CHARACTER(MAX_STRING_LENGTH)                   :: subgname
    CHARACTER(MAX_STRING_LENGTH)                   :: dset
    INTEGER(HID_T)                                 :: h5file_id
    INTEGER(HID_T)                                 :: group_id
    INTEGER(HID_T)                                 :: subgroup_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE    :: dims
    INTEGER                                        :: h5error
    REAL(rp)                                       :: rdatum

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fopen_f")')
    end if

    if (F%Bflux.OR.F%axisymmetric_fields) then
       dset = "/NR"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(1) = INT(rdatum)

       F%dims(2) = 0

       dset = "/NZ"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(3) = INT(rdatum)
       
    else
       dset = "/NR"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(1) = INT(rdatum)

       dset = "/NPHI"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(2) = INT(rdatum)

       dset = "/NZ"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(3) = INT(rdatum)
    end if


    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine load_dim_data_from_hdf5


  !> @brief Subroutine that queries the HDF5 file what data are present in the HDF5 input file (sanity check).
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param Bfield Logical variable that specifies if the magnetic field is present in the HDF5 file.
  !! @param Efield Logical variable that specifies if the electric field is present in the HDF5 file.
  !! @param Bflux Logical variable that specifies if the poloidal magnetic flux is present in the HDF5 file.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @param h5error HDF5 error status.
  subroutine which_fields_in_file(params,Bfield,Bflux)
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    LOGICAL, INTENT(OUT)               :: Bfield
    LOGICAL, INTENT(OUT)               :: Bflux
    CHARACTER(MAX_STRING_LENGTH)       :: filename
    CHARACTER(MAX_STRING_LENGTH)       :: gname
    CHARACTER(MAX_STRING_LENGTH)       :: subgname
    CHARACTER(MAX_STRING_LENGTH)       :: dset
    INTEGER(HID_T)                     :: h5file_id
    INTEGER(HID_T)                     :: group_id
    INTEGER(HID_T)                     :: subgroup_id
    INTEGER                            :: h5error

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
    end if

    gname = "BR"
    call h5lexists_f(h5file_id,TRIM(gname),Bfield,h5error)


    gname = "PSIp"
    call h5lexists_f(h5file_id,TRIM(gname),Bflux,h5error)


    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine which_fields_in_file


  !> @brief Subroutine that loads the fields data from the HDF5 input file.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @param h5error HDF5 error status.
  subroutine load_field_data_from_hdf5(params,F)
    TYPE(KORC_PARAMS), INTENT(IN)          :: params
    TYPE(FIELDS), INTENT(INOUT)            :: F
    CHARACTER(MAX_STRING_LENGTH)           :: filename
    CHARACTER(MAX_STRING_LENGTH)           :: gname
    CHARACTER(MAX_STRING_LENGTH)           :: subgname
    CHARACTER(MAX_STRING_LENGTH)           :: dset
    INTEGER(HID_T)                         :: h5file_id
    INTEGER(HID_T)                         :: group_id
    INTEGER(HID_T)                         :: subgroup_id
    INTEGER                                :: h5error
    LOGICAL :: Efield

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
    end if

    
    if (((.NOT.F%Bflux).AND.(.NOT.F%axisymmetric_fields)).OR. &
         F%Dim2x1t) then
       dset = "/PHI"
       call load_array_from_hdf5(h5file_id,dset,F%X%PHI)
    end if

    
    dset = "/R"
    call load_array_from_hdf5(h5file_id,dset,F%X%R)

    dset = "/Z"
    call load_array_from_hdf5(h5file_id,dset,F%X%Z)


    dset = '/Bo'
    call load_from_hdf5(h5file_id,dset,F%Bo)

    
    dset = '/Ro'
    call load_from_hdf5(h5file_id,dset,F%Ro)

    dset = '/Zo'
    call load_from_hdf5(h5file_id,dset,F%Zo)

    if ((F%Bflux.OR.F%axisymmetric_fields).AND.(.NOT.F%Dim2x1t)) then

       dset = "/FLAG"
       call load_array_from_hdf5(h5file_id,dset,F%FLAG2D)

    else
       dset = "/FLAG"
       call load_array_from_hdf5(h5file_id,dset,F%FLAG3D)
    end if
    
    if (F%Bflux) then

             
       dset = "/PSIp"
       gname = 'PSIp'

       call h5lexists_f(h5file_id,TRIM(gname),Efield,h5error)

       if (Efield) then
          call load_array_from_hdf5(h5file_id,dset,F%PSIp)
       else
          F%PSIp = 0.0_rp
       end if



    end if


    
    if (F%Bfield) then
       if (F%axisymmetric_fields) then
          dset = "/BR"
          call load_array_from_hdf5(h5file_id,dset,F%B_2D%R)

          dset = "/BPHI"
          call load_array_from_hdf5(h5file_id,dset,F%B_2D%PHI)

          dset = "/BZ"
          call load_array_from_hdf5(h5file_id,dset,F%B_2D%Z)
       else
          dset = "/BR"
          call load_array_from_hdf5(h5file_id,dset,F%B_3D%R)

          !write(6,*) 'BR(25,1,:)',F%B_3D%R(25,1,:)

          dset = "/BPHI"
          call load_array_from_hdf5(h5file_id,dset,F%B_3D%PHI)

          dset = "/BZ"
          call load_array_from_hdf5(h5file_id,dset,F%B_3D%Z)
       end if
    end if



    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine load_field_data_from_hdf5



  subroutine ALLOCATE_2D_FIELDS_ARRAYS(params,F,bfield,bflux)
    !! @note Subroutine that allocates the variables keeping the axisymmetric
    !! fields data. @endnote
    TYPE (KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)    :: F
    !! An instance of the KORC derived type FIELDS. In this variable we keep
    !! the loaded data.
    LOGICAL, INTENT(IN)            :: bfield
    LOGICAL, INTENT(IN)            :: bflux
    !! Logical variable that specifies if the variables that keep the poloidal
    !! magnetic flux data is allocated (bflux=T) or not (bflux=F).


    if (bfield.and.(.not.ALLOCATED(F%B_2D%R))) then
       call ALLOCATE_V_FIELD_2D(F%B_2D,F%dims)


    end if

    if (bflux.and.(.not.ALLOCATED(F%PSIp))) then
       ALLOCATE(F%PSIp(F%dims(1),F%dims(3)))
       F%PSIp=0._rp
    end if


    if (.NOT.ALLOCATED(F%FLAG2D)) ALLOCATE(F%FLAG2D(F%dims(1),F%dims(3)))

    if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
    if (.NOT.ALLOCATED(F%X%Z)) ALLOCATE(F%X%Z(F%dims(3)))
  end subroutine ALLOCATE_2D_FIELDS_ARRAYS


  !> @brief Subroutine that allocates the variables keeping the 3-D fields data.
  !!
  !! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
  !! @param[in] bfield Logical variable that specifies if the variables that keep the magnetic field data is allocated (bfield=T) or not (bfield=F).
  !! @param[in] efield Logical variable that specifies if the variables that keep the electric field data is allocated (efield=T) or not (efield=F).
  subroutine ALLOCATE_3D_FIELDS_ARRAYS(params,F,bfield)
    TYPE (KORC_PARAMS), INTENT(IN) 	:: params
    TYPE(FIELDS), INTENT(INOUT)    :: F
    LOGICAL, INTENT(IN)            :: bfield

    if (bfield) then
       call ALLOCATE_V_FIELD_3D(F%B_3D,F%dims)       
    end if


    if (.NOT.ALLOCATED(F%FLAG3D)) ALLOCATE(F%FLAG3D(F%dims(1),F%dims(2),F%dims(3)))

    if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
    if (.NOT.ALLOCATED(F%X%PHI)) ALLOCATE(F%X%PHI(F%dims(2)))
    if (.NOT.ALLOCATED(F%X%Z)) ALLOCATE(F%X%Z(F%dims(3)))
  end subroutine ALLOCATE_3D_FIELDS_ARRAYS


  !> @brief Subroutine that allocates the cylindrical components of an axisymmetric field.
  !!
  !! @param[in,out] F Vector field to be allocated.
  !! @param[in] dims Dimension of the mesh containing the field data.
  subroutine ALLOCATE_V_FIELD_2D(F,dims)
    TYPE(V_FIELD_2D), INTENT(INOUT)    :: F
    INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(3)))
  end subroutine ALLOCATE_V_FIELD_2D


  !> @brief Subroutine that allocates the cylindrical components of a 3-D field.
  !!
  !! @param[in,out] F Vector field to be allocated.
  !! @param[in] dims Dimension of the mesh containing the field data.
  subroutine ALLOCATE_V_FIELD_3D(F,dims)
    TYPE(V_FIELD_3D), INTENT(INOUT)    :: F
    INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(2),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(2),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(2),dims(3)))
  end subroutine ALLOCATE_V_FIELD_3D

  !> @brief Subroutine that deallocates all the variables of the electric and magnetic fields.
  !!
  !! @param[in,out] F An instance of the KORC derived type FIELDS.
  subroutine DEALLOCATE_FIELDS_ARRAYS(F)
    TYPE(FIELDS), INTENT(INOUT) :: F

    if (ALLOCATED(F%PSIp)) DEALLOCATE(F%PSIp)

    if (ALLOCATED(F%B_2D%R)) DEALLOCATE(F%B_2D%R)
    if (ALLOCATED(F%B_2D%PHI)) DEALLOCATE(F%B_2D%PHI)
    if (ALLOCATED(F%B_2D%Z)) DEALLOCATE(F%B_2D%Z)

    if (ALLOCATED(F%B_3D%R)) DEALLOCATE(F%B_3D%R)
    if (ALLOCATED(F%B_3D%PHI)) DEALLOCATE(F%B_3D%PHI)
    if (ALLOCATED(F%B_3D%Z)) DEALLOCATE(F%B_3D%Z)


    if (ALLOCATED(F%X%R)) DEALLOCATE(F%X%R)
    if (ALLOCATED(F%X%PHI)) DEALLOCATE(F%X%PHI)
    if (ALLOCATED(F%X%Z)) DEALLOCATE(F%X%Z)

    if (ALLOCATED(F%FLAG2D)) DEALLOCATE(F%FLAG2D)
    if (ALLOCATED(F%FLAG3D)) DEALLOCATE(F%FLAG3D)
  end subroutine DEALLOCATE_FIELDS_ARRAYS
end module field
