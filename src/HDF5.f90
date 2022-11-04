!! @note KORC module containing subroutines to read and write data in HDF5
!! files. @endnote
!! This module contains interfaces to use the HDF5 library in a more friendly
!! way. This module is intended to help developers to create new I/O
!! subroutines without having to deal with the sometimes cumbersome details
!! of the HDF5 API.
module PB_HDF5
  use PB_hpc
  use PB_types
  use PB_constants
  use HDF5

  IMPLICIT NONE

  INTEGER(HID_T), PRIVATE 	:: KORC_HDF5_REAL
  !! HDF5 real precision data type to be used in the simulation.
  INTEGER(SIZE_T), PRIVATE 	:: rp_hdf5
  !! Size of the HDF5 real precision data type used in the simulation.

  INTERFACE load_from_hdf5
     !! @note Fortran interface to subroutines loading a real or integer
     !! value from HDF5 files. @endnote
     module procedure iload_from_hdf5, rload_from_hdf5
  END INTERFACE load_from_hdf5


  INTERFACE load_array_from_hdf5
     !! @note Fortran interface to subroutines loading 2-D and 3-D arrays
     !! of real values from HDF5 files.
     module procedure rload_1d_array_from_hdf5, rload_3d_array_from_hdf5, rload_2d_array_from_hdf5
  END INTERFACE load_array_from_hdf5


  INTERFACE save_to_hdf5
     !! @note Fortran interface to subroutines saving real or integer
     !! values to HDF5 files.
     module procedure i1save_to_hdf5,i2save_to_hdf5,i4save_to_hdf5,i8save_to_hdf5,rsave_to_hdf5
  END INTERFACE save_to_hdf5

  !! @note Fortran interface to subroutines saving real and integer
  !! values to HDF5 files.
  INTERFACE save_1d_array_to_hdf5
     module procedure isave_1d_array_to_hdf5,rsave_1d_array_to_hdf5
  END INTERFACE save_1d_array_to_hdf5

  !! @note Fortran interface to subroutines saving 2-D arrays of real values to HDF5 files.
  !! @todo To code the corresponding subroutines for saving integer 2-D arrays.
  INTERFACE save_2d_array_to_hdf5
     module procedure rsave_2d_array_to_hdf5
  END INTERFACE save_2d_array_to_hdf5

  !! @note Fortran interface to subroutines saving 3-D arrays of real values to HDF5 files.
  !! @todo To include the corresponding subroutines for saving arrays of integers.
  INTERFACE save_3d_array_to_hdf5
     module procedure rsave_3d_array_to_hdf5
  END INTERFACE save_3d_array_to_hdf5

  !! @note Fortran interface to subroutines saving 1-D, 2-D or 3-D arrays of real values to HDF5 files.
  !! @todo To include the corresponding subroutines for saving arrays of integers.
  INTERFACE save_array_to_hdf5
     module procedure isave_1d_array_to_hdf5,rsave_1d_array_to_hdf5,rsave_2d_array_to_hdf5,rsave_3d_array_to_hdf5
  END INTERFACE save_array_to_hdf5

  PRIVATE :: rsave_to_hdf5,&
       isave_1d_array_to_hdf5,&
       rsave_1d_array_to_hdf5,&
       rsave_2d_array_to_hdf5,&
       iload_from_hdf5,&
       rload_from_hdf5,&
       rload_1d_array_from_hdf5,&
       rload_3d_array_from_hdf5,&
       rload_2d_array_from_hdf5,&
       i1save_to_hdf5,&
       i2save_to_hdf5,&
       i4save_to_hdf5,&
       i8save_to_hdf5

  PUBLIC :: initialize_HDF5,&
       finalize_HDF5,&
       save_simulation_parameters,&
       save_to_hdf5,&
       save_1d_array_to_hdf5,&
       save_2d_array_to_hdf5,&
       load_from_hdf5,&
       load_array_from_hdf5,&
       save_string_parameter

CONTAINS

  !! @note Initialization of HDF5 library.
  !!
  !! @param h5error HDF5 error status.
  subroutine initialize_HDF5()
    INTEGER :: h5error  ! Error flag
    call h5open_f(h5error)

#ifdef HDF5_DOUBLE_PRESICION
    call h5tcopy_f(H5T_NATIVE_DOUBLE, KORC_HDF5_REAL, h5error)
#elif HDF5_SINGLE_PRESICION
    call h5tcopy_f(H5T_NATIVE_REAL, KORC_HDF5_REAL, h5error)
#endif

    call h5tget_size_f(KORC_HDF5_REAL, rp_hdf5, h5error)
  end subroutine initialize_HDF5

  !! @note Finalization of HDF5 library.
  !!
  !! @param h5error HDF5 error status.
  subroutine finalize_HDF5()
    INTEGER :: h5error  ! Error flag
    call h5close_f(h5error)
  end subroutine finalize_HDF5

  !! @note Subroutine to load an integer datum from an HDF5 file.
  !!
  !! @todo Implement the reading of the attribute of idatum.
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[out] idatum Integer datum read from HDF5 file.
  !! @param[out] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  subroutine iload_from_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 				:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 		:: dset
    INTEGER, INTENT(OUT) 				:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(OUT) :: attr
    CHARACTER(4) 					:: aname = "Info"
    INTEGER(HID_T) 					:: dset_id
    INTEGER(HID_T) 					:: dspace_id
    INTEGER(HID_T) 					:: aspace_id
    INTEGER(HID_T) 					:: attr_id
    INTEGER(HID_T) 					:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 			:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 			:: adims = (/1/)
    INTEGER 						:: h5error

    ! * * * Read datum from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: iload_from_hdf5 &
            --> h5dopen_f")')
    end if

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, idatum, dims, h5error)

    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: iload_from_hdf5 &
            --> h5dread_f")')
    end if

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: iload_from_hdf5 &
            --> h5dclose_f")')
    end if

    if (PRESENT(attr)) then
       ! * * * Read attribute from file * * *

       ! * * * Read attribute from file * * *
    end if

    ! * * * Read datum from file * * *
  end subroutine iload_from_hdf5

  !! @note Subroutine to load a real datum from an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[out] rdatum Real datum read from HDF5 file and casted to
  !! KORC's real precision type.
  !! @param[out] attr Attribute of datum read from HDF5 file.
  !! @param raw_datum Datum read from HDF5 file.
  !! @param aname Name of rdatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attribute of rdatum.
  subroutine rload_from_hdf5(h5file_id,dset,rdatum,attr)
    INTEGER(HID_T), INTENT(IN) 				:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 		:: dset
    REAL(rp), INTENT(OUT) 				:: rdatum
    REAL 						:: raw_datum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(OUT) :: attr
    CHARACTER(4) 					:: aname = "Info"
    INTEGER(HID_T) 					:: dset_id
    INTEGER(HID_T) 					:: dspace_id
    INTEGER(HID_T) 					:: aspace_id
    INTEGER(HID_T) 					:: attr_id
    INTEGER(HID_T) 					:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 			:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 			:: adims = (/1/)
    INTEGER 						:: h5error

    ! * * * Read datum from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 &
            --> h5dopen_f")')
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_datum, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 &
            --> h5dread_f")')
    end if
    rdatum = REAL(raw_datum,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 &
            --> h5dclose_f")')
    end if

    if (PRESENT(attr)) then
       ! * * * Read attribute from file * * *

       ! * * * Read attribute from file * * *
    end if

    ! * * * Read datum from file * * *
  end subroutine rload_from_hdf5

  !! @note Subroutine to load a 1-D array of reals from an HDF5 file.
  !! @details The dimension of the 1-D array rdata is determined by the
  !! input-output array rdata.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[out] rdata 1-D array of real values read from HDF5 file and
  !! casted to KORC's real precision type.
  !! @param[out] attr 1-D array of attributes of rdata.
  !! @param raw_data 1-D array read from HDF5 file.
  !! @param aname Name of rdata attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attributes of rdata.
  subroutine rload_1d_array_from_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 				:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN)		:: dset
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: rdata
    REAL, DIMENSION(:), ALLOCATABLE 			:: raw_data
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT) 	:: attr
    CHARACTER(MAX_STRING_LENGTH) 			:: aname
    INTEGER(HID_T) 					:: dset_id
    INTEGER(HID_T) 					:: dspace_id
    INTEGER(HID_T) 					:: aspace_id
    INTEGER(HID_T) 					:: attr_id
    INTEGER(HID_T) 					:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 			:: dims
    INTEGER(HSIZE_T), DIMENSION(1) 			:: adims
    INTEGER 						:: h5error

    dims = (/ shape(rdata) /)

    ALLOCATE( raw_data(dims(1)) )

    ! * * * Read data from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
    end if
    rdata = REAL(raw_data,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
    end if

    DEALLOCATE( raw_data )

    if (PRESENT(attr)) then
       ! * * * Read data attribute(s) from file * * *
    end if

    ! * * * Read data from file * * *
  end subroutine rload_1d_array_from_hdf5

  !! @note Subroutine to load a 2-D array of reals from an HDF5 file.
  !! @details The dimensions of the 2-D array rdata is determined by the input-output array rdata.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[out] rdata 2-D array of real values read from HDF5 file and casted to KORC's real precision type.
  !! @param[out] attr 2-D array of attributes of rdata.
  !! @param raw_data 2-D array read from HDF5 file.
  !! @param aname Name of rdata attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attributes of rdata.
  subroutine rload_2d_array_from_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) 							:: rdata
    REAL, DIMENSION(:,:), ALLOCATABLE 												:: raw_data
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(MAX_STRING_LENGTH) 													:: aname
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(2) 													:: dims
    INTEGER(HSIZE_T), DIMENSION(2) 													:: adims
    INTEGER 																		:: h5error

    dims = shape(rdata)

    ALLOCATE( raw_data(dims(1),dims(2)) )

    ! * * * Read data from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
    end if
    rdata = REAL(raw_data,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
    end if

    DEALLOCATE( raw_data )

    if (PRESENT(attr)) then
       ! * * * Read data attribute(s) from file * * *
    end if

    ! * * * Read data from file * * *
  end subroutine rload_2d_array_from_hdf5

  !! @note Subroutine to load a 3-D array of reals from an HDF5 file.
  !! @details The dimensions of the 3-D array rdata is determined by the input-output array rdata.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[out] rdata 3-D array of real values read from HDF5 file and casted to KORC's real precision type.
  !! @param[out] attr 3-D array of attributes of rdata.
  !! @param raw_data 3-D array read from HDF5 file.
  !! @param aname Name of rdata attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attributes of rdata.
  subroutine rload_3d_array_from_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) 							:: rdata
    REAL, DIMENSION(:,:,:), ALLOCATABLE 											:: raw_data
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(MAX_STRING_LENGTH) 													:: aname
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(3) 													:: dims
    INTEGER(HSIZE_T), DIMENSION(3) 													:: adims
    INTEGER 																		:: h5error

    dims = shape(rdata)

    ALLOCATE( raw_data(dims(1),dims(2),dims(3)) )

    ! * * * Read data from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
    end if
    rdata = REAL(raw_data,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
    end if

    DEALLOCATE( raw_data )

    if (PRESENT(attr)) then
       ! * * * Read data attribute(s) from file * * *
    end if

    ! * * * Read data from file * * *
  end subroutine rload_3d_array_from_hdf5

  !! @note Subroutine to write a 1 byte (8 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i1save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=1), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idatum,idef), dims, h5error)

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i1save_to_hdf5

  !! @note Subroutine to write a 2 byte (16 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i2save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=2), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idatum,idef), dims, h5error)

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i2save_to_hdf5

  !! @note Subroutine to write a 4 byte (32 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i4save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=4), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idatum,idef), dims, h5error)

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i4save_to_hdf5

  !! @note Subroutine to write a 8 byte (64 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i8save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=8), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, REAL(idatum,8), dims, h5error)


    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i8save_to_hdf5

  !! @note Subroutine to write a 1-D array of integer values to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] idata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of idata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of idata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  !! @bug When using a 1-D array of attributes, only the first attribute is saved.
  subroutine isave_1d_array_to_hdf5(h5file_id,dset,idata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    INTEGER, DIMENSION(:), INTENT(IN) 												:: idata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER(SIZE_T) 																:: tmplen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(idata))
    ALLOCATE(dims(rank))
    dims = shape(idata)

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, idata, dims, h5error)

    if (PRESENT(attr)) then
       arank = size(shape(attr))
       ALLOCATE(adims(arank))
       adims = shape(attr)

       ! * * * Write attribute of data to file * * *
       tmplen = 0
       attrlen = 0
       do rr=1_idef,arank
          do dd=1_idef,adims(rr)
             tmplen = LEN_TRIM(attr(dd))
             if ( tmplen .GT. attrlen) then
                attrlen = tmplen
             end if
          end do
       end do

       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *

       DEALLOCATE(adims)
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine isave_1d_array_to_hdf5

  !! @note Subroutine to write a real to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] rdatum Real datum written to HDF5 file.
  !! @param[in] attr Attribute of datum written to HDF5 file.
  !! @param aname Name of rdatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data written to HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of rdatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of rdatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine rsave_to_hdf5(h5file_id,dset,rdatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    REAL(rp), INTENT(IN) 								:: rdatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdatum, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdatum,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine rsave_to_hdf5

  !! @note Subroutine to write a 1-D array of real values to an HDF5 file.
  !!
  !! @bug When using a 1-D array of attributes, only the first attribute is saved.
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] rdata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of rdata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of rdata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param tmplen Temporary length of rdata attribute's name.
  !! @param attrlen Lenght of rdata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  subroutine rsave_1d_array_to_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:), INTENT(IN) 												:: rdata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: tmplen
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(rdata))
    ALLOCATE(dims(rank))
    dims = shape(rdata)

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       arank = size(shape(attr))
       ALLOCATE(adims(arank))
       adims = shape(attr)

       ! * * * Write attribute of data to file * * *
       tmplen = 0
       attrlen = 0
       do rr=1_idef,arank
          do dd=1_idef,adims(rr)
             tmplen = LEN_TRIM(attr(dd))
             if ( tmplen .GT. attrlen) then
                attrlen = tmplen
             end if
          end do
       end do

       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *

       DEALLOCATE(adims)
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine rsave_1d_array_to_hdf5

  !! @note Subroutine to write a 2-D array of real values to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] rdata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of rdata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of rdata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of rdata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  !! @todo Implement the writting of attributes to HDF5 file.
  subroutine rsave_2d_array_to_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:,:), INTENT(IN) 											:: rdata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(rdata))
    ALLOCATE(dims(rank))
    dims = shape(rdata)

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine rsave_2d_array_to_hdf5

  !! @note Subroutine to write a 3-D array of real values to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] rdata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of rdata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of rdata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of rdata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  !! @todo Implement the writting of attributes to HDF5 file.
  subroutine rsave_3d_array_to_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:,:,:), INTENT(IN) 											:: rdata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(rdata))
    ALLOCATE(dims(rank))
    dims = shape(rdata)

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine rsave_3d_array_to_hdf5

  !! @note Subroutine to write an array of strings to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the array of strings.
  !! @param[in] string_array Array of characters containing the strings to be written to HDF5 file.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param dims Number of strings to be written to file.
  !! @param data_dims Dimensions of data written to HDF5 file. This is equal to (Maximum length of KORC string)x(Number of strings).
  !! @param str_len Size of strings to be written to file without blank spaces.
  !! @param string_type Native HDF5 string type.
  !! @param h5error HDF5 error status.
  subroutine save_string_parameter(h5file_id,dset,string_array)
    INTEGER(HID_T), INTENT(IN) 								:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 				:: dset
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), INTENT(IN) 	:: string_array
    INTEGER(HID_T) 											:: dset_id
    INTEGER(HID_T) 											:: dspace_id
    INTEGER(HSIZE_T), DIMENSION(1) 							:: dims
    INTEGER(HSIZE_T), DIMENSION(2) 							:: data_dims
    INTEGER(SIZE_T), DIMENSION(:), ALLOCATABLE 				:: str_len
    INTEGER(HID_T) 											:: string_type
    INTEGER 												:: h5error

    ALLOCATE(str_len(SIZE(string_array)))

    dims = (/SIZE(string_array)/)
    data_dims = (/MAX_STRING_LENGTH,SIZE(string_array)/)
    str_len = (/LEN_TRIM(string_array)/)

    call h5tcopy_f(H5T_STRING,string_type,h5error)
    call h5tset_strpad_f(string_type,H5T_STR_SPACEPAD_F,h5error)

    call h5screate_simple_f(1,dims,dspace_id,h5error)

    call h5dcreate_f(h5file_id,TRIM(dset),string_type,dspace_id,dset_id,h5error)

    call h5dwrite_vl_f(dset_id,string_type,string_array,data_dims,str_len,h5error,dspace_id)

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)

    DEALLOCATE(str_len)
  end subroutine save_string_parameter


  subroutine save_simulation_parameters(params,spp,F)
    !! @note Subroutine to save to a HDF5 file all the relevant simulation
    !! parameters. @endnote
    !! This subroutine saves to the HDF5 file "<a>simulation_parameters.h5</a>"
    !! all the relevant simulation parameters of KORC, most of them being part
    !! of the input file, but also including some derived quantities from the
    !! input parameters. This file is intended to facilitate the
    !! post-processing of KORC data using any software that supports
    !! the HDF5 software.
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !!Core KORC simulation parameters.
    TYPE(SPECIES), INTENT(IN) 	:: spp
    !! An instance of KORC's derived type SPECIES containing all
    !! the information of different electron species. See [[korc_types]].
    TYPE(FIELDS), INTENT(IN) 					:: F
    !! An instance of KORC's derived type FIELDS containing all the information
    !! about the fields used in the simulation. See [[korc_types]]
    !! and [[korc_fields]].
    CHARACTER(MAX_STRING_LENGTH) 				:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 				:: gname
    !! String containing the group name of a set of KORC parameters.
    CHARACTER(MAX_STRING_LENGTH) 				:: dset
    !! Name of data set to be saved to file.
    INTEGER(HID_T) 						:: h5file_id
    !!  HDF5 file identifier.
    INTEGER(HID_T) 						:: group_id
    !! HDF5 group identifier.
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 		:: dims
    !! Dimensions of data saved to HDF5 file.
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: rdata
    !! 1-D array of real data to be saved to HDF5 file.
    INTEGER, DIMENSION(:), ALLOCATABLE 				:: idata
    !! 1-D array of integer data to be saved to HDF5 file.
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
    !! An 1-D array with attributes of 1-D real or integer arrays that are
    !! passed to KORC interfaces of HDF5 I/O subroutines.
    CHARACTER(MAX_STRING_LENGTH) 				:: attr
    !!  A single attributes of real or integer data that is passed to KORC
    !! interfaces of HDF5 I/O subroutines.
    INTEGER 							:: h5error
    !! HDF5 error status.
    CHARACTER(19) 						:: tmp_str
    !! Temporary string used to manipulate various strings.
    REAL(rp) 							:: units
    !! Temporary variable used to add physical units to KORC parameters.

    ! * * * Error handling * * * !
    call h5eset_auto_f(params%HDF5_error_handling, h5error)
    ! Turn off: 0_idef. Turn on: 1_idef

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("Saving simulations parameters")')
    end if

    write(tmp_str,'(I18)') params%mpi_params%rank
    filename = TRIM(params%path_to_outputs) // "file_"  &
         // TRIM(ADJUSTL(tmp_str)) // ".h5"
    call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)
    call h5fclose_f(h5file_id, h5error)

    if (params%mpi_params%rank .EQ. 0) then
      if(params%output_orbit) then
         OPEN(UNIT=orbit_unit_write, &
         FILE=TRIM(params%path_to_outputs)//"orbit.out", &
         STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')
      end if
    end if

    if (params%mpi_params%rank .EQ. 0) then
       filename = TRIM(params%path_to_outputs) // "simulation_parameters.h5"

       call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

       ! Simulation parameters group
       gname = "simulation"
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

       ALLOCATE(attr_array(1))
       ALLOCATE(idata(1))

       dset = TRIM(gname) // "/field_model"
       call save_string_parameter(h5file_id,dset,(/params%field_model/))

       dset = TRIM(gname) // "/num_punctures"
       call save_to_hdf5(h5file_id,dset,params%num_punctures)

       dset = TRIM(gname) // "/dx"
       call save_to_hdf5(h5file_id,dset,params%dx)

       dset = TRIM(gname) // "/num_omp_threads"
       attr = "Number of omp threads"
       call save_to_hdf5(h5file_id,dset, params%num_omp_threads,attr)

       dset = TRIM(gname) // "/HDF5_error_handling"
       attr_array(1) = "Error handling option: 0=OFF, 1=ON"
       idata = params%HDF5_error_handling
       call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

       dset = TRIM(gname) // "/nmpi"
       attr = "Number of mpi processes"
       call save_to_hdf5(h5file_id,dset,params%mpi_params%nmpi,attr)

       dset = TRIM(gname) // "/phi_section"
       attr = "Phi of puncture plane"
       call save_to_hdf5(h5file_id,dset,params%phi_section,attr)

       DEALLOCATE(idata)
       DEALLOCATE(attr_array)

       call h5gclose_f(group_id, h5error)


       ! Plasma species group
       gname = "species"
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

       dset = TRIM(gname) // "/ppp"
       attr = "Particles per (mpi) process"
       call save_to_hdf5(h5file_id,dset,spp%ppp,attr)

       dset = TRIM(gname) // "/spatial_distribution"
       call save_string_parameter(h5file_id,dset,(/spp%spatial_distrib/))

       dset = TRIM(gname) // "/Y"
       call rsave_2d_array_to_hdf5(h5file_id, dset, spp%vars%Y)

       call h5gclose_f(group_id, h5error)


       ! Electromagnetic fields group

       gname = "fields"
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

       if (TRIM(params%field_model(1:10)) .EQ. 'ANALYTICAL') then
          dset = TRIM(gname) // "/Bo"
          attr = "Toroidal field at the magnetic axis in T"
          call save_to_hdf5(h5file_id,dset,F%Bo,attr)

          dset = TRIM(gname) // "/current_direction"
          call save_string_parameter(h5file_id,dset, &
               (/F%AB%current_direction/))

          dset = TRIM(gname) // "/a"
          attr = "Minor radius in m"
          call save_to_hdf5(h5file_id,dset,F%AB%a,attr)

          dset = TRIM(gname) // "/Ro"
          attr = "Magnetic axis radial position"
          call save_to_hdf5(h5file_id,dset,F%Ro,attr)

          dset = TRIM(gname) // "/Zo"
          attr = "Magnetic axis vertical position"
          call save_to_hdf5(h5file_id,dset,F%Zo,attr)

          dset = TRIM(gname) // "/qa"
          attr = "Safety factor at minor radius"
          call save_to_hdf5(h5file_id,dset,F%AB%qa,attr)

          dset = TRIM(gname) // "/qo"
          attr = "Safety factor at the magnetic axis"
          call save_to_hdf5(h5file_id,dset,F%AB%qo,attr)

          dset = TRIM(gname) // "/lambda"
          attr = "Parameter lamda in m"
          call save_to_hdf5(h5file_id,dset,F%AB%lambda,attr)

          dset = TRIM(gname) // "/Bpo"
          attr = "Poloidal magnetic field in T"
          call save_to_hdf5(h5file_id,dset,F%AB%Bpo,attr)

       else if ((params%field_model .EQ. 'PSPLINE').or. &
            (params%field_model .EQ. 'MARS')) then
          ALLOCATE(attr_array(1))

          dset = TRIM(gname) // "/dims"
          attr_array(1) = "Mesh dimension of the magnetic  &
               field (NR,NPHI,NZ)"
          call save_1d_array_to_hdf5(h5file_id,dset,F%dims,attr_array)

          dset = TRIM(gname) // "/R"
          attr_array(1) = "Radial position of the magnetic field grid nodes"
          call save_1d_array_to_hdf5(h5file_id,dset,F%X%R,attr_array)

          if (ALLOCATED(F%X%PHI)) then
             dset = TRIM(gname) // "/PHI"
             attr_array(1) = "Azimuthal angle of the magnetic &
                  field grid nodes"
             call save_1d_array_to_hdf5(h5file_id,dset,F%X%PHI,attr_array)
          end if

          dset = TRIM(gname) // "/Z"
          attr_array(1) = "Z position of the magnetic field grid nodes"
          call save_1d_array_to_hdf5(h5file_id,dset,F%X%Z,attr_array)

          dset = TRIM(gname) // "/Bo"
          attr = "Toroidal field at the magnetic axis in T"
          call save_to_hdf5(h5file_id,dset,F%Bo,attr)

          dset = TRIM(gname) // "/Ro"
          attr = "Radial position of magnetic axis"
          call save_to_hdf5(h5file_id,dset,F%Ro,attr)

          dset = TRIM(gname) // "/Zo"
          attr = "Radial position of magnetic axis"
          call save_to_hdf5(h5file_id,dset,F%Zo,attr)

          dset = TRIM(gname) // "/Axisymmetric"
          attr = "Axisymmetry"
          if(F%axisymmetric_fields) then
             call save_to_hdf5(h5file_id,dset,1_idef,attr)
          else
             call save_to_hdf5(h5file_id,dset,0_idef,attr)
          end if

          dset = TRIM(gname) // "/Dim2x1t"
          attr = "Dim2x1t"
          if(F%Dim2x1t) then
             call save_to_hdf5(h5file_id,dset,1_idef,attr)
          else
             call save_to_hdf5(h5file_id,dset,0_idef,attr)
          end if

          dset = TRIM(gname) // "/ind_2x1t"
          attr = "Index for 2x1t"
          call save_to_hdf5(h5file_id,dset,F%ind_2x1t,attr)

          dset = TRIM(gname) // "/psip_conv"
          attr = "Scaling factor for magnetic poloidal flux function"
          call save_to_hdf5(h5file_id,dset,F%psip_conv,attr)

          if (ALLOCATED(F%PSIp)) then
             dset = TRIM(gname) // "/psi_p"
             call rsave_2d_array_to_hdf5(h5file_id, dset,F%PSIp)
          end if

          if (ALLOCATED(F%PSIp3D)) then
             dset = TRIM(gname) // "/psi_p3D"
             call rsave_3d_array_to_hdf5(h5file_id, dset,F%PSIp3D)
          end if

          if (ALLOCATED(F%FLAG2D)) then
             dset = TRIM(gname) // "/Flag2D"
             call rsave_2d_array_to_hdf5(h5file_id, dset, &
                  F%FLAG2D)
          end if

          if (ALLOCATED(F%FLAG3D)) then
             dset = TRIM(gname) // "/Flag3D"
             call rsave_3d_array_to_hdf5(h5file_id, dset, &
                  F%FLAG3D)
          end if


          if  (F%Bfield) then

             if (F%axisymmetric_fields) then

                if (ALLOCATED(F%B_2D%R)) then
                   dset = TRIM(gname) // "/BR"
                   call rsave_2d_array_to_hdf5(h5file_id, dset,F%B_2D%R)
                end if

                if (ALLOCATED(F%B_2D%PHI)) then
                   dset = TRIM(gname) // "/BPHI"
                   call rsave_2d_array_to_hdf5(h5file_id, dset,F%B_2D%PHI)
                end if

                if (ALLOCATED(F%B_2D%Z)) then
                   dset = TRIM(gname) // "/BZ"
                   call rsave_2d_array_to_hdf5(h5file_id, dset,F%B_2D%Z)
                end if


             else

                dset = TRIM(gname) // "/BR"
                call rsave_3d_array_to_hdf5(h5file_id, dset,F%B_3D%R)

                dset = TRIM(gname) // "/BPHI"
                call rsave_3d_array_to_hdf5(h5file_id, dset,F%B_3D%PHI)

                dset = TRIM(gname) // "/BZ"
                call rsave_3d_array_to_hdf5(h5file_id, dset,F%B_3D%Z)

             end if

          endif

          if  (F%B1field) then

             dset = TRIM(gname) // "/AMP"
             attr = "Amplitude of perturbation field"
             call save_to_hdf5(h5file_id,dset,F%AMP,attr)

             if (ALLOCATED(F%B1Re_2D%R)) then
                dset = TRIM(gname) // "/B1Re_R"
                call rsave_2d_array_to_hdf5(h5file_id, dset,F%B1Re_2D%R)
             end if

             if (ALLOCATED(F%B1Re_2D%PHI)) then
                dset = TRIM(gname) // "/B1Re_PHI"
                call rsave_2d_array_to_hdf5(h5file_id, dset,F%B1Re_2D%PHI)
             end if

             if (ALLOCATED(F%B1Re_2D%Z)) then
                dset = TRIM(gname) // "/B1Re_Z"
                call rsave_2d_array_to_hdf5(h5file_id, dset,F%B1Re_2D%Z)
             end if

             if (ALLOCATED(F%B1Im_2D%R)) then
                dset = TRIM(gname) // "/B1Im_R"
                call rsave_2d_array_to_hdf5(h5file_id, dset,F%B1Im_2D%R)
             end if

             if (ALLOCATED(F%B1Im_2D%PHI)) then
                dset = TRIM(gname) // "/B1Im_PHI"
                call rsave_2d_array_to_hdf5(h5file_id, dset,F%B1Im_2D%PHI)
             end if

             if (ALLOCATED(F%B1Im_2D%Z)) then
                dset = TRIM(gname) // "/B1Im_Z"
                call rsave_2d_array_to_hdf5(h5file_id, dset,F%B1Im_2D%Z)
             end if

          endif

          DEALLOCATE(attr_array)

       end if

       call h5gclose_f(group_id, h5error)

       call h5fclose_f(h5file_id, h5error)
    end if

  end subroutine save_simulation_parameters


  subroutine save_simulation_outputs(params,spp,F)
    !! @note Subroutine that saves the electrons' variables specified in
    !! params::outputs_list to HDF5 files. @endnote
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !! Core KORC simulation parameters.
    TYPE(SPECIES), INTENT(IN) 	:: spp
    !! An instance of KORC's derived type SPECIES containing all
    !! the information
    !! of different electron species. See [[korc_types]].
    TYPE(FIELDS), INTENT(IN)                 :: F
    CHARACTER(MAX_STRING_LENGTH) 				:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 				:: gname
    !! String containing the group name of a set of KORC parameters.
    CHARACTER(MAX_STRING_LENGTH) 				:: subgname
    !! String containing the subgroup name of a set of KORC parameters.
    CHARACTER(MAX_STRING_LENGTH) 				:: dset
    !! Name of data set to be saved to file.
    INTEGER(HID_T) 						:: h5file_id
    !! HDF5 file identifier.
    INTEGER(HID_T) 						:: group_id
    !! HDF5 group identifier.
    INTEGER(HID_T) 						:: subgroup_id
    !! HDF5 subgroup identifier.
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 		:: dims
    !! Dimensions of data saved to HDF5 file.
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: rdata
    !! 1-D array of real data to be saved to HDF5 file.
    INTEGER, DIMENSION(:), ALLOCATABLE 				:: idata
    !!1-D array of integer data to be saved to HDF5 file.
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE       :: attr_array
    !! An 1-D array with attributes of 1-D real or integer arrays that are
    !! passed to KORC interfaces of HDF5 I/O subroutines.
    CHARACTER(MAX_STRING_LENGTH) 				:: attr
    !! A single attributes of real or integer data that is passed to KORC
    !! interfaces of HDF5 I/O subroutines.
    INTEGER 							:: h5error
    !!HDF5 error status.
    CHARACTER(19) 						:: tmp_str
    !!Temporary string used to manipulate various strings.
    REAL(rp) 						:: units
    !! Temporary variable used to add physical units to electrons' variables.
    INTEGER 						:: ss
    !! Electron species iterator.
    INTEGER 						:: jj
    !! Iterator for reading all the entried of params::outputs_list.
    LOGICAL 						:: object_exists
    !! Flag determining if a certain dataset is already present in
    !! the HDF5 output files.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE  ::YY
    !! Temporary variable get proper units on vars%Y(1,:) and vars%Y(3,:), which
    !! are lengths, while keeping vars%Y(2,:), which is an angle


    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("Saving simulations outputs")')
    end if

    write(tmp_str,'(I18)') params%mpi_params%rank
    filename = TRIM(params%path_to_outputs) // "file_" &
         // TRIM(ADJUSTL(tmp_str)) // ".h5"
    call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    gname = 'out'
    call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

    if (.NOT.object_exists) then ! Check if group does exist.
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)


       dset = "Punctures"
       call rsave_3d_array_to_hdf5(group_id, dset, &
            spp%vars%punct)

       dset = "ConnectionLength"
       call rsave_1d_array_to_hdf5(group_id, dset, &
            spp%vars%con_len)

       dset = "FlagCon"
       call save_1d_array_to_hdf5(group_id, dset, &
            INT(spp%vars%flagCon,idef))


       call h5gclose_f(group_id, h5error)
    end if ! Check if group does exist.

    call h5fclose_f(h5file_id, h5error)

  end subroutine save_simulation_outputs



end module PB_HDF5
