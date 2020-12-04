#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.
set (poinbabel_sources "")
list (APPEND poinbabel_sources
    interp.f90
    constants.f90
    field.f90
    finalize.f90
    ppusher.f90
    HDF5.f90
    types.f90
    hpc.f90
    initialize.f90
    fio_interface.f90
    input.f90
    coords.f90
    spatial_distribution.f90
    main.f90
)

add_executable (xpoinbabel ${poinbabel_sources})
set_property(TARGET xpoinbabel PROPERTY LINKER_LANGUAGE Fortran)


target_link_libraries (xpoinbabel ${MPI_Fortran_LIBRARIES})
target_link_libraries (xpoinbabel ${HDF5_LIBRARIES})

if(${USE_PSPLINE})
   target_link_libraries (xpoinbabel ${EZSpline_LIBRARIES})
endif()

if(${USE_FIO})
   target_link_libraries (xpoinbabel ${FIO_LIBRARIES} ${M3DC1_LIBRARIES})
   target_link_libraries (xpoinbabel ${NIMFIO_LIBRARIES} ${NIMBLK_LIBRARIES} ${NIMBTY_LIBRARIES} ${NIMLOC_LIBRARIES} ${NIMITER_LIBRARIES} ${NIMMAT_LIBRARIES} ${NIMPAR_LIBRARIES} ${NIMLIB_LIBRARIES} ${NIMMPI_LIBRARIES} ${NIMLAP_LIBRARIES} ${NIMSLU_LIBRARIES} ${NIMLUD_LIBRARIES})
endif()
