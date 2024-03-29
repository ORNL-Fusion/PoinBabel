cmake_minimum_required (VERSION 3.13)

cmake_policy(VERSION 3.13)

if(CMP0074)
   cmake_policy(SET CMP0074 NEW)
endif()

project (poinbabel_proj C CXX Fortran)

#  Set directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)

#  Fortran specific settings.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/modules)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/lib)	

#  Set definitions.
add_definitions (-DDOUBLE_PRECISION)

#  Configure OpenMP
find_package (OpenMP)
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")

#  Configure MPI
find_package (MPI REQUIRED)
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
include_directories (${MPI_Fortran_INCLUDE_PATH})

#  Configure HDF5
find_package (HDF5 REQUIRED COMPONENTS Fortran)
include_directories (${HDF5_INCLUDE_DIRS})

#  Configure EZspline
option (USE_PSPLINE "Use PSPLINE interpolation support" OFF)
if (${USE_PSPLINE})
    include_directories (${EZSpline_INCLUDE_PATH})
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EZSpline_LINK_FLAGS}")
    add_definitions (-DPSPLINE)
endif ()

#  Configure FIO interpolation
option (USE_FIO "Use FIO interpolation support" OFF)
if (${USE_FIO})
    include_directories (${FIO_INCLUDE_PATH} ${NIMROD_INCLUDE_PATH} ${NIMFIO_INCLUDE_PATH} ${NIMCORE_INCLUDE_PATH} ${NIMLIB_INCLUDE_PATH} ${NIMEXT_INCLUDE_PATH})
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${FIO_LINK_FLAGS} ${NIMFIO_LINK_FLAGS} ${NIMCORE_LINK_FLAGS} ${NIMLIB_LINK_FLAGS} ${NIMEXT_LINK_FLAGS}")
    add_definitions (-DFIO)
endif ()

add_subdirectory (src)

#  Add machine dependent definition for finding included scripts
option (CORI_DIR "Use Cori directories" OFF)
if (${CORI_DIR})
    add_definitions (-DCORI)
else ()
    add_definitions (-DMAC)	     
endif ()

