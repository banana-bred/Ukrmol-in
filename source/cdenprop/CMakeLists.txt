# ------------------------------------------------------------------------------------------------
# cdenprop library

add_library(libcdenprop
    "class_COOSparseMatrix_integer.f90"
    "class_COOSparseMatrix_real.f90"
    "class_namelists.f90"
    "class_molecular_properties_data.f90"
    "cdenprop_aux.f90"
    "cdenprop_defs.F90"
    "cdenprop_io.F90"
    "cdenprop_procs.F90"
)

target_include_directories(libcdenprop PUBLIC
    "${CMAKE_BINARY_DIR}/source/gbtolib/mod"
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/interfaces/mod"
)

set_target_properties(libcdenprop PROPERTIES
    OUTPUT_NAME cdenprop
    LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
)

target_compile_options(libcdenprop PUBLIC ${OpenMP_Fortran_FLAGS} ${GBTOlib_Fortran_FLAGS})

if(WITH_SCALAPACK OR SCALAPACK_LIBRARIES OR CMAKE_Fortran_FLAGS MATCHES "-q?mkl(-ilp64)?=cluster")
    message(STATUS "CDENPROP will be built with ScaLAPACK")
    add_definitions(-Dscalapack)
endif()

target_link_libraries(libcdenprop
    libukplus_interfaces
    libGBTO
    libukplus_global
    ${SCALAPACK_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

# ------------------------------------------------------------------------------------------------
# cdenprop

add_executable(cdenprop "cdenprop.f90")

target_include_directories(cdenprop PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

set_target_properties(cdenprop PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})

target_link_libraries(cdenprop
    libcdenprop
    libukplus_utilities
)

# ------------------------------------------------------------------------------------------------
# cdenprop_all

add_executable(cdenprop_all "cdenprop_all.f90")

target_include_directories(cdenprop_all PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_link_libraries(cdenprop_all
    libcdenprop
    libukplus_utilities
    libGBTO
)

set_target_properties(cdenprop_all PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})

# ------------------------------------------------------------------------------------------------
# installation targets

install(TARGETS libcdenprop cdenprop cdenprop_all
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cdenprop")
