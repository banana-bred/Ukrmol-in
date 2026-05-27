# ------------------------------------------------------------------------------------------------
# libscatci library

set(libscatci-sources
    "davidson.f"
    "m_2_matrix_sparse.f90"
    "scatci_data.f"
    "scatci_diag.F"
    "scatci_routines.f"
)

if(ARPACK_LIBRARIES)
    message(STATUS "SCATCI will be built with Arpack")
    add_definitions(-Darpack)
    set(libscatci-sources
        ${libscatci-sources}
        "m_3_arpack_diagonalizers.f90"
        "m_4_mkarp.f90")
endif()

add_library(libscatci ${libscatci-sources})

target_compile_options(libscatci PUBLIC ${OpenMP_Fortran_FLAGS})

target_link_libraries(libscatci
    libukplus_interfaces
    libukplus_utilities
    libGBTO
    libukplus_global
    ${ARPACK_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

target_include_directories(libscatci PUBLIC
    "${CMAKE_BINARY_DIR}/source/gbtolib/mod"
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/interfaces/mod"
)

set_target_properties(libscatci PROPERTIES
    OUTPUT_NAME scatci
    LINK_FLAGS ${OpenMP_Fortran_FLAGS}
)

# ------------------------------------------------------------------------------------------------
# scatci program

add_executable(scatci "scatci_build.f")

target_include_directories(scatci PUBLIC
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_link_libraries(scatci
    libscatci
)

set_target_properties(scatci PROPERTIES
    LINK_FLAGS ${OpenMP_Fortran_FLAGS}
)

# ------------------------------------------------------------------------------------------------
# installation targets

install(TARGETS libscatci scatci
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ci-diag")
