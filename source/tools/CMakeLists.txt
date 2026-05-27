# --------------------------------------------------------------------------------
# RADDEN

add_executable(radden
    "radden/source/dqags.f"
    "radden/source/radden_modules.f90"
    "radden/source/radden.f90"
)

target_link_libraries(radden
    libukplus_utilities
    libukplus_global
)

target_include_directories(radden PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

# --------------------------------------------------------------------------------
# PHASE_MATCH (two versions)

add_executable(phase_match          "phase_match/phase_match.f90")
add_executable(phase_match_orbitals "phase_match/phase_match_orbitals.f90")

target_compile_options(phase_match          PUBLIC ${OpenMP_Fortran_FLAGS})
target_compile_options(phase_match_orbitals PUBLIC ${OpenMP_Fortran_FLAGS})

set_target_properties(phase_match          PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})
set_target_properties(phase_match_orbitals PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})

target_link_libraries(phase_match
    libcdenprop
    libscatci
    libGBTO
    libukplus_global
    libukplus_interfaces
    libukplus_utilities
    ${SCALAPACK_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

target_link_libraries(phase_match_orbitals
    libGBTO
    libukplus_global
    ${BLAS_LIBRARIES}
)

target_include_directories(phase_match PUBLIC
    "${CMAKE_BINARY_DIR}/source/ci-diag/mod"
    "${CMAKE_BINARY_DIR}/source/cdenprop/mod"
    "${CMAKE_BINARY_DIR}/source/gbtolib/mod"
    "${CMAKE_BINARY_DIR}/source/global/mod"
)

target_include_directories(phase_match_orbitals PUBLIC
    "${CMAKE_BINARY_DIR}/source/gbtolib/mod"
    "${CMAKE_BINARY_DIR}/source/interfaces/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
    "${CMAKE_BINARY_DIR}/source/global/mod"
)

# --------------------------------------------------------------------------------
# installation targets

install(TARGETS phase_match phase_match_orbitals radden
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
