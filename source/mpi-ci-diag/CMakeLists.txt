# -------------------------------------------------------------------------------------------------------#
# UKRmol+ suite - MPI-SCATCI                                                                             #
# -------------------------------------------------------------------------------------------------------#

# -------------------------------------------------------------------------------------------------------#
# MPI-SCATCI source files                                                                                #
# -------------------------------------------------------------------------------------------------------#

set(mpi-scatci-sources
    "Modules/Diagonalizers/DavidsonDiagonalizer_module.f90"
    "Modules/Diagonalizers/Diagonalizer_module.f90"
    "Modules/Diagonalizers/DiagonalizerResult_module.f90"
    "Modules/Diagonalizers/LapackDiagonalizer_module.f90"
    "Modules/Dispatcher/Dispatcher_module.F90"
    "Modules/Globals/consts_mpi_ci.f90"
    "Modules/Hamiltonian/CI_Hamiltonian_module.f90"
    "Modules/Hamiltonian/Contracted_Hamiltonian_module.f90"
    "Modules/Hamiltonian/Contracted_Symbolic_Module.f90"
    "Modules/Hamiltonian/CSF_module.F90"
    "Modules/Hamiltonian/Hamiltonian_module.f90"
    "Modules/Hamiltonian/Orbital_module.f90"
    "Modules/Hamiltonian/SolutionHandler_module.F90"
    "Modules/Hamiltonian/Symbolic_Module.f90"
    "Modules/Hamiltonian/Target_RMat_CI_module.f90"
    "Modules/Hamiltonian/Uncontracted_Hamiltonian_module.f90"
    "Modules/Integral/ALCHEMY_Module.f90"
    "Modules/Integral/BaseIntegral_module.f90"
    "Modules/Integral/SWEDEN_Module.F90"
    "Modules/Integral/UKRMOL_module.f90"
    "Modules/Matrix/BaseMatrix_module.f90"
    "Modules/Matrix/DistributedMatrix_module.f90"
    "Modules/Matrix/MatrixCache_module.f90"
    "Modules/Matrix/MatrixElement_module.f90"
    "Modules/Matrix/WriterMatrix_module.f90"
    "Modules/MPI-Objects/Parallelization_module.F90"
    "Modules/Options/Options_module.f90"
    "Modules/Utilities/MemoryManager_module.f90"
    "Modules/Utilities/Postprocessing_module.F90"
    "Modules/Utilities/Timing_Module.f90"
    "Modules/Utilities/Utility_module.f90"
    "mpi_scatci.f90"
)

# -------------------------------------------------------------------------------------------------------#
# Use optional components                                                                                #
# -------------------------------------------------------------------------------------------------------#

if(ARPACK_LIBRARIES)
    message(STATUS "MPI-SCATCI will be built with Arpack")
    add_definitions(-Darpack)
    list(APPEND mpi-scatci-sources
        "Modules/Diagonalizers/ARPACKDiagonalizer_module.f90")
endif()

if(WITH_SCALAPACK OR SCALAPACK_LIBRARIES OR CMAKE_Fortran_FLAGS MATCHES "-q?mkl(-ilp64)?=cluster")
    message(STATUS "MPI-SCATCI will be built with ScaLAPACK")
    add_definitions(-Dscalapack)
    list(APPEND mpi-scatci-sources
        "Modules/Diagonalizers/SCALAPACKDiagonalizer_module.f90"
        "Modules/Matrix/SCALAPACKMatrix_module.f90")
    if(ELPA_LIBRARIES)
        message(STATUS "MPI-SCATCI will be built with ELPA")
        add_definitions(-Duseelpa)
        list(APPEND mpi-scatci-sources
            "Modules/Diagonalizers/ELPADiagonalizer_module.f90"
            "Modules/Matrix/ELPAMatrix_module.f90")
    endif()
    if(SLEPC_LIBRARIES)
        message(STATUS "MPI-SCATCI will be built with SLEPc")
        add_definitions(-Dslepc)
        list(APPEND mpi-scatci-sources
            "Modules/Diagonalizers/SLEPCDiagonalizer_module.F90"
            "Modules/Diagonalizers/SLEPCDavidsonDiagonalizer_module.F90"
            "Modules/Matrix/SLEPCMatrix_module.F90")
    endif()
endif()

# -------------------------------------------------------------------------------------------------------#
# Define the MPI-SCATCI executable target                                                                #
# -------------------------------------------------------------------------------------------------------#

add_executable(mpi-scatci ${mpi-scatci-sources})

target_compile_options(mpi-scatci PUBLIC ${OpenMP_Fortran_FLAGS} ${GBTOlib_Fortran_FLAGS})

target_include_directories(mpi-scatci PUBLIC
    "${SLEPC_INCLUDE_DIRS}"
    "${ELPA_INCLUDE_DIRS}"
    "${CMAKE_BINARY_DIR}/source/gbtolib/mod"
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/ci-diag/mod"
    "${CMAKE_BINARY_DIR}/source/cdenprop/mod"
    "${CMAKE_BINARY_DIR}/source/interfaces/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

set_target_properties(mpi-scatci
    PROPERTIES
        LINK_FLAGS ${OpenMP_Fortran_FLAGS}
        LINKER_LANGUAGE Fortran
)

target_link_libraries(mpi-scatci
    libcdenprop
    libscatci
    libukplus_interfaces
    libukplus_utilities
    libGBTO
    libukplus_global
    ${ELPA_LIBRARIES}
    ${SLEPC_LIBRARIES}
    ${SCALAPACK_LIBRARIES}
    ${ARPACK_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

# -------------------------------------------------------------------------------------------------------#
# Installation targets                                                                                   #
# -------------------------------------------------------------------------------------------------------#

install(TARGETS mpi-scatci
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ci-diag")

# -------------------------------------------------------------------------------------------------------#
# Build Doxygen documentation                                                                            #
# -------------------------------------------------------------------------------------------------------#

if(DOXYGEN_FOUND)
    set(DOXYGEN_IN "${CMAKE_CURRENT_LIST_DIR}/doc/Doxyfile")

    add_custom_target(mpi-scatci-doc ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/doc"
        COMMENT "Generating MPI-SCATCI documentation with Doxygen"
        VERBATIM)

   install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/doc/html/"
           DESTINATION "${CMAKE_INSTALL_DOCDIR}/mpi-scatci-doxygen")
endif()
