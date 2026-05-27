# -------------------------------------------------------------------------------------------------------#
# UKRmol+ suite - GBTO library                                                                           #
# -------------------------------------------------------------------------------------------------------#
#                                                                                                        #
# Environment variables                                                                                  #
# ---------------------                                                                                  #
#                                                                                                        #
# - BLA_VENDOR              Defines the BLAS/LAPACK vendor, typically "Intel10_64ilp" for ILP64 MKL.     #
#                           In case of Intel Compiler this needs to be supplemented with the option      #
#                           "-qmkl". (Only available with CMake 3.13+)                                   #
#                                                                                                        #
# - BLA_STATIC              Look for static BLAS/LAPACK libraries. In case of Intel Compiler, this needs #
#                           to be consistent with usage or absence of the flag "-static-intel".          #
#                           (Only available with CMake 3.13+)                                            #
#                                                                                                        #
# CMake options                                                                                          #
# -------------                                                                                          #
#                                                                                                        #
# - WITH_MPI                Build with MPI support (default: ON).                                        #
#                                                                                                        #
# - GBTOlib_Fortran_FLAGS   Force specific MPI-related compiler flags (default: determined by CMake).    #
#                                                                                                        #
# - BUILD_DOC               Generate development documentation (default: OFF).                           #
#                                                                                                        #
# - BUILD_TESTING           Enable the calculation test suite (default: ON).                             #
#                                                                                                        #
# - BUILD_UNIT_TESTS        Enable development unit tests (default: OFF).                                #
#                                                                                                        #
# - BUILD_SHARED_LIBS       Switch to building shared libraries instead of static ones (default: OFF).   #
#                                                                                                        #
# - MPIEXEC_PREFLAGS        Additional flags to pass to MPI launcher when running the test suite.        #
#                                                                                                        #
# - WITH_MOLPRO             Use Molpro with test suite (default: OFF).                                   #
#                                                                                                        #
# - WITH_PSI4               Use Psi4 with test suite (default: OFF; Molpro has priority if both set).    #
#                                                                                                        #
# -------------------------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.5...3.31)

project(GBTOlib C Fortran)

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(CMAKE_Fortran_MODULE_DIRECTORY "mod")

# -------------------------------------------------------------------------------------------------------#
# User switches; use e.g. cmake -D WITH_MPI=OFF to change                                                #
# -------------------------------------------------------------------------------------------------------#

option(WITH_MPI       "Use MPI if available."     ON)
option(BUILD_DOC      "Build documentation"       OFF)
option(BUILD_TESTING  "Enable the test suite"     ON)
option(BUILD_UNIT_TEST "Enable the development unit tests" OFF)

option(WITH_PSI4      "Use the Psi4 quantum chemistry package for tests." OFF)
option(WITH_MOLPRO    "Use Molpro quantum chemistry package for tests."   OFF)

set(GBTOlib_Fortran_FLAGS "${GBTOlib_Fortran_FLAGS}"
    CACHE STRING "MPI-related compiler parameters for GBTOlib")

# -------------------------------------------------------------------------------------------------------#
# Find the required libraries                                                                            #
# -------------------------------------------------------------------------------------------------------#

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# -------------------------------------------------------------------------------------------------------#
# Check if MPI is available, and what are its characteristics                                            #
# -------------------------------------------------------------------------------------------------------#

if(WITH_MPI)
    set(MPI_DETERMINE_LIBRARY_VERSION ON)

    find_package(MPI)

    if(MPI_Fortran_FOUND AND "${GBTOlib_Fortran_FLAGS}" STREQUAL "")
        message(STATUS "GBTOlib will be compiled with ${MPI_Fortran_LIBRARY_VERSION_STRING}")
        set(GBTOlib_Fortran_FLAGS -Dusempi -Dsplitreduce)

        if(MPI_Fortran_VERSION_MAJOR)
            if(MPI_Fortran_VERSION_MAJOR GREATER_EQUAL 3)
                message(STATUS "GBTOlib will use MPI-3 features")
                set(GBTOlib_Fortran_FLAGS ${GBTOlib_Fortran_FLAGS} -Dmpithree)
            endif()
        endif()

        string(CONCAT testMPIIntSize_src
            "program testMPIIntSize\n"
            "    use iso_fortran_env, only: int32, int64\n"
            "    use mpi\n"
            "    select case (kind(MPI_COMM_WORLD))\n"
            "        case (int32); write (*, '(\"4\")')\n"
            "        case (int64); write (*, '(\"8\")')\n"
            "        case default; write (*, '(\"Unknown\")')\n"
            "    end select\n"
            "end program\n")

        file(WRITE "${CMAKE_BINARY_DIR}/testMPIIntSize.f90" "${testMPIIntSize_src}")

        set(CMAKE_CROSSCOMPILING_EMULATOR "${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS}")
        try_run(_RUN_RESULT _COMPILE_RESULT
                "${CMAKE_BINARY_DIR}" SOURCES "${CMAKE_BINARY_DIR}/testMPIIntSize.f90"
                LINK_LIBRARIES MPI::MPI_Fortran
                COMPILE_OUTPUT_VARIABLE _MPI_INTEGER_SIZE_TEST_OUTPUT
                RUN_OUTPUT_VARIABLE _MPI_INTEGER_SIZE)
        unset(CMAKE_CROSSCOMPILING_EMULATOR)

        if(_MPI_INTEGER_SIZE MATCHES "^4")
            message(STATUS "GBTOlib will use 32-bit MPI integers")
        elseif(_MPI_INTEGER_SIZE MATCHES "^8")
            message(STATUS "GBTOlib will use 64-bit MPI integers")
        else()
            message(FATAL_ERROR "Failed to detemine MPI integer interface.\n"
                                "Maybe explicit cmake flag -DMPIEXEC_EXECUTABLE= is needed?\n"
                                "Anyway, this is the output of the compiler...\n"
                                "${_MPI_INTEGER_SIZE_TEST_OUTPUT}\n"
                                "... and this is the output from the failed test execution...\n"
                                "${_MPI_INTEGER_SIZE}")
        endif()

        string(CONCAT testMPIQuadReduce_src
            "program testMPIQuadReduce\n"
            "    use mpi\n"
            "    integer(kind=kind(MPI_COMM_WORLD)) :: ierr, pid, cnt\n"
            "    real(kind=selected_real_kind(19)), allocatable :: src(:), dst(:)\n"
            "    call MPI_Init(ierr)\n"
            "    call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)\n"
            "    call MPI_Comm_size(MPI_COMM_WORLD, cnt, ierr)\n"
            "    allocate(src(cnt), dst(cnt))\n"
            "    src(:) = 0; src(pid + 1) = 1; dst(:) = 0\n"
            "    call MPI_Allreduce(src, dst, cnt, MPI_REAL16, MPI_SUM, MPI_COMM_WORLD, ierr)\n"
            "    if (pid == 0) write(*,'(A)') merge('Y', 'N', nint(sum(dst)) == cnt)\n"
            "    call MPI_Finalize(ierr)\n"
            "end program\n")

        file(WRITE "${CMAKE_BINARY_DIR}/testQuadReduce.f90" "${testMPIQuadReduce_src}")

        set(CMAKE_CROSSCOMPILING_EMULATOR "${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};${MPIEXEC_MAX_NUMPROCS};${MPIEXEC_PREFLAGS}")
        message(STATUS "Testing MPI using: ${CMAKE_CROSSCOMPILING_EMULATOR}")
        try_run(_RUN_RESULT _COMPILE_RESULT
                "${CMAKE_BINARY_DIR}" SOURCES "${CMAKE_BINARY_DIR}/testQuadReduce.f90"
                LINK_LIBRARIES MPI::MPI_Fortran
                RUN_OUTPUT_VARIABLE _QUAD_REDUCE_WORKS)
        string(STRIP "${_QUAD_REDUCE_WORKS}" _QUAD_REDUCE_WORKS)
        unset(CMAKE_CROSSCOMPILING_EMULATOR)

        if(_QUAD_REDUCE_WORKS MATCHES "^Y")
            message(STATUS "GBTOlib will use native quad MPI reductions")
            set(GBTOlib_Fortran_FLAGS ${GBTOlib_Fortran_FLAGS} -Dquadreduceworks)
        else()
            message(STATUS "GBTOlib will use custom quad MPI reductions")
        endif()
    endif()
endif()

if(MPI_Fortran_FOUND)
    message(STATUS "GBTOlib_Fortran_FLAGS: ${GBTOlib_Fortran_FLAGS}")
else()
    message(STATUS "GBTOlib will be compiled without MPI")
endif()

# -------------------------------------------------------------------------------------------------------#
# Determine BLAS & LAPACK integer interface                                                              #
# -------------------------------------------------------------------------------------------------------#

if(NOT "${GBTOlib_Fortran_FLAGS}" MATCHES "Dblas..bit")

    # This test uses the LAPACK subroutine SLASWP to twice interchange rows in a 3-by-3 matrix to mimic
    # overall row permutation [3, 1, 2]. To achieve that, one needs to pass two-component array to the
    # subroutine. However, actual position of the byte where the second integer is supposed to start
    # depends on the integer interface of the library and one or the other case will necessarily fail (or
    # not give the same result as the straightforward in-code permutation). If both tests fail for some
    # reason, the script will set neither -Dblas64bit, nor -Dblas32bit and use blas integer equal to the
    # default integer.

    string(CONCAT testLAPACKIntSize_src
        "program mkl_interface\n"
        "    use iso_fortran_env, only: INT, real32\n"
        "    real(real32) :: A(3,3) = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), (/ 3, 3 /))\n"
        "    real(real32) :: M(3,3)\n"
        "    integer(INT) :: P(3) = (/ 3, 1, 2 /), Q(2) = (/ 3, 3 /)\n"
        "    integer(INT) :: one = 1, two = 2, three = 3\n"
        "    M(:,:) = A(:,:)\n"
        "    call slaswp(three, M, three, one, two, Q, one)\n"
        "    write (*, '(9I0)') nint(M(:,:) - A(P,:))\n"
        "end program mkl_interface\n")

    file(WRITE "${CMAKE_BINARY_DIR}/testLAPACKIntSize.F90" "${testLAPACKIntSize_src}")

    try_run(_RUN_RESULT _COMPILE_RESULT
        "${CMAKE_BINARY_DIR}" SOURCES "${CMAKE_BINARY_DIR}/testLAPACKIntSize.F90"
        LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
        COMPILE_DEFINITIONS -DINT=int32
        COMPILE_OUTPUT_VARIABLE _MESSAGES4
        RUN_OUTPUT_VARIABLE _LAPACK_INTEGER4)
    try_run(_RUN_RESULT _COMPILE_RESULT
        "${CMAKE_BINARY_DIR}" SOURCES "${CMAKE_BINARY_DIR}/testLAPACKIntSize.F90"
        LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
        COMPILE_DEFINITIONS -DINT=int64
        COMPILE_OUTPUT_VARIABLE _MESSAGES8
        RUN_OUTPUT_VARIABLE _LAPACK_INTEGER8)

    if(_LAPACK_INTEGER8 MATCHES "^000000000")
        message(STATUS "Detected 8-byte integer LAPACK interface")
        set(GBTOlib_Fortran_FLAGS ${GBTOlib_Fortran_FLAGS} -Dblas64bit)
    elseif(_LAPACK_INTEGER4 MATCHES "^000000000")
        message(STATUS "Detected 4-byte integer LAPACK interface")
        set(GBTOlib_Fortran_FLAGS ${GBTOlib_Fortran_FLAGS} -Dblas32bit)
    else()
        message("${_MESSAGES4}")
        message("${_LAPACK_INTEGER4}")
        message("${_MESSAGES8}")
        message("${_LAPACK_INTEGER8}")
        message(STATUS "Unknown LAPACK integer interface, assuming default integer")
    endif()
endif()

# -------------------------------------------------------------------------------------------------------#
# Set up Doxygen documentation                                                                           #
# -------------------------------------------------------------------------------------------------------#

if(BUILD_DOC)
    find_package(Doxygen)

    if(NOT DOXYGEN_FOUND)
        message("Doxygen needs to be installed to generate the doxygen documentation")
    endif(NOT DOXYGEN_FOUND)
endif(BUILD_DOC)

# -------------------------------------------------------------------------------------------------------#
# Set up the target and its properties                                                                   #
# -------------------------------------------------------------------------------------------------------#

set(integral_library_sources
    "source/BB_shell_mixed_integrals_mod.f90"
    "source/GG_shell_mixed_integrals_mod.f90"
    "source/atomic_basis_mod.f90"
    "source/basis_data_generic_mod.f90"
    "source/blas_lapack.F90"
    "source/bspline_base.f90"
    "source/bspline_grid_mod.f90"
    "source/bto_gto_integrals_mod.f90"
    "source/bto_integrals_mod.f90"
    "source/cgto_hgp.f90"
    "source/cgto_integrals.f90"
    "source/cgto_pw_expansions_mod.f90"
    "source/common_obj.f90"
    "source/const.f90"
    "source/coupling_obj.f90"
    "source/data_file.f90"
    "source/eri_sph_coord.f90"
    "source/file_mapping.c"
    "source/file_mapping.F90"
    "source/free_scattering_mod.f90"
    "source/function_integration.f90"
    "source/general_quadrature.f90"
    "source/grid_mod.f90"
    "source/gto_routines.f90"
    "source/gxg_integrals_mod.f90"
    "source/integral_indexing.f90"
    "source/integral_storage_mod.f90"
    "source/lag_cfs.f90"
    "source/lebedev.f90"
    "source/molden_const.f90"
    "source/molden_mod.f90"
    "source/molecular_basis_mod.f90"
    "source/mpi_memory_mod.F90"
    "source/mpi_mod.F90"
    "source/orbital_routines.f90"
    "source/orthogonalization.f90"
    "source/parallel_arrays.f90"
    "source/phys_const.f90"
    "source/precisn.F90"
    "source/quadrature_module.f90"
    "source/sort.f90"
    "source/sparse_block_storage.f90"
    "source/special_functions.f90"
    "source/symmetry.f90"
    "source/ukrmol_interface.f90"
    "source/utils.f90"
    "source/wigner_cf.f90"
    "source/pco_mod.f90"
)

add_library(libGBTO ${integral_library_sources})

target_compile_options(libGBTO PUBLIC
    $<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>
    $<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>
    $<$<COMPILE_LANGUAGE:Fortran>:${GBTOlib_Fortran_FLAGS}>
)

list(REMOVE_ITEM CMAKE_C_IMPLICIT_LINK_LIBRARIES mpifort mpi) # let mpi_ilp64 come first if present

target_link_libraries(libGBTO
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

if(MPI_Fortran_FOUND)
    target_link_libraries(libGBTO MPI::MPI_Fortran)
endif()

set_target_properties(libGBTO PROPERTIES
    OUTPUT_NAME "GBTO"
    LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
)

add_executable(get_integrals       "source/programs/get_integrals_template.f90")

target_link_libraries(get_integrals       libGBTO)

set_target_properties(get_integrals       PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

set(programs
    basis_read
    orbitals_to_cube
    read_fock_blocks
    scatci_integrals
)

foreach(program ${programs})
    add_executable(${program} "source/programs/${program}.f90")
    target_link_libraries(${program} libGBTO)
    set_target_properties(${program} PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
endforeach()

install(TARGETS libGBTO get_integrals ${programs}
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gbtolib")

# -------------------------------------------------------------------------------------------------------#
# Set up CTest                                                                                           #
# -------------------------------------------------------------------------------------------------------#

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
    if(BUILD_UNIT_TESTS)
        add_subdirectory(source/tests)
    endif()
endif()

# -------------------------------------------------------------------------------------------------------#
# Build up Doxygen documentation                                                                         #
# -------------------------------------------------------------------------------------------------------#

if (DOXYGEN_FOUND)
    set(DOXYGEN_IN "${CMAKE_CURRENT_LIST_DIR}/doc/Doxyfile")

    add_custom_target(libGBTO-doc ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/doc"
        COMMENT "Generating GBTOlib documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/doc/html/"
            DESTINATION "${CMAKE_INSTALL_DOCDIR}/gbtolib-doxygen")
endif (DOXYGEN_FOUND)

install(FILES "${CMAKE_CURRENT_LIST_DIR}/doc/orbitals_to_cube"
        DESTINATION "${CMAKE_INSTALL_DOCDIR}")
