# -------------------------------------------------------------------------------------------------------#
# UKRmol+ suite - inner region                                                                           #
# -------------------------------------------------------------------------------------------------------#
#                                                                                                        #
# Environment variables                                                                                  #
# ---------------------                                                                                  #
#                                                                                                        #
# - BLA_VENDOR              Defines the BLAS/LAPACK vendor, typically "Intel10_64ilp" for ILP64 MKL.     #
#                           In case of Intel Compiler this needs to be supplemented with the option      #
#                           "-mkl". (Only available with CMake 3.13+)                                    #
#                                                                                                        #
# - BLA_STATIC              Look for static BLAS/LAPACK libraries. In case of Intel Compiler, this needs #
#                           to be consistent with usage or absence of the flag "-static-intel".          #
#                           (Only available with CMake 3.13+)                                            #
#                                                                                                        #
# CMake options                                                                                          #
# -------------                                                                                          #
#                                                                                                        #
# - BLA_F95                 Use external blas95/lapack95 wrappers from Intel MKL (default: OFF).         #
#                                                                                                        #
# - BUILD_DOC               Generate development documentation (default: OFF).                           #
#                                                                                                        #
# - BUILD_TESTING           Enable the test suite (default: ON).                                         #
#                                                                                                        #
# - MPIEXEC_PREFLAGS        Additional flags to pass to MPI launcher when running the test suite.        #
#                                                                                                        #
# - GBTOlib_Fortran_FLAGS   Force specific compiler flags to be used by GBTOlib and the rest of the      #
#                           suite (default: determined by CMake).                                        #
#                                                                                                        #
# - WITH_GIT                Attempt to retrieve code version from the local Git repository clone         #
#                           (default: ON).                                                               #
#                                                                                                        #
# - WITH_MPI                Build with MPI support (default: ON).                                        #
#                                                                                                        #
# - WITH_SCALAPACK          Build with ScaLAPACK support (default: OFF). Use this switch when the        #
#                           ScaLAPACK library is automatically linked by the compiler.                   #
#                                                                                                        #
# - SCALAPACK_LIBRARIES     Path to ScaLAPACK libraries if not implicitly linked.                        #
#                                                                                                        #
# - ARPACK_LIBRARIES        Path to Arpack libraries.                                                    #
#                                                                                                        #
# - SLEPC_INCLUDE_DIRS      Path to SLEPc (and PetSc) include dirs. May need to include several paths    #
#                           separated by semicolons.                                                     #
#                                                                                                        #
# - SLEPC_LIBRARIES         Path to SLEPc libraries (and PetSc, if needed for linking).                  #
#                                                                                                        #
# - WITH_MOLPRO             Use Molpro with test suite (default: OFF).                                   #
#                                                                                                        #
# - WITH_PSI4               Use Psi4 with test suite (default: OFF; Molpro has priority if both set).    #
#                                                                                                        #
# - WITH_NUMDIFF            Use numdiff program for machine comparison of test suite (default: OFF).     #
#                                                                                                        #
# - NUMDIFF_EXECUTABLE      Path to the "numdiff" program for use in the test suite (optional).          #
#                                                                                                        #
# -------------------------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.5...3.31)

project(UKRmol+ C Fortran)

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

option(WITH_GIT         "Add Git revision to program header (requires Git)" ON)
option(WITH_PSI4        "Use the Psi4 quantum chemistry package for tests." OFF)
option(WITH_MOLPRO      "Use Molpro quantum chemistry package for tests."   OFF)
option(WITH_NUMDIFF     "Use numdiff for machine comparison of tests."      OFF)
option(BUILD_DOC        "Build documentation"                               OFF)
option(BUILD_TESTING    "Enable the test suite"                             ON)

# define the GBTO flags cache variable that can be forced by the user
set(GBTOlib_Fortran_FLAGS "${GBTOlib_Fortran_FLAGS}"
    CACHE STRING "MPI-related compiler parameters for GBTOlib")

# -------------------------------------------------------------------------------------------------------#
# Set up Doxygen documentation                                                                           #
# -------------------------------------------------------------------------------------------------------#

# check that Doxygen exists when documentation is requested
if(BUILD_DOC)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        message(STATUS "Doxygen needs to be installed to generate the doxygen documentation")
    endif(NOT DOXYGEN_FOUND)
endif(BUILD_DOC)

# -------------------------------------------------------------------------------------------------------#
# Find the required libraries                                                                            #
# -------------------------------------------------------------------------------------------------------#

# most codes can make use of OpenMP
find_package(OpenMP REQUIRED)

# find BLAS + LAPACK
if(NOT "$ENV{BLA_STATIC}" STREQUAL "")
    set(BLA_STATIC ON)
endif()
find_package(LAPACK REQUIRED)

# -------------------------------------------------------------------------------------------------------#
# Project configuration                                                                                  #
# -------------------------------------------------------------------------------------------------------#

# collect infomation about version of the code if possible
if (WITH_GIT)
    find_package(Git)
    if (Git_FOUND)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:"%an"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE UKRMOL_GIT_AUTH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:"%h"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE UKRMOL_GIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:"%ad"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE UKRMOL_GIT_DATE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" describe --tags HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE UKRMOL_GIT_TAGS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif (Git_FOUND)
endif (WITH_GIT)

# override the data from the .release tag file if present
if(EXISTS ${CMAKE_SOURCE_DIR}/.release)
    file(STRINGS ${CMAKE_SOURCE_DIR}/.release RELEASE_INFO)
    if(RELEASE_INFO)
        list(GET RELEASE_INFO 0 UKRMOL_GIT_HASH)
        list(GET RELEASE_INFO 1 UKRMOL_GIT_AUTH)
        list(GET RELEASE_INFO 2 UKRMOL_GIT_DATE)
    endif()
endif()

# communicate the version to the user
message(STATUS "Last UKRmol-in commit author: ${UKRMOL_GIT_AUTH}")
message(STATUS "Last UKRmol-in commit hash: ${UKRMOL_GIT_HASH}")
message(STATUS "Last UKRmol-in commit date: ${UKRMOL_GIT_DATE}")

# include UKRmol+ library targets
add_subdirectory(source/interfaces)
add_subdirectory(source/global)
add_subdirectory(source/gbtolib)
add_subdirectory(source/utilities)

# pull in the GBTO library flags if not set by user
if ("${GBTOlib_Fortran_FLAGS}" STREQUAL "")
    get_directory_property(GBTOlib_Fortran_FLAGS
        DIRECTORY source/gbtolib
        DEFINITION GBTOlib_Fortran_FLAGS)
endif()

# include UKRmol+ executable targets
add_subdirectory(source/bas)
add_subdirectory(source/cdenprop)
add_subdirectory(source/congen)
add_subdirectory(source/ci-diag)
add_subdirectory(source/den)
add_subdirectory(source/mpi-ci-diag)
add_subdirectory(source/tools)

# activate test suite
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests/suite)
endif()

# include also UKRmol-out targets if path provided (useful for joint testing)
if(UKRMOL_OUT_DIR)
    set(GBTOLIB_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/source/gbtolib/mod")
    set(GBTOLIB_LIBRARIES    "libGBTO")
    add_subdirectory("${UKRMOL_OUT_DIR}" "${CMAKE_BINARY_DIR}/ukrmol-out")
endif()

# install documentation
file(GLOB inner_docs "${PROJECT_SOURCE_DIR}/doc/*")
install(FILES ${inner_docs} DESTINATION "${CMAKE_INSTALL_DOCDIR}")

# -------------------------------------------------------------------------------------------------------#
# Define distribution release target (make dist)                                                         #
# -------------------------------------------------------------------------------------------------------#

add_custom_command(OUTPUT ukrmol-in.tar
    COMMAND git archive --format tar --prefix=ukrmol-in-${UKRMOL_GIT_TAGS}/ -o ${CMAKE_BINARY_DIR}/ukrmol-in.tar HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_command(OUTPUT gbtolib.tar
    COMMAND git archive --format tar --prefix=ukrmol-in-${UKRMOL_GIT_TAGS}/source/gbtolib/ -o ${CMAKE_BINARY_DIR}/gbtolib.tar HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/source/gbtolib/
)

file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/.release-in
    CONTENT "\"${UKRMOL_GIT_TAGS}\"\n${UKRMOL_GIT_AUTH}\n${UKRMOL_GIT_DATE}\n"
)

add_custom_target(dist-in
    DEPENDS ukrmol-in.tar gbtolib.tar .release-in
    COMMAND tar --extract --file ukrmol-in.tar
    COMMAND tar --extract --file gbtolib.tar
    COMMAND cp .release-in ukrmol-in-${UKRMOL_GIT_TAGS}/.release
    COMMAND tar --create --gzip --file ukrmol-in-${UKRMOL_GIT_TAGS}.tar.gz ukrmol-in-${UKRMOL_GIT_TAGS}
    COMMAND rm -rf ukrmol-in-${UKRMOL_GIT_TAGS}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

if(NOT TARGET dist)
    add_custom_target(dist)  # may be already defined in GBTOlib / UKRmol-out
endif()

add_dependencies(dist dist-in)
