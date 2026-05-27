if(BUILD_TESTING)
    enable_testing()

    # find numdiff for test suite
    if(WITH_NUMDIFF)
        find_program(NUMDIFF_EXECUTABLE numdiff
                     DOC "Path to the numdiff program for use in the test suite.")
        if(NUMDIFF_EXECUTABLE)
            message(STATUS "Found numdiff: ${NUMDIFF_EXECUTABLE}")
        else()
            message(FATAL_ERROR "Numdiff requested but not found!")
        endif()
    else()
        unset(NUMDIFF_EXECUTABLE CACHE)
    endif()
endif()

# -------------------------------------------------------------------------------------------------------#
# Helper macro for addition of tests                                                                     #
# -------------------------------------------------------------------------------------------------------#

macro(ukrmolin_add_test testname)
    if(NOT "-Dusempi" IN_LIST GBTOlib_Fortran_FLAGS)
        # add serial test for not MPI-enabled executables
        add_test(NAME "ukrmolin_${testname}_serial"
                 COMMAND cmake -D "CONGEN_PROGRAM=$<TARGET_FILE:congen>"
                               -D "SCATCI_PROGRAM=$<TARGET_FILE:scatci>"
                               -D "DENPROP_PROGRAM=$<TARGET_FILE:denprop>"
                               -D "CDENPROP_PROGRAM=$<TARGET_FILE:cdenprop>"
                               -D "CDENPROP_ALL_PROGRAM=$<TARGET_FILE:cdenprop_all>"
                               -D "NUMDIFF_EXECUTABLE=${NUMDIFF_EXECUTABLE}"
                               -D "UKRmolIn_TEST_SRC=${CMAKE_CURRENT_LIST_DIR}/${testname}"
                               -D "UKRmolIn_TEST_TGT=${CMAKE_BINARY_DIR}/Testing/${testname}_serial"
                               -P "${CMAKE_CURRENT_LIST_DIR}/TestDriver.cmake")
        set_property(TEST "ukrmolin_${testname}_serial" APPEND PROPERTY DEPENDS "integrals_${testname}_serial")
    else()
        # add serial & parallel test for MPI-enabled executables
        add_test(NAME "ukrmolin_${testname}_serial"
                 COMMAND cmake -D "CONGEN_PROGRAM=$<TARGET_FILE:congen>"
                               -D "SCATCI_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:scatci>;${MPIEXEC_POSTFLAGS}"
                               -D "DENPROP_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:denprop>;${MPIEXEC_POSTFLAGS}"
                               -D "CDENPROP_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:cdenprop>;${MPIEXEC_POSTFLAGS}"
                               -D "CDENPROP_ALL_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:cdenprop_all>;${MPIEXEC_POSTFLAGS}"
                               -D "NUMDIFF_EXECUTABLE=${NUMDIFF_EXECUTABLE}"
                               -D "UKRmolIn_TEST_SRC=${CMAKE_CURRENT_LIST_DIR}/${testname}"
                               -D "UKRmolIn_TEST_TGT=${CMAKE_BINARY_DIR}/Testing/${testname}_serial"
                               -P "${CMAKE_CURRENT_LIST_DIR}/TestDriver.cmake")
        add_test(NAME "ukrmolin_${testname}_parallel"
                 COMMAND cmake -D "CONGEN_PROGRAM=$<TARGET_FILE:congen>"
                               -D "SCATCI_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};2;${MPIEXEC_PREFLAGS};$<TARGET_FILE:mpi-scatci>;${MPIEXEC_POSTFLAGS}"
                               -D "DENPROP_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:denprop>;${MPIEXEC_POSTFLAGS}"
                               -D "CDENPROP_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:cdenprop>;${MPIEXEC_POSTFLAGS}"
                               -D "CDENPROP_ALL_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:cdenprop_all>;${MPIEXEC_POSTFLAGS}"
                               -D "NUMDIFF_EXECUTABLE=${NUMDIFF_EXECUTABLE}"
                               -D "UKRmolIn_TEST_SRC=${CMAKE_CURRENT_LIST_DIR}/${testname}"
                               -D "UKRmolIn_TEST_TGT=${CMAKE_BINARY_DIR}/Testing/${testname}_parallel"
                               -P "${CMAKE_CURRENT_LIST_DIR}/TestDriver.cmake")
        set_property(TEST "ukrmolin_${testname}_serial"   APPEND PROPERTY DEPENDS "integrals_${testname}_serial")
        set_property(TEST "ukrmolin_${testname}_parallel" APPEND PROPERTY DEPENDS "integrals_${testname}_parallel")
    endif()
endmacro()

# -------------------------------------------------------------------------------------------------------#
# Basic scattering test suite                                                                            #
# -------------------------------------------------------------------------------------------------------#

foreach(group C1 Cs Ci C2 C2v C2h D2h)
    foreach(model HF CAS)
        foreach(orientation 1 2 3)
            ukrmolin_add_test(${group}_target_${model}${orientation})
        endforeach()
    endforeach()
endforeach()

foreach(group C1 Cs Ci C2 C2v C2h D2h)
    foreach(model SEP CC)
        ukrmolin_add_test(${group}_scattering_${model})
    endforeach()
endforeach()

# -------------------------------------------------------------------------------------------------------#
# Additional feature tests                                                                               #
# -------------------------------------------------------------------------------------------------------#

ukrmolin_add_test(C2v_scattering_PCCHF_A)
ukrmolin_add_test(C2v_scattering_PCCHF_B)
ukrmolin_add_test(C2v_photoionization_PCCHF_A)
ukrmolin_add_test(C2v_photoionization_PCCHF_B)
ukrmolin_add_test(C2v_photoionization_CC_GTO)
ukrmolin_add_test(C2v_photoionization_CC_Mixed)
ukrmolin_add_test(D2h_photoionization_CC)
ukrmolin_add_test(D2h_photoionization_CC_atomic)
ukrmolin_add_test(D2h_rmt_data_CC)
ukrmolin_add_test(D2h_scattering_ECP)
ukrmolin_add_test(D2h_scattering_ECP_atomic_BTO)
ukrmolin_add_test(D2h_scattering_Positron)
ukrmolin_add_test(D2h_scattering_Positron_PCO)
ukrmolin_add_test(C2v_scattering_PCO)
