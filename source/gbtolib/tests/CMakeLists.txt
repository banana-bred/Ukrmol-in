if(BUILD_TESTING)
    enable_testing()

    if(WITH_MOLPRO)
        find_program(MOLPRO_PROGRAM molpro)
        if(MOLPRO_PROGRAM)
            message(STATUS "Found Molpro for testing: ${MOLPRO_PROGRAM}")
        else()
            message(FATAL_ERROR "Molpro requested but not found!")
        endif()
    elseif(WITH_PSI4)
        find_program(PSI4_PROGRAM psi4)
        if(PSI4_PROGRAM)
            message(STATUS "Found Psi4 for testing: ${PSI4_PROGRAM}")
        else()
            message(FATAL_ERROR "Psi4 requested but not found!")
        endif()
    else()
        message(STATUS "No quantum chemistry program will be used for testing.")
    endif()
endif()

# -------------------------------------------------------------------------------------------------------#
# Helper macro for addition of tests                                                                     #
# -------------------------------------------------------------------------------------------------------#

macro(gbtolib_add_scatci_integrals_test testname)
    if(NOT "-Dusempi" IN_LIST GBTOlib_Fortran_FLAGS)
        # add test for not MPI-enabled executable
        add_test(NAME "integrals_${testname}_serial"
                 COMMAND cmake -D "SCATCI_INTEGRALS_PROGRAM=$<TARGET_FILE:scatci_integrals>"
                               -D "MOLPRO_PROGRAM=${MOLPRO_PROGRAM}"
                               -D "PSI4_PROGRAM=${PSI4_PROGRAM}"
                               -D "GBTOlib_TEST_SRC=${CMAKE_CURRENT_LIST_DIR}/${testname}"
                               -D "GBTOlib_TEST_TGT=${CMAKE_BINARY_DIR}/Testing/${testname}_serial"
                               -P "${CMAKE_CURRENT_LIST_DIR}/TestDriver.cmake")
    else()
        # add tests for MPI-enabled executables
        add_test(NAME "integrals_${testname}_serial"
                 COMMAND cmake -D "SCATCI_INTEGRALS_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS};$<TARGET_FILE:scatci_integrals>;${MPIEXEC_POSTFLAGS}"
                               -D "MOLPRO_PROGRAM=${MOLPRO_PROGRAM}"
                               -D "PSI4_PROGRAM=${PSI4_PROGRAM}"
                               -D "GBTOlib_TEST_SRC=${CMAKE_CURRENT_LIST_DIR}/${testname}"
                               -D "GBTOlib_TEST_TGT=${CMAKE_BINARY_DIR}/Testing/${testname}_serial"
                               -P "${CMAKE_CURRENT_LIST_DIR}/TestDriver.cmake")
        add_test(NAME "integrals_${testname}_parallel"
                 COMMAND cmake -D "SCATCI_INTEGRALS_PROGRAM=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};2;${MPIEXEC_PREFLAGS};$<TARGET_FILE:scatci_integrals>;${MPIEXEC_POSTFLAGS}"
                               -D "MOLPRO_PROGRAM=${MOLPRO_PROGRAM}"
                               -D "PSI4_PROGRAM=${PSI4_PROGRAM}"
                               -D "GBTOlib_TEST_SRC=${CMAKE_CURRENT_LIST_DIR}/${testname}"
                               -D "GBTOlib_TEST_TGT=${CMAKE_BINARY_DIR}/Testing/${testname}_parallel"
                               -P "${CMAKE_CURRENT_LIST_DIR}/TestDriver.cmake")
    endif()
endmacro()

# -------------------------------------------------------------------------------------------------------#
# Basic scattering test suite                                                                            #
# -------------------------------------------------------------------------------------------------------#

foreach(group C1 Cs Ci C2 C2v C2h D2h)
    foreach(model HF CAS)
        foreach(orientation 1 2 3)
            gbtolib_add_scatci_integrals_test(${group}_target_${model}${orientation})
        endforeach()
    endforeach()
endforeach()

foreach(group C1 Cs Ci C2 C2v C2h D2h)
    foreach(model SEP CC)
        gbtolib_add_scatci_integrals_test(${group}_scattering_${model})
    endforeach()
endforeach()

# -------------------------------------------------------------------------------------------------------#
# Additional feature tests                                                                               #
# -------------------------------------------------------------------------------------------------------#

gbtolib_add_scatci_integrals_test(C2v_scattering_PCCHF_A)
gbtolib_add_scatci_integrals_test(C2v_scattering_PCCHF_B)
gbtolib_add_scatci_integrals_test(C2v_photoionization_PCCHF_A)
gbtolib_add_scatci_integrals_test(C2v_photoionization_PCCHF_B)
gbtolib_add_scatci_integrals_test(C2v_photoionization_CC_GTO)
gbtolib_add_scatci_integrals_test(C2v_photoionization_CC_Mixed)
gbtolib_add_scatci_integrals_test(D2h_photoionization_CC)
gbtolib_add_scatci_integrals_test(D2h_photoionization_CC_atomic)
gbtolib_add_scatci_integrals_test(D2h_rmt_data_CC)
gbtolib_add_scatci_integrals_test(D2h_scattering_ECP)
gbtolib_add_scatci_integrals_test(D2h_scattering_ECP_atomic_BTO)
gbtolib_add_scatci_integrals_test(D2h_scattering_Positron)
gbtolib_add_scatci_integrals_test(D2h_scattering_Positron_PCO)
gbtolib_add_scatci_integrals_test(C2v_scattering_PCO)
gbtolib_add_scatci_integrals_test(D2h_minimal_scattering_integrals)
