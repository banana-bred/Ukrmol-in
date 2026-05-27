cmake_minimum_required(VERSION 3.5...3.31)

# Copy the test sources to the working directory
file(MAKE_DIRECTORY "${GBTOlib_TEST_TGT}")
file(MAKE_DIRECTORY "${GBTOlib_TEST_TGT}/logs")
file(GLOB_RECURSE sources "${GBTOlib_TEST_SRC}/*")
file(COPY ${sources} DESTINATION "${GBTOlib_TEST_TGT}")

# Rename input file
if(EXISTS "${GBTOlib_TEST_TGT}/target.integrals.inp")
    file(RENAME "${GBTOlib_TEST_TGT}/target.integrals.inp" "${GBTOlib_TEST_TGT}/inp")
endif()
if(EXISTS "${GBTOlib_TEST_TGT}/scattering.integrals.inp")
    file(RENAME "${GBTOlib_TEST_TGT}/scattering.integrals.inp" "${GBTOlib_TEST_TGT}/inp")
endif()

# Remove potentially existing molecular integrals file (when re-running a test)
file(REMOVE "${GBTOlib_TEST_TGT}/moints")

# Run quantum chemistry program (if available)
if(MOLPRO_PROGRAM)
    execute_process(COMMAND ${MOLPRO_PROGRAM}
                    INPUT_FILE target.molpro.inp
                    OUTPUT_FILE logs/target.molpro.out
                    ERROR_FILE logs/target.molpro.err
                    RESULT_VARIABLE errcode
                    WORKING_DIRECTORY "${GBTOlib_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Molpro failed!")
    endif()
elseif(PSI4_PROGRAM)
    execute_process(COMMAND ${PSI4_PROGRAM} --input target.psi4.inp --output logs/target.psi4.out
                    ERROR_FILE logs/target.psi4.err
                    RESULT_VARIABLE errcode
                    WORKING_DIRECTORY "${GBTOlib_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Psi4 failed!")
    endif()
endif()

# Run scatci_integrals
execute_process(COMMAND ${SCATCI_INTEGRALS_PROGRAM}
                RESULT_VARIABLE errcode
                OUTPUT_FILE logs/scatci_integrals.out
                ERROR_FILE logs/scatci_integrals.err
                WORKING_DIRECTORY "${GBTOlib_TEST_TGT}")

# When the exit code is all right, test that the log file ends with 'done:mpi_mod:mpi_mod_finalize'
if("${errcode}" STREQUAL "0")
    file(READ "${GBTOlib_TEST_TGT}/log_file.0" log)
    file(RENAME "${GBTOlib_TEST_TGT}/log_file.0" "${GBTOlib_TEST_TGT}/logs/scatci_integrals.log")
    string(FIND "${log}" "done:mpi_mod:mpi_mod_finalize" verdict)
    if(${verdict} EQUAL -1)
        set(errcode 1)
    endif()
endif()

# Throw error on failure
if(NOT "${errcode}" STREQUAL "0")
    message(FATAL_ERROR "Execution failed!")
endif()

