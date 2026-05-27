cmake_minimum_required(VERSION 3.5...3.31)

# Extract test name
get_filename_component(testname "${UKRmolIn_TEST_TGT}" NAME)

# Make the directory for extracted outputs
file(MAKE_DIRECTORY "${UKRmolIn_TEST_SRC}/outputs")

# Get all congen inputs
file(GLOB congen_target     "${UKRmolIn_TEST_SRC}/inputs/target.congen.*.inp")
file(GLOB congen_scattering "${UKRmolIn_TEST_SRC}/inputs/scattering.congen.*.inp")

# Sort inputs alphabetically
list(SORT congen_target)
list(SORT congen_scattering)

# Move integrals to the expected scratch file
if(EXISTS "${UKRmolIn_TEST_TGT}/moints")
    file(RENAME "${UKRmolIn_TEST_TGT}/moints" "${UKRmolIn_TEST_TGT}/fort.16")
endif()

# Target run
foreach(input ${congen_target})

    string(REGEX REPLACE "(.*\\.congen\\.)(.*)(\\.inp)" "\\2" spinsym ${input})

    # congen
    execute_process(COMMAND ${CONGEN_PROGRAM}
                    RESULT_VARIABLE errcode
                    INPUT_FILE "${UKRmolIn_TEST_SRC}/inputs/target.congen.${spinsym}.inp"
                    OUTPUT_FILE "logs/target.congen.${spinsym}.out"
                    ERROR_FILE "logs/target.congen.${spinsym}.err"
                    WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Target congen failed for ${spinsym} with error ${errcode}")
    endif()

    # scatci
    execute_process(COMMAND ${SCATCI_PROGRAM} "${UKRmolIn_TEST_SRC}/inputs/target.scatci.${spinsym}.inp"
                    RESULT_VARIABLE errcode
                    INPUT_FILE "${UKRmolIn_TEST_SRC}/inputs/target.scatci.${spinsym}.inp"
                    OUTPUT_FILE "logs/target.scatci.${spinsym}.out"
                    ERROR_FILE "logs/target.scatci.${spinsym}.err"
                    WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Target scatci failed for ${spinsym} with error ${errcode}")
    endif()

endforeach()

# Run denprop
execute_process(COMMAND ${DENPROP_PROGRAM}
                RESULT_VARIABLE errcode
                INPUT_FILE "${UKRmolIn_TEST_SRC}/inputs/target.denprop.inp"
                OUTPUT_FILE "logs/target.denprop.out"
                ERROR_FILE "logs/target.denprop.err"
                WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
if(NOT "${errcode}" STREQUAL "0")
    message(FATAL_ERROR "Denprop failed with error ${errcode}")
endif()

file(READ "${UKRmolIn_TEST_TGT}/fort.24" props)
string(REPLACE "D" "E" props "${props}")  # Sorry, dysprosium, dubnium and darmstadtium...
file(WRITE "${UKRmolIn_TEST_TGT}/outputs/target.properties.dat" "${props}")
file(APPEND "${UKRmolIn_TEST_TGT}/../output.target.properties.dat" "${testname}\n---->\n${props}<----\n")

# Compare properties
if(NUMDIFF_EXECUTABLE)
    execute_process(COMMAND "${NUMDIFF_EXECUTABLE}" -D "\\s \\s- \\n" -a 0 -r 1e-3 "${UKRmolIn_TEST_SRC}/outputs/target.properties.dat" "${UKRmolIn_TEST_TGT}/outputs/target.properties.dat"
                    RESULT_VARIABLE errcode
                    OUTPUT_FILE "logs/target.denprop.numdiff.out"
                    ERROR_FILE "logs/target.denprop.numdiff.err"
                    WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Target properties wrong!")
    endif()
endif()

# Scattering run
foreach(input ${congen_scattering})

    string(REGEX REPLACE "(.*\\.congen\\.)(.*)(\\.inp)" "\\2" spinsym ${input})

    # congen
    execute_process(COMMAND ${CONGEN_PROGRAM}
                    RESULT_VARIABLE errcode
                    INPUT_FILE "${UKRmolIn_TEST_SRC}/inputs/scattering.congen.${spinsym}.inp"
                    OUTPUT_FILE "logs/scattering.congen.${spinsym}.out"
                    ERROR_FILE "logs/scattering.congen.${spinsym}.err"
                    WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Scattering congen failed for ${spinsym} with error ${errcode}")
    endif()

    # scatci (or mpi-scatci)
    execute_process(COMMAND ${SCATCI_PROGRAM} "${UKRmolIn_TEST_SRC}/inputs/scattering.scatci.${spinsym}.inp"
                    RESULT_VARIABLE errcode
                    INPUT_FILE "${UKRmolIn_TEST_SRC}/inputs/scattering.scatci.${spinsym}.inp"
                    OUTPUT_FILE "logs/scattering.scatci.${spinsym}.out"
                    ERROR_FILE "logs/scattering.scatci.${spinsym}.err"
                    WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Scattering scatci failed for ${spinsym} with error ${errcode}")
    endif()

    # extract eigen-energies for comparison
    file(READ "${UKRmolIn_TEST_TGT}/logs/scattering.scatci.${spinsym}.out" logfile)
    string(REGEX REPLACE "(.*[\n\r])( ?EIGEN-ENERGIES[\n\r]+[-0-9 \\.\t\n\r]*)(.*)" "\\2" eigs "${logfile}")
    string(STRIP "${eigs}" eigs)
    file(WRITE "${UKRmolIn_TEST_TGT}/outputs/scattering.eigs.${spinsym}.dat" "${eigs}")
    file(APPEND "${UKRmolIn_TEST_TGT}/../output.scattering.eigs.dat" "${testname}: ${spinsym}\n---->\n${eigs}\n<----\n")
    if(NUMDIFF_EXECUTABLE)
        execute_process(COMMAND "${NUMDIFF_EXECUTABLE}" -a 0 -r 1e-3 "${UKRmolIn_TEST_SRC}/outputs/scattering.eigs.${spinsym}.dat" "${UKRmolIn_TEST_TGT}/outputs/scattering.eigs.${spinsym}.dat"
                        RESULT_VARIABLE errcode
                        OUTPUT_FILE "logs/scattering.eigs.${spinsym}.numdiff.out"
                        ERROR_FILE "logs/scattering.eigs.${spinsym}.numdiff.err"
                        WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
        if(NOT "${errcode}" STREQUAL "0")
            message(FATAL_ERROR "Scattering R-matrix poles wrong for spin-symmetry ${spinsym}!")
        endif()
    endif()

    # cdenprop (possibly parallel)
    if(EXISTS "${UKRmolIn_TEST_SRC}/inputs/scattering.cdenprop.${spinsym}.inp")
        file(COPY "${UKRmolIn_TEST_SRC}/inputs/scattering.cdenprop.${spinsym}.inp" DESTINATION "${UKRmolIn_TEST_TGT}")
        file(RENAME "${UKRmolIn_TEST_TGT}/scattering.cdenprop.${spinsym}.inp" "${UKRmolIn_TEST_TGT}/inp")
        execute_process(COMMAND ${CDENPROP_PROGRAM}
                        RESULT_VARIABLE errcode
                        INPUT_FILE "${UKRmolIn_TEST_SRC}/inputs/scattering.cdenprop.${spinsym}.inp"
                        OUTPUT_FILE "logs/scattering.cdenprop.${spinsym}.out"
                        ERROR_FILE "logs/scattering.cdenprop.${spinsym}.err"
                        WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
        if(NOT "${errcode}" STREQUAL "0")
            message(FATAL_ERROR "Scattering cdenprop failed for ${spinsym} with error ${errcode}")
        endif()
    endif()

endforeach()

# Postprocessing run
if(EXISTS "${UKRmolIn_TEST_SRC}/inputs/scattering.cdenprop_all.inp")
    execute_process(COMMAND ${CDENPROP_ALL_PROGRAM}
                    RESULT_VARIABLE errcode
                    INPUT_FILE "${UKRmolIn_TEST_SRC}/inputs/scattering.cdenprop_all.inp"
                    OUTPUT_FILE "logs/scattering.cdenprop_all.out"
                    ERROR_FILE "logs/scattering.cdenprop_all.err"
                    WORKING_DIRECTORY "${UKRmolIn_TEST_TGT}")
    if(NOT "${errcode}" STREQUAL "0")
        message(FATAL_ERROR "Scattering cdenprop_all failed with error ${errcode}")
    endif()
endif()
