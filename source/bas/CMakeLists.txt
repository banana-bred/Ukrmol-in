set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(CMAKE_Fortran_MODULE_DIRECTORY "mod")

# ---------------------------------------------------------------------------- #
# Programs GTOBAS, NUMCBAS, NUMSBAS

add_library(numcbas_data STATIC "numcbas_data.f")

add_executable(gtobas "gtobas.f")
add_executable(numcbas "numcbas.f")
add_executable(numsbas "numsbas.f90")

target_include_directories(numcbas_data PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_include_directories(gtobas PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_include_directories(numcbas PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_include_directories(numsbas PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_link_libraries(numcbas_data libukplus_global)

target_link_libraries(gtobas  numcbas_data libukplus_utilities libukplus_global)
target_link_libraries(numcbas numcbas_data libukplus_utilities libukplus_global)
target_link_libraries(numsbas numcbas_data libukplus_utilities libukplus_global)

install(TARGETS gtobas numcbas numsbas
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/bas")
