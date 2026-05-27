# ------------------------------------------------------------------------------------------------
# interface library

add_library(libukplus_interfaces "maths_interface.F90")

if(BLA_F95)
    message(STATUS "UKRmol-in will be compiled with external blas95/lapack95 wrappers")
    set_source_files_properties(maths_interface.F90 PROPERTIES COMPILE_DEFINITIONS BLA_F95=1)
endif()

set_target_properties(libukplus_interfaces PROPERTIES
    OUTPUT_NAME ukplus_interfaces
)

target_include_directories(libukplus_interfaces PUBLIC
    "${CMAKE_BINARY_DIR}/source/gbtolib/mod"
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_link_libraries(libukplus_interfaces PUBLIC
    libukplus_utilities
    libukplus_global
    libGBTO
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

# ------------------------------------------------------------------------------------------------
# installation targets

install(TARGETS libukplus_interfaces
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/interfaces")
