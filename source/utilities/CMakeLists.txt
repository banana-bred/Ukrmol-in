# --------------------------------------------------------------------------------
# global utils

add_library(libukplus_utilities "global_utils.f90" "integer_packing.f90")

target_include_directories(libukplus_utilities PUBLIC "${CMAKE_BINARY_DIR}/source/global/mod")

set_target_properties(libukplus_utilities PROPERTIES OUTPUT_NAME ukplus_utilities)

target_link_libraries(libukplus_utilities libukplus_global)

# --------------------------------------------------------------------------------
# installation targets

install(TARGETS libukplus_utilities
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/utilities")
