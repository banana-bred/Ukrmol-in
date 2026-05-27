# ------------------------------------------------------------------------------------------------
# libukplus_global library

add_library(libukplus_global
    "algorithms.f90"
    "containers.f90"
    "global_data.f90"
    "precisn.f90"
    "swmol3_data.f90"
    "version_control.F90"
)

set_target_properties(libukplus_global PROPERTIES
    OUTPUT_NAME ukplus_global
)

if(UKRMOL_GIT_AUTH)
    set_property(SOURCE "version_control.F90" APPEND PROPERTY COMPILE_DEFINITIONS "GIT_AUTH=${UKRMOL_GIT_AUTH}")
endif()

if(UKRMOL_GIT_HASH)
    set_property(SOURCE "version_control.F90" APPEND PROPERTY COMPILE_DEFINITIONS "GIT_HASH=${UKRMOL_GIT_HASH}")
endif()

if(UKRMOL_GIT_DATE)
    set_property(SOURCE "version_control.F90" APPEND PROPERTY COMPILE_DEFINITIONS "GIT_DATE=${UKRMOL_GIT_DATE}")
endif()

# ------------------------------------------------------------------------------------------------
# installation targets

install(TARGETS libukplus_global
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/global")
