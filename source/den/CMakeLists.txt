# ------------------------------------------------------------------------------------------------
# denprop progeam

add_executable(denprop "denprop.f")

target_include_directories(denprop PUBLIC
    "${CMAKE_BINARY_DIR}/source/cdenprop/mod"
    "${CMAKE_BINARY_DIR}/source/gbtolib/mod"
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_link_libraries(denprop
    libcdenprop
    libukplus_utilities
    libGBTO
    libukplus_global
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

set_target_properties(denprop PROPERTIES
    LINK_FLAGS ${OpenMP_Fortran_FLAGS}
)

# ------------------------------------------------------------------------------------------------
# installation targets

install(TARGETS denprop
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/den")
