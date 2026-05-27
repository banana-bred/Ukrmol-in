# ------------------------------------------------------------------------------------------------
# congen program

add_executable(congen
    "source/congen_bstree.f90"
    "source/congen_data.f90"
    "source/congen_distribution.f90"
    "source/congen_driver.f90"
    "source/congen_pagecontrol.f90"
    "source/congen_projection.f90"
    "source/congen.f90"
)

target_include_directories(congen PUBLIC
    "${CMAKE_BINARY_DIR}/source/global/mod"
    "${CMAKE_BINARY_DIR}/source/utilities/mod"
)

target_link_libraries(congen PUBLIC
    libukplus_utilities
    libukplus_global
)

# ------------------------------------------------------------------------------------------------
# congen documentation

if (DOXYGEN_FOUND)
    set(DOXYGEN_IN "${CMAKE_CURRENT_LIST_DIR}/doc/Doxyfile")

    add_custom_target(congen-doc ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/doc"
        COMMENT "Generating CONGEN documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/doc/html/"
            DESTINATION "${CMAKE_INSTALL_DOCDIR}/congen-doxygen")
endif (DOXYGEN_FOUND)

# ------------------------------------------------------------------------------------------------
# installation targets

install(TARGETS congen
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_Fortran_MODULE_DIRECTORY}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/congen")
