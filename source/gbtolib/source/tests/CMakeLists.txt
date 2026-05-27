set(unit_tests
    is_centre_symmetric
    multiway_merge_float_int
    xlm_mgvn
)

foreach(test_name ${unit_tests})

    add_executable(unit_test_${test_name} ${test_name}.f90)

    target_include_directories(unit_test_${test_name}
        PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/../../${CMAKE_Fortran_MODULE_DIRECTORY}"
    )

    set_target_properties(unit_test_${test_name} PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
    target_link_libraries(unit_test_${test_name} libGBTO)

    add_test(NAME unit_test_${test_name} COMMAND unit_test_${test_name})

endforeach()
