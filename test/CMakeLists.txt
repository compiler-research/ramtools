list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

function(add_ramtools_unittest test_name)
    add_executable(${test_name} ${ARGN})
    
    target_link_libraries(${test_name} PRIVATE
        ramtools
        gtest
        gtest_main
    )
    
    add_test(NAME ${test_name} COMMAND ${test_name})
    
    set_tests_properties(${test_name} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tools
        ENVIRONMENT "ROOT_INCLUDE_PATH=${CMAKE_SOURCE_DIR}/inc;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}"
    )
endfunction()

add_ramtools_unittest(CountEntriesTest CountEntriesTest.cxx)
add_ramtools_unittest(RamviewCountsTest RamviewCountsTest.cxx)

configure_file(${CMAKE_SOURCE_DIR}/data/samexample.sam
               ${CMAKE_SOURCE_DIR}/tools/samexample.sam COPYONLY)

