list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

add_executable(CountEntriesTest CountEntriesTest.cxx)
target_link_libraries(CountEntriesTest PRIVATE
    ramtools_tree
    ramtools_ntuple
    gtest
    gtest_main
)
add_test(NAME CountEntriesTest COMMAND CountEntriesTest)
set_tests_properties(CountEntriesTest PROPERTIES
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tools
    ENVIRONMENT "ROOT_INCLUDE_PATH=${CMAKE_SOURCE_DIR}/inc;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}"
)

add_executable(RamviewCountsTest RamviewCountsTest.cxx)
target_link_libraries(RamviewCountsTest PRIVATE
    ramtools_tree
    ramtools_ntuple
    gtest
    gtest_main
)
add_test(NAME RamviewCountsTest COMMAND RamviewCountsTest)
set_tests_properties(RamviewCountsTest PROPERTIES
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tools
    ENVIRONMENT "ROOT_INCLUDE_PATH=${CMAKE_SOURCE_DIR}/inc;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}"
)

configure_file(${CMAKE_SOURCE_DIR}/data/samexample.sam
               ${CMAKE_SOURCE_DIR}/tools/samexample.sam COPYONLY)