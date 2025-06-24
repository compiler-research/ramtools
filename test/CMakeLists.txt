include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(CountEntriesTest CountEntriesTest.cxx)

target_link_libraries(CountEntriesTest PRIVATE 
    ${ROOT_LIBRARIES}
    ROOT::Core 
    ROOT::RIO 
    ROOT::Tree 
    ROOT::ROOTNTuple
    GTest::gtest 
    GTest::gtest_main
)

add_test(NAME CountEntries COMMAND CountEntriesTest)

add_executable(RamviewCountsTest RamviewCountsTest.cxx)
target_link_libraries(RamviewCountsTest PRIVATE 
    ${ROOT_LIBRARIES}
    ROOT::Core 
    ROOT::RIO 
    ROOT::Tree 
    ROOT::ROOTNTuple
    GTest::gtest 
    GTest::gtest_main
)

add_test(NAME RamviewCounts COMMAND RamviewCountsTest)

set_tests_properties(CountEntries PROPERTIES 
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tools
)

set_tests_properties(RamviewCounts PROPERTIES 
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tools
)

configure_file(${CMAKE_SOURCE_DIR}/data/samexample.sam 
               ${CMAKE_SOURCE_DIR}/tools/samexample.sam COPYONLY)

