list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

function(add_ramcore_test test_name)
    ROOT_EXECUTABLE(${test_name}
        ${ARGN}
        LIBRARIES
            ramcore
            ROOT::Core
            ROOT::Tree
            ROOT::ROOTNTuple
            gtest
            gtest_main
    )

    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_ramcore_test(ramcoretests ramcoretests.cxx)

if(EXISTS ${CMAKE_SOURCE_DIR}/test/samexample.sam)
    configure_file(${CMAKE_SOURCE_DIR}/test/samexample.sam
               ${CMAKE_CURRENT_BINARY_DIR}/samexample.sam
               COPYONLY)
endif()

