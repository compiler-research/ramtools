list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

add_executable(generate_test_sam generate_test_sam.cxx)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/samexample.sam
    COMMAND generate_test_sam ${CMAKE_CURRENT_BINARY_DIR}/samexample.sam 100
    DEPENDS generate_test_sam
    COMMENT "Generating test SAM file"
)

add_custom_target(generate_sam_data ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/samexample.sam
)

function(add_ramcore_test test_name)
    ROOT_EXECUTABLE(${test_name}
        ${ARGN}
        LIBRARIES
            ramcore
            ROOT::Core
            ROOT::Tree
            ROOT::ROOTNTuple
            ROOT::RIO
            gtest
            gtest_main
    )
    
    add_dependencies(${test_name} generate_sam_data)
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_ramcore_test(ramcoretests ramcoretests.cxx)

