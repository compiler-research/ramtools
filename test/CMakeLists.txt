list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/samexample.sam
    COMMAND $<TARGET_FILE:generate_sam_benchmark> --generate ${CMAKE_CURRENT_BINARY_DIR}/samexample.sam 100
    DEPENDS generate_sam_benchmark
    COMMENT "Generating test SAM file using benchmark generator"
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
    target_include_directories(${test_name}
        PRIVATE
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/inc>
    )
    add_dependencies(${test_name} generate_sam_data)
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_ramcore_test(ramcoretests ramcoretests.cxx)

install(TARGETS ramcoretests
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

