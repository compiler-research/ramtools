list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

function(add_ramcore_test test_name)
    ROOT_EXECUTABLE(${test_name}
        ${ARGN}
        LIBRARIES
            ramcore
            sam_generator
            benchmark::benchmark
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
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/benchmark>
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_ramcore_test(ramcoretests ramcoretests.cxx)

install(TARGETS ramcoretests
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

