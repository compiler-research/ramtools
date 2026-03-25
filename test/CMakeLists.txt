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
            ROOT::TreePlayer
            ROOT::RIO
            sam_generator  
            gtest
            gtest_main
    )

    if(RAMTOOLS_BUILD_BENCHMARKS)
    
        target_link_libraries(${test_name} benchmark::benchmark)
    endif()

    target_include_directories(${test_name}
        PRIVATE
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/inc>
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/benchmark>
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_ramcore_test(ramcoretests ramcoretests.cxx)
add_ramcore_test(chromosome_split_test chromosome_split_test.cxx)

install(TARGETS ramcoretests chromosome_split_test
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_ramcore_test(samparsertest samparsertest.cxx)
