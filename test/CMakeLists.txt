list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

function(add_ramcore_test test_name)
    ROOT_EXECUTABLE(${test_name}
        ${ARGN}
        LIBRARIES
            ramcore
            ROOT::Core
            ROOT::ROOTNTuple
            ROOT::RIO
            sam_generator
            gtest
            gtest_main
    )

    if(RAMTOOLS_BUILD_BENCHMARKS)
        target_link_libraries(${test_name} PUBLIC benchmark::benchmark)
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
add_ramcore_test(bam_conversion bam_conversion.cxx)
add_ramcore_test(parallel_conversion_test parallel_conversion_test.cxx)
add_ramcore_test(quality_blocks_test quality_blocks_test.cxx)

if(TARGET samtoramntuple)
    add_test(NAME samtoramntuple_rejects_bad_compression COMMAND samtoramntuple in.sam -compression asd)
    add_test(NAME bamtoramntuple_rejects_bad_compression COMMAND bamtoramntuple in.bam -compression 999)
    set_tests_properties(samtoramntuple_rejects_bad_compression bamtoramntuple_rejects_bad_compression
        PROPERTIES PASS_REGULAR_EXPRESSION "invalid -compression value")
    add_test(NAME samtoramntuple_rejects_bad_threads COMMAND samtoramntuple in.sam -threads 0)
    add_test(NAME bamtoramntuple_rejects_bad_threads COMMAND bamtoramntuple in.bam -threads x)
    # A number followed by garbage must not be silently truncated to the number.
    add_test(NAME samtoramntuple_rejects_trailing_threads COMMAND samtoramntuple in.sam -threads 3sdfd)
    add_test(NAME bamtoramntuple_rejects_trailing_threads COMMAND bamtoramntuple in.bam -threads 3sdfd)
    set_tests_properties(samtoramntuple_rejects_bad_threads bamtoramntuple_rejects_bad_threads
        samtoramntuple_rejects_trailing_threads bamtoramntuple_rejects_trailing_threads
        PROPERTIES PASS_REGULAR_EXPRESSION "invalid -threads value")
endif()

if(TARGET ramdump)
    add_ramcore_test(ramdump_test ramdump_test.cxx)
    add_dependencies(ramdump_test samtoramntuple ramdump)
    target_compile_definitions(ramdump_test PRIVATE
        SAMTORAMNTUPLE_BIN="$<TARGET_FILE:samtoramntuple>"
        RAMDUMP_BIN="$<TARGET_FILE:ramdump>")
endif()

install(TARGETS ramcoretests chromosome_split_test bam_conversion parallel_conversion_test quality_blocks_test
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
