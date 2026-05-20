list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

if(RAMTOOLS_BUILD_TESTS OR RAMTOOLS_BUILD_BENCHMARKS)
    add_library(sam_generator STATIC 
        generate_sam.cxx  
    )

    target_include_directories(sam_generator PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

if(RAMTOOLS_BUILD_BENCHMARKS)
    include(GoogleBenchmark)
    
    add_library(ramtools_views STATIC
        ${CMAKE_SOURCE_DIR}/tools/ramview.cxx
        ${CMAKE_SOURCE_DIR}/tools/ramntupleview.cxx
    )

    target_link_libraries(ramtools_views PRIVATE 
        ramcore 
        ROOT::TreePlayer
    )

    ROOT_EXECUTABLE(sam_to_ram_benchmark
        sam_to_ram_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
    )

    ROOT_EXECUTABLE(conversion_time_benchmark
        conversion_time_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
    )

    ROOT_EXECUTABLE(region_query_benchmark
        region_query_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
            ramtools_views
    )

    ROOT_EXECUTABLE(chromosome_split_benchmark
        chromosome_split_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
    )

    install(TARGETS sam_to_ram_benchmark conversion_time_benchmark region_query_benchmark chromosome_split_benchmark
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    # Data for region_query_benchmark: the 1000G chr11 sample converted to both
    # RAM formats. Absolute paths so the benchmark finds them regardless of CWD.
    set(RAMBENCH_DATA_DIR "${CMAKE_SOURCE_DIR}/data"
        CACHE PATH "Directory for region_query_benchmark input data")
    set(RAMBENCH_TTREE "${RAMBENCH_DATA_DIR}/input_ttree.root")
    set(RAMBENCH_RNTUPLE "${RAMBENCH_DATA_DIR}/input_rntuple.root")

    if(TARGET samtoram AND TARGET samtoramntuple)
        add_custom_target(benchmark_data
            COMMAND ${CMAKE_COMMAND} -E env
                BUILD_DIR=${CMAKE_BINARY_DIR}
                DATA_DIR=${RAMBENCH_DATA_DIR}
                bash ${CMAKE_SOURCE_DIR}/scripts/prepare_region_query_data.sh
            DEPENDS samtoram samtoramntuple
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            USES_TERMINAL
            COMMENT "Preparing 1000G chr11 data for region_query_benchmark"
        )
    else()
        # Tools disabled: keep the target name valid but explain why it no-ops.
        add_custom_target(benchmark_data
            COMMAND ${CMAKE_COMMAND} -E echo
                "benchmark_data needs RAMTOOLS_BUILD_TOOLS=ON (samtoram/samtoramntuple)"
        )
    endif()

    add_custom_target(benchmark_region_query
        COMMAND ${CMAKE_COMMAND} -E env
            RAMBENCH_TTREE=${RAMBENCH_TTREE}
            RAMBENCH_RNTUPLE=${RAMBENCH_RNTUPLE}
            $<TARGET_FILE:region_query_benchmark>
        DEPENDS region_query_benchmark benchmark_data
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        USES_TERMINAL
        COMMENT "Running region query benchmark on 1000G chr11 data"
    )

    add_custom_target(benchmark
        COMMAND ${CMAKE_COMMAND} -E echo "=== SAM to RAM Benchmark ==="
        COMMAND sam_to_ram_benchmark
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "=== Conversion Time Benchmark ==="
        COMMAND conversion_time_benchmark
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "=== Region Query Benchmark ==="
        COMMAND ${CMAKE_COMMAND} -E env
            RAMBENCH_TTREE=${RAMBENCH_TTREE}
            RAMBENCH_RNTUPLE=${RAMBENCH_RNTUPLE}
            $<TARGET_FILE:region_query_benchmark>
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "=== Chromosome Split Benchmark ==="
        COMMAND chromosome_split_benchmark
        DEPENDS sam_to_ram_benchmark conversion_time_benchmark
                region_query_benchmark chromosome_split_benchmark benchmark_data
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        USES_TERMINAL
        COMMENT "Running all RAM tools benchmarks"
    )
endif()

