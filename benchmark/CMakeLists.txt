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

    ROOT_EXECUTABLE(columnar_read_benchmark
        columnar_read_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
            ROOT::ROOTNTuple
            ROOT::RIO
    )

    install(TARGETS sam_to_ram_benchmark conversion_time_benchmark region_query_benchmark
            chromosome_split_benchmark columnar_read_benchmark
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    # Single-command runner: executes every benchmark and writes one JSON result file per
    # binary into the build directory (no source-tree output folder). Pass dataset paths to
    # individual binaries directly, or use scripts/run_benchmarks.sh, for custom datasets.
    add_custom_target(benchmark
        COMMAND ${CMAKE_COMMAND} -E echo "=== SAM to RAM Benchmark ==="
        COMMAND sam_to_ram_benchmark --benchmark_out=${CMAKE_CURRENT_BINARY_DIR}/sam_to_ram_benchmark.json --benchmark_out_format=json
        COMMAND ${CMAKE_COMMAND} -E echo "=== Conversion Time Benchmark ==="
        COMMAND conversion_time_benchmark --benchmark_out=${CMAKE_CURRENT_BINARY_DIR}/conversion_time_benchmark.json --benchmark_out_format=json
        COMMAND ${CMAKE_COMMAND} -E echo "=== Region Query Benchmark ==="
        COMMAND region_query_benchmark --benchmark_out=${CMAKE_CURRENT_BINARY_DIR}/region_query_benchmark.json --benchmark_out_format=json
        COMMAND ${CMAKE_COMMAND} -E echo "=== Chromosome Split Benchmark ==="
        COMMAND chromosome_split_benchmark --benchmark_out=${CMAKE_CURRENT_BINARY_DIR}/chromosome_split_benchmark.json --benchmark_out_format=json
        COMMAND ${CMAKE_COMMAND} -E echo "=== Columnar Read Benchmark ==="
        COMMAND columnar_read_benchmark --benchmark_out=${CMAKE_CURRENT_BINARY_DIR}/columnar_read_benchmark.json --benchmark_out_format=json
        COMMAND ${CMAKE_COMMAND} -E echo "JSON results written to ${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS sam_to_ram_benchmark conversion_time_benchmark region_query_benchmark chromosome_split_benchmark columnar_read_benchmark
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running all RAM tools benchmarks (JSON output in build dir)"
    )
endif()

