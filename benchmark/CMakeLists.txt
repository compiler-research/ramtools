list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(GoogleBenchmark)

add_library(sam_generator STATIC 
    generate_sam_benchmark.cxx
)

target_include_directories(sam_generator PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(sam_generator PRIVATE 
    benchmark::benchmark
)

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

install(TARGETS sam_to_ram_benchmark conversion_time_benchmark region_query_benchmark
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_custom_target(benchmark
    COMMAND ${CMAKE_COMMAND} -E echo "=== SAM to RAM Benchmark ==="
    COMMAND sam_to_ram_benchmark
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "=== Conversion Time Benchmark ==="
    COMMAND conversion_time_benchmark
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "=== Region Query Benchmark ==="
    COMMAND region_query_benchmark
    DEPENDS sam_to_ram_benchmark conversion_time_benchmark region_query_benchmark
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running all RAM tools benchmarks"
)

