list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(GoogleBenchmark)

add_library(sam_generator STATIC
    generate_sam_benchmark.cxx
)

target_include_directories(sam_generator PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(sam_generator PRIVATE benchmark::benchmark)
add_dependencies(sam_generator benchmark::benchmark)

ROOT_EXECUTABLE(generate_sam_benchmark
    generate_sam_benchmark.cxx
    LIBRARIES
        benchmark::benchmark
)

install(TARGETS generate_sam_benchmark
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_custom_target(sam-benchmark
    COMMAND generate_sam_benchmark
    DEPENDS generate_sam_benchmark
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running SAM generation benchmarks"
)

