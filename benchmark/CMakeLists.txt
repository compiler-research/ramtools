list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(GoogleBenchmark)

ROOT_EXECUTABLE(generate_sam_benchmark
    generate_sam_benchmark.cxx
    LIBRARIES
        benchmark::benchmark
)

install(TARGETS generate_sam_benchmark
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

