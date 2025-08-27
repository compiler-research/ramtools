list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(GoogleBenchmark)

add_library(sam_generator STATIC
    generate_sam_benchmark.cxx
)

target_include_directories(sam_generator PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

ROOT_EXECUTABLE(generate_sam_benchmark
    generate_sam_benchmark.cxx
    LIBRARIES
        benchmark::benchmark
)

install(TARGETS generate_sam_benchmark
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

