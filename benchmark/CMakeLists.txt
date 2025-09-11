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

install(TARGETS sam_to_ram_benchmark conversion_time_benchmark
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_custom_target(benchmark
    COMMAND sam_to_ram_benchmark
    COMMAND conversion_time_benchmark
    DEPENDS sam_to_ram_benchmark conversion_time_benchmark
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

