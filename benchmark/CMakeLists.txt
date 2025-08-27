find_package(benchmark REQUIRED)

ROOT_EXECUTABLE(generate_sam_benchmark
    generate_sam_benchmark.cxx
    LIBRARIES
        benchmark::benchmark
)

install(TARGETS generate_sam_benchmark
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
) 