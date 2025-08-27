include(FetchContent)

# Try to find benchmark package first
find_package(benchmark QUIET)

if(NOT benchmark_FOUND)
    message(STATUS "Google Benchmark not found, downloading and building with FetchContent")
    
    # Set options for Google Benchmark
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable gtest tests" FORCE)
    set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON CACHE BOOL "Download benchmark dependencies" FORCE)
    
    FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.8.4
        GIT_SHALLOW TRUE
    )
    
    FetchContent_MakeAvailable(googlebenchmark)
    
    # Create alias if needed
    if(TARGET benchmark AND NOT TARGET benchmark::benchmark)
        add_library(benchmark::benchmark ALIAS benchmark)
    endif()
    
    if(TARGET benchmark_main AND NOT TARGET benchmark::benchmark_main)
        add_library(benchmark::benchmark_main ALIAS benchmark_main)
    endif()
endif()

# Verify benchmark target exists
if(NOT TARGET benchmark::benchmark)
    message(FATAL_ERROR "Google Benchmark target benchmark::benchmark not found")
endif() 