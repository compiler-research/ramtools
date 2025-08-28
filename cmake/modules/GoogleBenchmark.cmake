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
endif()

# Resolve canonical target names to variables without creating aliases
if(TARGET benchmark::benchmark)
    set(BENCHMARK_LIB_TARGET benchmark::benchmark CACHE INTERNAL "Google Benchmark library target")
elseif(TARGET benchmark)
    set(BENCHMARK_LIB_TARGET benchmark CACHE INTERNAL "Google Benchmark library target")
endif()

# Verify benchmark target exists
if(NOT DEFINED BENCHMARK_LIB_TARGET)
    message(FATAL_ERROR "Google Benchmark target not found (neither benchmark::benchmark nor benchmark)")
endif() 

