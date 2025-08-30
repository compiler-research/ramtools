include(ExternalProject)

set(GBENCHMARK_PREFIX "${CMAKE_BINARY_DIR}/googlebenchmark-prefix")

ExternalProject_Add(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.4
    GIT_SHALLOW 1
    EXCLUDE_FROM_ALL 1
    CMAKE_ARGS 
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBENCHMARK_ENABLE_TESTING=OFF
        -DBENCHMARK_ENABLE_GTEST_TESTS=OFF
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS "${GBENCHMARK_PREFIX}/src/googlebenchmark-build/src/${CMAKE_STATIC_LIBRARY_PREFIX}benchmark${CMAKE_STATIC_LIBRARY_SUFFIX}"
)

ExternalProject_Get_Property(googlebenchmark source_dir binary_dir)
file(MAKE_DIRECTORY ${source_dir}/include)

add_library(benchmark::benchmark IMPORTED STATIC GLOBAL)
set_target_properties(benchmark::benchmark PROPERTIES
    IMPORTED_LOCATION ${binary_dir}/src/${CMAKE_STATIC_LIBRARY_PREFIX}benchmark${CMAKE_STATIC_LIBRARY_SUFFIX}
    INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
)
add_dependencies(benchmark::benchmark googlebenchmark)

