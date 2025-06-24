cmake_minimum_required(VERSION 3.16)
project(ramtools VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(RAMTOOLS_BUILD_TESTS "Build unit tests" ON)
option(RAMTOOLS_BUILD_BENCHMARKS "Build benchmarks" ON)

find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO Tree ROOTNTuple ROOTNTupleUtil)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)

add_library(ramtools_tree STATIC
    src/ttree/RAMRecord.cxx
)
target_include_directories(ramtools_tree PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(ramtools_tree PUBLIC 
    ROOT::Core 
    ROOT::RIO 
    ROOT::Tree
)

add_library(ramtools_ntuple STATIC
    src/rntuple/RAMNTupleRecord.cxx
)
target_include_directories(ramtools_ntuple PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(ramtools_ntuple PUBLIC 
    ROOT::Core 
    ROOT::RIO 
    ROOT::Tree 
    ROOT::ROOTNTuple 
    ROOT::ROOTNTupleUtil
)
install(TARGETS ramtools_tree ramtools_ntuple
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY inc/ DESTINATION include)

if(RAMTOOLS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

if(RAMTOOLS_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()

