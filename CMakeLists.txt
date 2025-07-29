cmake_minimum_required(VERSION 3.16)
project(ramtools VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(RAMTOOLS_BUILD_TESTS "Build unit tests" ON)
option(RAMTOOLS_BUILD_BENCHMARKS "Build benchmarks" ON)
option(RAMTOOLS_BUILD_TOOLS "Build tools" ON)

find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO Tree ROOTNTuple ROOTNTupleUtil)

include(${ROOT_USE_FILE})

add_library(ramtools_tree SHARED
    src/ttree/RAMRecord.cxx
)

target_include_directories(ramtools_tree PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/inc)

ROOT_GENERATE_DICTIONARY(G__ramtools_tree
    ttree/RAMRecord.h
    ttree/Utils.h
    ttree/CigarOps.h
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/ttree/LinkDef.h
    MODULE ramtools_tree
)

target_link_libraries(ramtools_tree PUBLIC 
    ROOT::Core 
    ROOT::RIO 
    ROOT::Tree
)

add_library(ramtools_ntuple SHARED
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

install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/libramtools_tree_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libramtools_tree.rootmap
    DESTINATION lib
)

if(RAMTOOLS_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

if(RAMTOOLS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

if(RAMTOOLS_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()