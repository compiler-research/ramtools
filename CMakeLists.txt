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

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_INCLUDEDIR include)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/inc)

ROOT_STANDARD_LIBRARY_PACKAGE(RAMCore
    HEADERS
        ttree/RAMRecord.h
        ttree/Utils.h
        ttree/CigarOps.h
        rntuple/RAMNTupleRecord.h
        RAMCore/SamParser.h
        RAMCore/SamToTTree.h
        RAMCore/SamToNTuple.h
    SOURCES
        src/ttree/RAMRecord.cxx
        src/rntuple/RAMNTupleRecord.cxx
        src/RAMCore/SamParser.cxx
        src/RAMCore/SamToTTree.cxx
        src/RAMCore/SamToNTuple.cxx
    LINKDEF
        inc/ttree/LinkDef.h
    DEPENDENCIES
        ROOT::Core
        ROOT::RIO
        ROOT::Tree
        ROOT::ROOTNTuple
        ROOT::ROOTNTupleUtil
    INSTALL_OPTIONS
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY inc/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(RAMTOOLS_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

if(RAMTOOLS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

