cmake_minimum_required(VERSION 3.16)
project(ramtools VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(RAMTOOLS_BUILD_TESTS "Build unit tests" ON)
option(RAMTOOLS_BUILD_BENCHMARKS "Build benchmarks" ON)
option(RAMTOOLS_BUILD_TOOLS "Build tools" ON)
option(ENABLE_COVERAGE "Enable code coverage reporting" OFF)

if(ENABLE_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Code coverage enabled")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
    else()
        message(WARNING "Code coverage requires GCC or Clang")
    endif()
endif()

set(ROOT_MIN_VERSION 6.26)
find_package(ROOT ${ROOT_MIN_VERSION} REQUIRED COMPONENTS Core RIO Tree ROOTNTuple ROOTNTupleUtil)

include(${ROOT_USE_FILE})

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_INCLUDEDIR include)

ROOT_STANDARD_LIBRARY_PACKAGE(ramcore
    HEADERS
        inc/ttree/RAMRecord.h
        inc/ttree/Utils.h
        inc/ttree/CigarOps.h
        inc/rntuple/RAMNTupleRecord.h
        inc/ramcore/SamParser.h
        inc/ramcore/SamToTTree.h
        inc/ramcore/SamToNTuple.h
        inc/ramcore/RAMNTupleView.h
    SOURCES
        src/ttree/RAMRecord.cxx
        src/rntuple/RAMNTupleRecord.cxx
        src/ramcore/SamParser.cxx
        src/ramcore/SamToTTree.cxx
        src/ramcore/SamToNTuple.cxx
        src/ramcore/RAMNTupleView.cxx
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

target_include_directories(ramcore
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(DIRECTORY inc/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(RAMTOOLS_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

if(RAMTOOLS_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()

if(RAMTOOLS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
    
    if(ENABLE_COVERAGE)
        add_custom_target(coverage
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
            COMMAND gcovr -r ${CMAKE_SOURCE_DIR} --html --html-details -o ${CMAKE_BINARY_DIR}/coverage.html
            COMMAND gcovr -r ${CMAKE_SOURCE_DIR} --xml-pretty --xml coverage.xml
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Running tests and generating code coverage report..."
        )
    endif()
endif()

