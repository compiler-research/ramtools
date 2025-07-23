cmake_minimum_required(VERSION 3.16)
<<<<<<< HEAD
project(ramtools VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
=======

project(ramtools LANGUAGES CXX)
>>>>>>> develop

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

<<<<<<< HEAD
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
=======
find_package(ROOT 6.26 REQUIRED
  COMPONENTS Core RIO Tree ROOTNTuple)

set(RNTUPLE_LIBS ROOT::ROOTNTuple ROOT::ROOTNTupleUtil)

include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/rntuple)

ROOT_STANDARD_LIBRARY_PACKAGE(ramrecord
  SOURCES  ramrecord.C
  HEADERS  ramrecord.h
  LIBRARIES ROOT::Core ROOT::RIO ROOT::Tree
)

set(RAMTOOLS_ROOT_LIBS ROOT::Core ROOT::RIO ROOT::Tree)

ROOT_STANDARD_LIBRARY_PACKAGE(ramntuplerecord
  SOURCES  rntuple/RAMntuplerecord.C
  HEADERS  rntuple/RAMNTupleRecord.h
  LIBRARIES ROOT::Core ROOT::RIO ROOT::Tree ${RNTUPLE_LIBS}
)

ROOT_EXECUTABLE(samtoramntuple
  rntuple/samtoramntuple.C
  LIBRARIES ramntuplerecord ${RAMTOOLS_ROOT_LIBS} ${RNTUPLE_LIBS})

include(CTest)
file(GLOB RAMTESTS tests/test_*.C)
foreach(testfile ${RAMTESTS})
  get_filename_component(testname ${testfile} NAME_WE)
  add_test(NAME ${testname}
           COMMAND root -l -b -q ${testfile}+)
  set_tests_properties(${testname} PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach() 
>>>>>>> develop

