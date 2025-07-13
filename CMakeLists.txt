cmake_minimum_required(VERSION 3.16)

project(ramtools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

