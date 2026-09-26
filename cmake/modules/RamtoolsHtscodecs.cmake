# ramtools::fqzcomp: the system htscodecs if pkg-config finds it, otherwise the
# two source files fqzcomp needs, built from a pinned release.

pkg_check_modules(htscodecs QUIET IMPORTED_TARGET htscodecs)

if(htscodecs_FOUND)
  add_library(ramtools_fqzcomp INTERFACE)
  target_link_libraries(ramtools_fqzcomp INTERFACE PkgConfig::htscodecs)
  add_library(ramtools::fqzcomp ALIAS ramtools_fqzcomp)
  message(STATUS "fqzcomp: system htscodecs ${htscodecs_VERSION}")
  return()
endif()

include(FetchContent)
FetchContent_Declare(htscodecs
  URL https://github.com/samtools/htscodecs/releases/download/v1.6.7/htscodecs-1.6.7.tar.gz
  URL_HASH SHA256=5100e4b27646a27042a00b07cada35ce9fa3c7fc5aacd0de9874cd5d16c35fbf)
FetchContent_GetProperties(htscodecs)
if(NOT htscodecs_POPULATED)
  FetchContent_Populate(htscodecs)
endif()

# The sources include "config.h", which the codec's two files need nothing from.
file(WRITE ${CMAKE_BINARY_DIR}/htscodecs-config/config.h "")

find_package(Threads REQUIRED)
add_library(ramtools_fqzcomp STATIC
  ${htscodecs_SOURCE_DIR}/htscodecs/fqzcomp_qual.c
  ${htscodecs_SOURCE_DIR}/htscodecs/utils.c)
set_target_properties(ramtools_fqzcomp PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 99)
target_include_directories(ramtools_fqzcomp
  PRIVATE ${CMAKE_BINARY_DIR}/htscodecs-config
  PUBLIC $<BUILD_INTERFACE:${htscodecs_SOURCE_DIR}>)
target_link_libraries(ramtools_fqzcomp PRIVATE Threads::Threads m)
add_library(ramtools::fqzcomp ALIAS ramtools_fqzcomp)
message(STATUS "fqzcomp: htscodecs 1.6.7 built from source")
