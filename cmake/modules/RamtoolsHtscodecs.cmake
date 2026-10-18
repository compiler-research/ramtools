# ramtools::htscodecs: the system htscodecs if pkg-config finds it, otherwise the
# source files of fqzcomp and the name tokeniser, built from a pinned release.

pkg_check_modules(htscodecs QUIET IMPORTED_TARGET htscodecs)

if(htscodecs_FOUND)
  add_library(ramtools_htscodecs INTERFACE)
  target_link_libraries(ramtools_htscodecs INTERFACE PkgConfig::htscodecs)
  add_library(ramtools::htscodecs ALIAS ramtools_htscodecs)
  message(STATUS "htscodecs: system htscodecs ${htscodecs_VERSION}")
  return()
endif()

include(FetchContent)
FetchContent_Declare(htscodecs
  URL https://github.com/samtools/htscodecs/releases/download/v1.6.7/htscodecs-1.6.7.tar.gz
  URL_HASH SHA256=5100e4b27646a27042a00b07cada35ce9fa3c7fc5aacd0de9874cd5d16c35fbf)
# htscodecs builds with autotools and has no CMakeLists.txt, so this only
# downloads and unpacks it.
FetchContent_MakeAvailable(htscodecs)

# The sources include "config.h". Left empty, it builds the portable rANS coder
# without SIMD, and without bzip2, which the name tokeniser uses only with
# arithmetic coding at its high levels.
file(WRITE ${CMAKE_BINARY_DIR}/htscodecs-config/config.h "")

find_package(Threads REQUIRED)
add_library(ramtools_htscodecs STATIC
  ${htscodecs_SOURCE_DIR}/htscodecs/fqzcomp_qual.c
  ${htscodecs_SOURCE_DIR}/htscodecs/tokenise_name3.c
  ${htscodecs_SOURCE_DIR}/htscodecs/arith_dynamic.c
  ${htscodecs_SOURCE_DIR}/htscodecs/rANS_static4x16pr.c
  ${htscodecs_SOURCE_DIR}/htscodecs/rANS_static32x16pr.c
  ${htscodecs_SOURCE_DIR}/htscodecs/pack.c
  ${htscodecs_SOURCE_DIR}/htscodecs/rle.c
  ${htscodecs_SOURCE_DIR}/htscodecs/utils.c)
set_target_properties(ramtools_htscodecs PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 99)
target_include_directories(ramtools_htscodecs
  PRIVATE ${CMAKE_BINARY_DIR}/htscodecs-config
  PUBLIC $<BUILD_INTERFACE:${htscodecs_SOURCE_DIR}>)
target_link_libraries(ramtools_htscodecs PRIVATE Threads::Threads m)
add_library(ramtools::htscodecs ALIAS ramtools_htscodecs)
message(STATUS "htscodecs: 1.6.7 built from source")
