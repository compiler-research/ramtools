# ramtools::libdeflate: the system libdeflate if pkg-config finds it (htslib is
# usually built with it), otherwise a pinned release built from source.

pkg_check_modules(libdeflate QUIET IMPORTED_TARGET libdeflate)

if(libdeflate_FOUND)
  add_library(ramtools_libdeflate INTERFACE)
  target_link_libraries(ramtools_libdeflate INTERFACE PkgConfig::libdeflate)
  add_library(ramtools::libdeflate ALIAS ramtools_libdeflate)
  message(STATUS "libdeflate: system libdeflate ${libdeflate_VERSION}")
  return()
endif()

include(FetchContent)
FetchContent_Declare(libdeflate
  URL https://github.com/ebiggers/libdeflate/releases/download/v1.26/libdeflate-1.26.tar.gz
  URL_HASH SHA256=125856d4656e0feab660f94842f835923410c9281fedbcee64598c918da42b5a)
set(LIBDEFLATE_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
set(LIBDEFLATE_BUILD_GZIP OFF CACHE BOOL "" FORCE)
set(LIBDEFLATE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(libdeflate)
set_target_properties(libdeflate_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(ramtools::libdeflate ALIAS libdeflate_static)
message(STATUS "libdeflate: 1.26 built from source")
