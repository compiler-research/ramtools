set(_gtest_byproduct_binary_dir
  ${CMAKE_BINARY_DIR}/test/googletest-prefix/src/googletest-build)
set(_gtest_byproducts
  ${_gtest_byproduct_binary_dir}/lib/libgtest.a
  ${_gtest_byproduct_binary_dir}/lib/libgtest_main.a
  ${_gtest_byproduct_binary_dir}/lib/libgmock.a
  ${_gtest_byproduct_binary_dir}/lib/libgmock_main.a
  )

if(MSVC)
  set(EXTRA_GTEST_OPTS
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${_gtest_byproduct_binary_dir}/lib/
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL:PATH=${_gtest_byproduct_binary_dir}/lib/
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${_gtest_byproduct_binary_dir}/lib/
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO:PATH=${_gtest_byproduct_binary_dir}/lib/
    -Dgtest_force_shared_crt=ON)
elseif(APPLE)
  set(EXTRA_GTEST_OPTS -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
                       -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES})
endif()

include(ExternalProject)
ExternalProject_Add(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  EXCLUDE_FROM_ALL 1
  GIT_SHALLOW 1
  GIT_TAG v1.14.0
  UPDATE_COMMAND ""
  CMAKE_ARGS -G ${CMAKE_GENERATOR}
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                -DCMAKE_AR=${CMAKE_AR}
                -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                ${EXTRA_GTEST_OPTS}
  BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config Release
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS ${_gtest_byproducts}
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON
  TIMEOUT 600
  )

ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIR ${source_dir}/googletest/include)
set(GMOCK_INCLUDE_DIR ${source_dir}/googlemock/include)
file(MAKE_DIRECTORY ${GTEST_INCLUDE_DIR} ${GMOCK_INCLUDE_DIR})

ExternalProject_Get_Property(googletest binary_dir)
set(_G_LIBRARY_PATH ${binary_dir}/lib/)

foreach(lib gtest gtest_main gmock gmock_main)
  add_library(${lib} IMPORTED STATIC GLOBAL)
  set_target_properties(${lib} PROPERTIES
    IMPORTED_LOCATION "${_G_LIBRARY_PATH}${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}"
    INTERFACE_INCLUDE_DIRECTORIES "${GTEST_INCLUDE_DIR}"
    )
  add_dependencies(${lib} googletest)
endforeach()

target_include_directories(gtest INTERFACE ${GTEST_INCLUDE_DIR})
target_include_directories(gmock INTERFACE ${GMOCK_INCLUDE_DIR})

add_library(GTest::gtest ALIAS gtest)
add_library(GTest::gtest_main ALIAS gtest_main)
add_library(GTest::gmock ALIAS gmock)
add_library(GTest::gmock_main ALIAS gmock_main)

