find_program(SPHINX_BUILD sphinx-build)
if(NOT SPHINX_BUILD)
    message(FATAL_ERROR "RAMTOOLS_BUILD_DOCS is ON but sphinx-build was not found")
endif()

add_custom_target(docs
    COMMAND ${SPHINX_BUILD} -q -W -b html ${CMAKE_CURRENT_SOURCE_DIR}/source ${CMAKE_BINARY_DIR}/docs/html
    COMMENT "Building the documentation into ${CMAKE_BINARY_DIR}/docs/html"
)
