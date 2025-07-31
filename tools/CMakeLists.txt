add_executable(samtoram samtoram.cxx)
target_link_libraries(samtoram PRIVATE
    RAMCore
    ROOT::Core
    ROOT::Tree
)

add_executable(samtoramntuple samtoramntuple.cxx)
target_link_libraries(samtoramntuple PRIVATE
    RAMCore
    ROOT::Core
    ROOT::ROOTNTuple
)

install(TARGETS samtoram samtoramntuple
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

set(ROOT_MACROS
    checkindex.cxx
    rammerge.cxx
    ramreader.cxx
    ramrandom.cxx
    ramview.cxx
    ramview_no_index.cxx
    ramntupleview.cxx
    parsetreestats.cxx
)

install(FILES ${ROOT_MACROS}
    DESTINATION share/ramtools/macros
)

