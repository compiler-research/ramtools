ROOT_EXECUTABLE(samtoram
    samtoram.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::Tree
)

ROOT_EXECUTABLE(samtoramntuple
    samtoramntuple.cxx
    LIBRARIES
        ramcore
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

