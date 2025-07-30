ROOT_EXECUTABLE(samtoramntuple samtoramntuple.cxx
    LIBRARIES ramtools_ntuple ROOT::Core ROOT::RIO ROOT::ROOTNTuple ROOT::ROOTNTupleUtil
)

configure_file(${CMAKE_SOURCE_DIR}/data/samexample.sam 
               ${CMAKE_CURRENT_BINARY_DIR}/samexample.sam COPYONLY)

set(ROOT_MACROS
    checkindex.cxx
    rammerge.cxx
    ramreader.cxx
    ramrandom.cxx
    ramview.cxx
    ramview_no_index.cxx
    ramntupleview.cxx
    samtoram.cxx
    parsetreestats.cxx
)

install(FILES ${ROOT_MACROS}
    DESTINATION share/ramtools/macros
)