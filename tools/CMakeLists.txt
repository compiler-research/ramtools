ROOT_EXECUTABLE(samtoramntuple
    samtoramntuple.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::ROOTNTuple
)

ROOT_EXECUTABLE(ramntupleview
    ramntupleview.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::ROOTNTuple
)

ROOT_EXECUTABLE(ramdump
    ramdump.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::RIO
        ROOT::ROOTNTuple
)

ROOT_EXECUTABLE(bamtoramntuple
    bamtoramntuple.cxx
    LIBRARIES
        ramcore
        PkgConfig::htslib
        ROOT::Core
        ROOT::ROOTNTuple
)

install(TARGETS samtoramntuple ramntupleview ramdump bamtoramntuple
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
