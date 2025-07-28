function(add_ram_tool tool_name)
    add_executable(${tool_name} ${tool_name}.cxx)
    
    if(${tool_name} MATCHES "ntuple")
        target_link_libraries(${tool_name} PRIVATE
            ramtools_ntuple
            ROOT::Core
            ROOT::RIO
            ROOT::ROOTNTuple
            ROOT::ROOTNTupleUtil
        )
    else()
        target_link_libraries(${tool_name} PRIVATE
            ramtools_tree
            ROOT::Core
            ROOT::RIO
            ROOT::Tree
        )
    endif()
    
    install(TARGETS ${tool_name} RUNTIME DESTINATION bin)
endfunction()

add_ram_tool(samtoramntuple)

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

