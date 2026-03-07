#ifndef RAMCORE_RDF_RAMNTUPLEVIEW_H
#define RAMCORE_RDF_RAMNTUPLEVIEW_H
#include <RtypesCore.h>

ULong64_t rdf_ramntupleview(const char *file, const char *query = "", bool cache = true, bool perfstats = false,
                       const char *perfstatsfilename = "perf.root");

#endif  //RAMCORE_RDF_RAMNTUPLEVIEW_H

