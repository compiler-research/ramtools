#ifndef RAMCORE_RAMNTUPLEVIEW_H
#define RAMCORE_RAMNTUPLEVIEW_H
#include <Rtypes.h>

Long64_t ramntupleview(const char *file, const char *query = "", bool cache = true, bool perfstats = false,
                       const char *perfstatsfilename = "perf.root");
ULong64_t mt_ramntupleview(int numthreads, const char *file, const char *query = "", bool cache = true,
                           bool perfstats = false, const char *perfstatsfilename = "perf.root");
#endif // RAMCORE_RAMNTUPLEVIEW_H
