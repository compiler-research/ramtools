#ifndef RAMCORE_RAMNTUPLEVIEW_H
#define RAMCORE_RAMNTUPLEVIEW_H

void ramntupleview(const char *file, const char *query = "", bool cache = true, bool perfstats = false,
                   const char *perfstatsfilename = "perf.root");

#endif // RAMCORE_RAMNTUPLEVIEW_H
