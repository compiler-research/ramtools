#ifndef RAMCORE_RAMNTUPLEVIEW_H
#define RAMCORE_RAMNTUPLEVIEW_H
#include <Rtypes.h>
#include <functional>
#include <string>

class RAMFile;

struct RAMNTupleViewOpts {
   bool fCache = true;
   bool fPerfStats = false;
   std::string perfStatsFilename = "perf.root";
};

/// Calls on_row for every record overlapping query, in file order, and returns
/// how many there were. The overlap rule is samtools': same reference, and
/// [POS, POS+refspan-1] meets the query interval. An empty query or "*" visits
/// the whole file. on_row may be null.
///
/// On a coordinate-sorted file the scan starts at the first record at or after
/// the query start minus the longest span in the file, found by binary search
/// over the refid and pos columns, and stops at the first record past the end.
/// An unsorted file is read end to end with the same overlap rule.
Long64_t ramntuplescan(RAMFile &file, const char *query, const std::function<void(Long64_t)> &on_row);

/// Counts the records overlapping query.
Long64_t ramntupleview(const char *file, const char *query = "", const RAMNTupleViewOpts & = RAMNTupleViewOpts());

#endif // RAMCORE_RAMNTUPLEVIEW_H
