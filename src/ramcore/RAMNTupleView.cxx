#include "ramcore/RAMNTupleView.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <TStopwatch.h>
#include <TString.h>
#include "rntuple/RAMNTupleRecord.h"
#include <Rtypes.h>

Long64_t ramntupleview(const char *file, const char *query, bool cache, bool perfstats, const char *perfstatsfilename)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      printf("ramntupleview: failed to open file %s\n", file);
      return 0;
   }

   std::string region = query;
   if (region.empty() || region == "*") {
      stopwatch.Print();
      return reader->GetNEntries();
   }

   TString rname;
   Int_t range_start = 1;
   Int_t range_end = std::numeric_limits<Int_t>::max();

   int chrDelimiterPos = region.find(":");
   if (chrDelimiterPos == std::string::npos) {
      rname = region;
   } else {
      rname = region.substr(0, chrDelimiterPos);
      int rangeDelimiterPos = region.find("-", chrDelimiterPos);
      if (rangeDelimiterPos == std::string::npos) {
         try {
            range_start = std::stoi(region.substr(chrDelimiterPos + 1));
            range_end = range_start;
         } catch (...) {
            std::cerr << "Invalid region format. Use rname:start-end or rname:position\n";
            return 0;
         }
      } else {
         range_start = std::stoi(region.substr(chrDelimiterPos + 1, rangeDelimiterPos - chrDelimiterPos - 1));
         range_end = std::stoi(region.substr(rangeDelimiterPos + 1));
      }
   }

   auto refs = RAMNTupleRecord::GetRnameRefs();
   auto refid = refs->GetRefId(rname.Data());

   if (refid < 0) {
      if (rname.BeginsWith("chr")) {
         TString stripped_rname = rname(3, rname.Length() - 3);
         refid = refs->GetRefId(stripped_rname.Data());
      }
      if (refid < 0 && (rname == "chrM" || rname == "M")) {
         refid = refs->GetRefId("MT");
      }
   }

   if (refid < 0) {
      std::cerr << "Error: Reference name '" << rname.Data() << "' not found in file.\n";
      return 0;
   }

   auto index = RAMNTupleRecord::GetIndex();
   auto flagView = reader->GetView<uint16_t>("record.flag");
   auto refidView = reader->GetView<int32_t>("record.refid");
   auto posView = reader->GetView<int32_t>("record.pos");
   auto cigarView = reader->GetView<std::vector<uint32_t>>("record.cigar");

   Long64_t count = 0;
   const Long64_t totalEntries = reader->GetNEntries();

   const int FLAG_FILTER = 0x904;
   const Int_t rs0 = (range_start > 0) ? (range_start - 1) : 0;
   const Int_t re0 = (range_end > 0) ? (range_end - 1) : 0;
   constexpr int kMaxRefSpanHeuristic = 10000;

   auto computeRefSpan = [](const std::vector<uint32_t> &cigarOps) -> int {
      int span = 0;
      for (uint32_t op : cigarOps) {
         int len = static_cast<int>(op >> 4);
         int code = static_cast<int>(op & 0xF);

         if (code == 0 || code == 2 || code == 3 || code == 7 || code == 8) {
            span += len;
         }
      }
      return span;
   };

   if (!index || index->Size() == 0) {

      bool seenRef = false;
      for (auto i : reader->GetEntryRange()) {
         const auto flag = flagView(i);
         if (flag & FLAG_FILTER)
            continue;

         const auto curRef = refidView(i);
         if (curRef == refid) {
            seenRef = true;
         } else {
            if (seenRef && curRef > refid)
               break;
            continue;
         }
         const auto curPos = posView(i);
         if (curPos > re0)
            break;

         int read_start = curPos;
         if (read_start >= rs0) {
            count++;
            continue;
         }
         if (read_start + kMaxRefSpanHeuristic < rs0) {
            continue;
         }
         const auto &cigarOps = cigarView(i);
         int refSpan = computeRefSpan(cigarOps);
         int read_end = (refSpan > 0) ? (read_start + refSpan - 1) : read_start;

         if (read_start <= re0 && read_end >= rs0) {
            count++;
         }
      }
   } else {

      auto start_entry = index->GetRow(refid, rs0);

      if (start_entry < 0)
         start_entry = index->GetRow(refid, 0);
      if (start_entry < 0)
         start_entry = 0;

      for (Long64_t j = start_entry; j < totalEntries; j++) {
         const auto flag = flagView(j);
         if (flag & FLAG_FILTER)
            continue;

         const auto curRef = refidView(j);
         if (curRef < refid)
            continue;
         if (curRef > refid)
            break;

         const auto curPos = posView(j);
         if (curPos > re0)
            break;

         int read_start = curPos;
         if (read_start >= rs0) {
            count++;
            continue;
         }
         if (read_start + kMaxRefSpanHeuristic < rs0) {
            continue;
         }
         const auto &cigarOps = cigarView(j);
         int refSpan = computeRefSpan(cigarOps);
         int read_end = (refSpan > 0) ? (read_start + refSpan - 1) : read_start;

         if (read_end >= rs0) {
            count++;
         }
      }
   }

   stopwatch.Print();
   std::cout << "Found " << static_cast<long long>(count) << " records in region " << (query ? query : "") << std::endl;
   return count;
}
