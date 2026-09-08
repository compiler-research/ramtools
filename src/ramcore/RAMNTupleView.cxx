#include "ramcore/RAMNTupleView.h"
#include <algorithm>

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <Rtypes.h>
#include <TStopwatch.h>
#include <TString.h>

#include "rntuple/RAMNTupleRecord.h"

namespace {

// Computes how many reference bases an alignment covers using CIGAR.
//
// CIGAR (Compact Idiosyncratic Gapped Alignment Report) describes how a read
// aligns to the reference. It's a sequence of operations like "50M2D48M" meaning:
// 50 bases match, 2 bases deleted from reference, 48 bases match.
//
// BAM stores only alignment start position. To find the end position (needed
// for overlap queries), we compute: end = start + refSpan - 1
//
// In BAM, each CIGAR op is packed into a uint32: (length << 4) | opcode
// Only some ops consume reference bases:
//   M(0)=match/mismatch, D(2)=deletion, N(3)=skip, =(7)=seq match, X(8)=seq mismatch
// Ops like I(insertion), S(soft-clip), H(hard-clip) don't consume reference.
//
// Example: CIGAR "50M2D48M" -> refSpan = 50 + 2 + 48 = 100
// (the read is 98bp but covers 100 reference bases due to 2bp deletion)
//
// See SAM spec section 1.4.6: https://samtools.github.io/hts-specs/SAMv1.pdf
int computeRefSpan(const std::vector<uint32_t> &cigarOps)
{
   int span = 0;
   for (uint32_t op : cigarOps) {
      int code = op & 0xF;
      if (code == 0 || code == 2 || code == 3 || code == 7 || code == 8) {
         span += static_cast<int>(op >> 4);
      }
   }
   // A placed record with no CIGAR still occupies the base it sits on; htslib's
   // bam_endpos says the same, so region queries agree with samtools.
   return span > 0 ? span : 1;
}

int resolveRefId(const char *name)
{
   auto refs = RAMNTupleRecord::GetRnameRefs();
   if (!refs)
      return -1;

   const auto &refVec = refs->GetRefs();
   for (size_t i = 0; i < refVec.size(); i++) {
      if (refVec[i] == name)
         return static_cast<int>(i);
   }
   return -1;
}

bool parseRegion(const std::string &region, TString &rname, Int_t &start, Int_t &end)
{
   // Default: entire chromosome
   start = 0;
   end = std::numeric_limits<Int_t>::max() - 1;

   // No colon means chromosome-only query e.g. "chr1". Split on the last colon:
   // reference names may contain colons themselves (GRCh38 has "HLA-A*01:01:01:01").
   const std::size_t colonPos = region.rfind(':');
   if (colonPos == std::string::npos) {
      rname = region;
      return true;
   }

   // Split "chr1:1000-2000" into rname="chr1", rest="1000-2000"
   rname = region.substr(0, colonPos);
   const std::string rest = region.substr(colonPos + 1);

   const std::size_t dashPos = rest.find('-');
   if (dashPos == std::string::npos) {
      // Single position e.g. "chr1:500"
      if (rest.empty() || !std::all_of(rest.begin(), rest.end(), ::isdigit))
         return false;
      start = end = std::stoi(rest) - 1; // SAM 1-based → 0-based
   } else {
      // Range e.g. "chr1:1000-2000"
      const std::string startStr = rest.substr(0, dashPos);
      const std::string endStr = rest.substr(dashPos + 1);

      // Reject empty or non-numeric strings before calling stoi
      if (startStr.empty() || endStr.empty() || !std::all_of(startStr.begin(), startStr.end(), ::isdigit) ||
          !std::all_of(endStr.begin(), endStr.end(), ::isdigit))
         return false;

      start = std::stoi(startStr) - 1; // SAM 1-based → 0-based
      end = std::stoi(endStr) - 1;
   }

   // Clamp negative start (e.g. user passed position 0)
   if (start < 0)
      start = 0;
   return true;
}

} // namespace

// NOLINTNEXTLINE(misc-use-internal-linkage)
Long64_t ramntuplescan(ROOT::RNTupleReader &reader, const char *query, const std::function<void(Long64_t)> &on_row)
{
   const std::string region = query ? query : "";
   const auto total = static_cast<Long64_t>(reader.GetNEntries());

   // No region means every record, the way `samtools view` with no region does.
   if (region.empty() || region == "*") {
      if (on_row) {
         for (Long64_t i = 0; i < total; i++)
            on_row(i);
      }
      return total;
   }

   TString rname;
   Int_t rs = 0;
   Int_t re = std::numeric_limits<Int_t>::max() - 1;

   // Try the whole string as a reference name before reading its last colon as a
   // coordinate separator, the same order samtools uses.
   int refid = resolveRefId(region.c_str());
   if (refid < 0) {
      if (!parseRegion(region, rname, rs, re)) {
         std::cerr << "Invalid region format. Use rname[:start[-end]]\n";
         return 0;
      }
      refid = resolveRefId(rname.Data());
      if (refid < 0) {
         std::cerr << "Reference '" << rname.Data() << "' not found\n";
         return 0;
      }
   }

   auto refidView = reader.GetView<int32_t>("record.refid");
   auto posView = reader.GetView<int32_t>("record.pos");
   auto cigarView = reader.GetView<std::vector<uint32_t>>("record.cigar");

   // The index entry at or before the region start is not enough on its own: a
   // read beginning earlier can still reach into the region, and starting there
   // would step over it. Backing off by the longest span in the file is exact.
   // A file that does not record the span (0) is scanned from the reference's
   // first entry instead.
   // Seeking and stopping early both assume coordinate order. A file without
   // it is read end to end with the same overlap test.
   const bool sorted = RAMNTupleRecord::IsCoordinateSorted();
   auto index = RAMNTupleRecord::GetIndex();
   Long64_t start = 0;
   if (sorted && index && index->Size() > 0) {
      const Int_t maxSpan = static_cast<Int_t>(RAMNTupleRecord::GetMaxRefSpan());
      const Int_t seekPos = (maxSpan > 0 && rs > maxSpan) ? rs - maxSpan : 0;
      start = index->GetRow(refid, seekPos);
      if (start < 0)
         start = 0;
   }

   Long64_t count = 0;

   for (Long64_t i = start; i < total; i++) {
      const int curRef = refidView(i);
      if (curRef != refid) {
         if (sorted && curRef > refid)
            break;
         continue;
      }

      const int pos = posView(i);
      if (pos > re) {
         if (sorted)
            break;
         continue;
      }

      if (pos < rs && pos + computeRefSpan(cigarView(i)) - 1 < rs)
         continue;

      count++;
      if (on_row)
         on_row(i);
   }

   return count;
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
Long64_t ramntupleview(const char *file, const char *query, const RAMNTupleViewOpts &)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      std::cerr << "ramntupleview: failed to open file " << file << std::endl;
      return 0;
   }

   const Long64_t count = ramntuplescan(*reader, query, nullptr);

   stopwatch.Print();
   std::cout << "Found " << count << " records in region " << (query ? query : "") << std::endl;
   return count;
}
