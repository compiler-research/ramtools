#include "ramcore/RAMNTupleView.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RNTupleDS.hxx>
#include <algorithm>

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <Rtypes.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TString.h>

#include "rntuple/RAMNTupleRecord.h"

namespace {

constexpr int FLAG_UNMAPPED = 0x4;
constexpr int FLAG_SECONDARY = 0x100;
constexpr int FLAG_SUPPLEMENTARY = 0x800;
constexpr int FLAG_FILTER = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY;

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
   return span;
}
int GetRefId(ROOT::RDataFrame &df, const std::string &rname)
{
   if (rname == "*")
      return -1;

   auto refs = df.Take<std::vector<std::string>>("rname_refs");
   const auto &refids = refs.GetValue()[0];

   auto it = std::find(refids.begin(), refids.end(), rname);
   return (it == refids.end()) ? -1 : std::distance(refids.begin(), it);
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

   // No colon means chromosome-only query e.g. "chr1"
   const std::size_t colonPos = region.find(':');
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
Long64_t ramntupleview(const char *file, const char *query, bool /*cache*/, bool /*perfstats*/,
                       const char * /*perfstatsfilename*/)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      std::cerr << "ramntupleview: failed to open file " << file << std::endl;
      return 0;
   }

   std::string region = query ? query : "";
   if (region.empty() || region == "*") {
      stopwatch.Print();
      return reader->GetNEntries();
   }

   TString rname;
   Int_t rs, re;
   if (!parseRegion(region, rname, rs, re)) {
      std::cerr << "Invalid region format. Use rname[:start[-end]]\n";
      return 0;
   }

   int refid = resolveRefId(rname.Data());
   if (refid < 0) {
      std::cerr << "Reference '" << rname.Data() << "' not found\n";
      return 0;
   }

   auto flagView = reader->GetView<uint16_t>("record.flag");
   auto refidView = reader->GetView<int32_t>("record.refid");
   auto posView = reader->GetView<int32_t>("record.pos");
   auto cigarView = reader->GetView<std::vector<uint32_t>>("record.cigar");

   auto index = RAMNTupleRecord::GetIndex();
   Long64_t start = (index && index->Size() > 0) ? index->GetRow(refid, rs) : 0;
   if (start < 0)
      start = 0;

   Long64_t count = 0;
   const Long64_t total = reader->GetNEntries();

   for (Long64_t i = start; i < total; i++) {
      if (flagView(i) & FLAG_FILTER)
         continue;
      int curRef = refidView(i);
      if (curRef < refid)
         continue;
      if (curRef > refid)
         break;

      int pos = posView(i);
      if (pos > re)
         break;

      // Overlap: read_end >= rs (pos <= re guaranteed by break above)
      if (pos >= rs) {
         count++;
      } else {
         int readEnd = pos + computeRefSpan(cigarView(i)) - 1;
         if (readEnd >= rs)
            count++;
      }
   }

   stopwatch.Print();
   std::cout << "Found " << count << " records in region " << region << std::endl;
   return count;
}
// NOLINTNEXTLINE(misc-use-internal-linkage)
ULong64_t mt_ramntupleview(const int numthreads, const char *file, const char *query, bool /*cache*/,
                           bool /*perfstats*/, const char * /*perfstatsfilename*/)
{
   TStopwatch st;
   st.Start();
   TString rname;
   std::string region = query;
   Int_t start = 0;
   Int_t end = 0;
   if (!parseRegion(region, rname, start, end)) {
      std::cerr << "Invalid region format. Use rname[:start[-end]]\n";
      return 0;
   }
   auto metadata = ROOT::RDF::FromRNTuple("METADATA", file);
   const int refid = GetRefId(metadata, rname.Data());
   if (refid < 0) {
      std::cerr << "Reference" << rname.Data() << " not found\n";
   }
   ROOT::EnableImplicitMT(numthreads);
   auto ram = ROOT::RDF::FromRNTuple("RAM", file);
   auto filterfunc = [refid, start, end](int32_t refidentry, int32_t pos) {
      return (refid == refidentry) && (pos >= start) && (pos <= end);
   };

   auto filtered = ram.Filter(filterfunc, {"record.refid", "record.pos"});
   auto count = filtered.Count();
   *count;
   st.Print();
   return *count;
}
