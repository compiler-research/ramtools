#include <ROOT/RDataFrame.hxx>
#include <ramcore/RDF_RAMNTupleView.h>
using namespace ROOT::RDF;
namespace rdf {

static int GetRefId(ROOT::RDataFrame &df, const std::string &rname)
{
   if (rname == "*")
      return -1;

   auto refs = df.Take<std::vector<std::string>>("rname_refs");
   const auto &vec = refs.GetValue()[0];

   auto it = std::find(vec.begin(), vec.end(), rname);
   return (it == vec.end()) ? -1 : std::distance(vec.begin(), it);
}
} // namespace rdf

Long64_t rdf_ramntupleview(const int num_threads, const char *file, const char *query, bool cache, bool perfstats,
                           const char *perfstatsfilename)
{
   TStopwatch ts;
   ts.Start();
   std::string region = query;
   int chrDelimiterPos = region.find(":");
   if (chrDelimiterPos == std::string::npos) {
      std::cerr << "Invalid region format. Use rname:start-end\n";
      return 0;
   }
   const TString rname = region.substr(0, chrDelimiterPos);
   int rangeDelimiterPos = region.find("-", chrDelimiterPos);
   if (rangeDelimiterPos == std::string::npos) {
      std::cerr << "Invalid region format. Use rname:start-end\n";
      return 0;
   }
   const Int_t range_start = std::stoi(region.substr(chrDelimiterPos + 1, rangeDelimiterPos - chrDelimiterPos - 1));
   const Int_t range_end = std::stoi(region.substr(rangeDelimiterPos + 1));

   //  auto index = FromRNTuple("INDEX", file);
   auto metadata = FromRNTuple("METADATA", file);
   const auto ref = rdf::GetRefId(metadata, rname.Data());
   if (ref == -1) {
      std::cerr << "Chromosome not found:" << rname << "\n";
   }
   ROOT::EnableImplicitMT(num_threads);
   auto ram = FromRNTuple("RAM", file);
   auto check = [ref, range_start, range_end](int32_t refid, int32_t pos) {
      return (refid == ref) && (pos >= range_start - 1) && (pos <= range_end - 1);
   };
   auto filtered = ram.Filter(check, {"record.refid", "record.pos"});
   auto num = filtered.Count();
   *num;
   ts.Print();
   return *num;
}
