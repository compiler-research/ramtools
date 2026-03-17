#include <ROOT/RDataFrame.hxx>
#include <ramcore/RDF_RAMNTupleView.h>
using namespace ROOT::RDF;
namespace rdf{

static int GetRefId(ROOT::RDataFrame& rname_refs, const std::string& rname){
   if (rname == "*"){
         return -1;
   }
   int result = -1;
   auto find_index = [rname, &result](ROOT::VecOps::RVec<std::string> refs){
      for(int i = 0; i < refs.size(); i++){
         if(refs[i] == rname){
           result = i;
         }
      }
   };
   rname_refs.Foreach(find_index, {"rname_refs"});
   return result;
}

static void Display(ROOT::RDataFrame ram, ROOT::RDataFrame index, ROOT::RDataFrame metadata){
   
   auto inde = index.Describe();
   inde.Print();
   auto description = ram.Describe();
   description.Print();
   auto des = metadata.Describe();
   des.Print();
   auto display = ram.Display({"record.qname", "record.refid", "record.refnext", "record.seq", "record.pos", "record.cigar"});
   display->Print();  
   auto dis = metadata.Display({"rname_refs", "rnext_refs"});
   dis->Print();
   auto d = index.Display({"index_entries.entry", "index_entries.pos", "index_entries.refid"});
   d->Print();
}


} // namespace rdf


ULong64_t
rdf_ramntupleview(const char *file, const char *query, bool cache, bool perfstats, const char *perfstatsfilename)
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
   ROOT::EnableImplicitMT();
   auto ram = FromRNTuple("RAM", file);
   
   auto check = [&range_start, &range_end, &ref](int32_t refid, int32_t pos){
      return refid == ref && pos >= range_start - 1 && pos <= range_end - 1;
   };
   auto filtered = ram.Filter(check, {"record.refid", "record.pos"});
   auto num = filtered.Count();
   *num;
   ts.Print();
   return *num; 
}
