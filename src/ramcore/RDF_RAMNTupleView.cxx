#include <ROOT/RDataFrame.hxx>
#include <ramcore/RDF_RAMNTupleView.h>

ULong64_t
rdf_ramntupleview(const char *file, const char *query, bool cache, bool perfstats, const char *perfstatsfilename)
{

   TStopwatch ts;
   ts.Start();
   auto reader = ROOT::RDF::FromRNTuple("RAM", file);
   auto description = reader.Describe();
   description.Print();
   ts.Print();
   return 0;
}
