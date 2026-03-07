

#include <iostream>
#include <ramcore/RDF_RAMNTupleview.h>
int main(int argc, char **argv){

if (argc < 2){
      std::cerr << "Usage: " << argv[0] << " <file.root> [rname:start-end]\n";
      std::cerr << "Example: " << argv[0] << " output.root chr1:1000-2000\n";
      return 1;
}

const char* file = argv[1];

const char* query = argv[2];
ULong64_t reads = rdf_ramntupleview(file, query);

}
