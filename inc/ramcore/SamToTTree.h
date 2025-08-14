#ifndef RAMCORE_SAMTOTREE_H
#define RAMCORE_SAMTOTREE_H

#include <Rtypes.h> 

void samtoram(const char *datafile,
              const char *treefile,
              bool index, bool split, bool cache,
              Int_t compression_algorithm,
              UInt_t quality_policy);

#endif // RAMCORE_SAMTOTREE_H 

