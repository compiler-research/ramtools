#ifndef RAMCORE_SAMTONTUPLE_H
#define RAMCORE_SAMTONTUPLE_H

#include <cstdint> 

void samtoramntuple(const char *datafile,
                    const char *treefile,
                    bool index, bool split, bool cache,
                    int compression_algorithm,
                    uint32_t quality_policy);

#endif // RAMCORE_SAMTONTUPLE_H 

