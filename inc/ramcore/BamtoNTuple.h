#ifndef RAMCORE_BAMTONTUPLE_H
#define RAMCORE_BAMTONTUPLE_H

#include <cstdint>

void bamtoramntuple(const char *bamfile, const char *treefile, bool index, bool split, bool cache,
                    int compression_algorithm, uint32_t quality_policy);

#endif // RAMCORE_BAMTONTUPLE_H