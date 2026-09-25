#ifndef RAMCORE_BAMTONTUPLE_H
#define RAMCORE_BAMTONTUPLE_H

#include <cstdint>

/// Converts a BAM file to one RAM file through htslib, producing the same
/// file samtoramntuple() would from the equivalent SAM. \a split and \a cache
/// are accepted but not used; the other parameters are as for samtoramntuple().
void bamtoramntuple(const char *bamfile, const char *treefile, bool split, bool cache, int compression_algorithm,
                    uint32_t quality_policy);

#endif // RAMCORE_BAMTONTUPLE_H