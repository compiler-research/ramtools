#ifndef RAMCORE_BAMTONTUPLE_H
#define RAMCORE_BAMTONTUPLE_H

#include <cstddef>
#include <cstdint>

/// Converts a BAM file to one RAM file on \p threads threads, each filling
/// blocks of about \p block_bytes of records; the records keep their input
/// order. htslib inflates the input on as many threads again.
void bamtoramntuple(const char *bamfile, const char *treefile, bool split, bool cache, int compression_algorithm,
                    uint32_t quality_policy, int threads = 1, size_t block_bytes = 64U << 20);

#endif // RAMCORE_BAMTONTUPLE_H
