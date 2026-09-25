#ifndef RAMCORE_SAMTONTUPLE_H
#define RAMCORE_SAMTONTUPLE_H

#include <cstddef>
#include <cstdint>
#include <string>

/// Converts a SAM file to one RAM file.
///
/// - `datafile`: input SAM
/// - `treefile`: output RAM file, created or overwritten
/// - `compression_algorithm`: ROOT compression code, algorithm * 100 + level (505 is ZSTD level 5)
/// - `quality_policy`: one of RAMNTupleRecord::EQualCompressionBits
/// - `threads`: threads that parse, encode and compress, each a block of its own
/// - `block_bytes`: size of the blocks of whole lines the input is read in
///
/// The records keep their input order whatever the thread count. Malformed
/// records are reported on standard error and skipped. Returns false when the
/// input cannot be opened or read or the output cannot be created.
bool samtoramntuple(const char *datafile, const char *treefile, int compression_algorithm, uint32_t quality_policy,
                    int threads = 1, size_t block_bytes = 64U << 20);

/// Writes one RAM file per reference, named `<output_prefix>_<rname>.root`.
/// Records stream to their file as they are parsed, in the order they arrive.
void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy);

#endif
