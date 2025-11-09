#ifndef RAMCORE_SAMTONTUPLE_H
#define RAMCORE_SAMTONTUPLE_H

#include <cstdint>
#include <string>
#include <vector>

void samtoramntuple(const char *datafile,
                    const char *treefile,
                    bool index, bool split, bool cache,
                    int compression_algorithm,
                    uint32_t quality_policy);

void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy, int num_threads = 4,
                                        const std::vector<std::string> &only_chromosomes = {});

#endif
