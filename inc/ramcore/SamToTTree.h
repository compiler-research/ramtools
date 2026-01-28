#ifndef RAMCORE_SAMTOTREE_H
#define RAMCORE_SAMTOTREE_H

void samtoram(const char *datafile, const char *treefile, bool index, bool split, bool cache, int compression_algorithm,
              unsigned int quality_policy);

#endif // RAMCORE_SAMTOTREE_H
