#ifndef RAMCORE_RAMSORT_H
#define RAMCORE_RAMSORT_H

/// Sort a RAM (RNTuple) file by coordinate (refid, pos) or by QNAME.
/// \param inputFile   Path to input .root RAM file
/// \param outputFile  Path to output .root RAM file
/// \param byName      If true, sort by QNAME; otherwise sort by (refid, pos)
/// \return 0 on success, 1 on error
int ramsortntuple(const char *inputFile, const char *outputFile, bool byName = false);

#endif // RAMCORE_RAMSORT_H
