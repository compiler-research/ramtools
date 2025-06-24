#ifndef CIGAR_OPS_H
#define CIGAR_OPS_H

// CIGAR operation codes (from BAM format)
const unsigned char RAM_CIGAR_M = 0;
const unsigned char RAM_CIGAR_I = 1;
const unsigned char RAM_CIGAR_D = 2;
const unsigned char RAM_CIGAR_N = 3;
const unsigned char RAM_CIGAR_S = 4;
const unsigned char RAM_CIGAR_H = 5;
const unsigned char RAM_CIGAR_P = 6;
const unsigned char RAM_CIGAR_EQUAL = 7;
const unsigned char RAM_CIGAR_X = 8;

#endif // CIGAR_OPS_H