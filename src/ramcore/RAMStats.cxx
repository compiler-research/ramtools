#include "ramcore/RAMStats.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <TFile.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <map>

namespace ramcore {

static constexpr uint16_t FLAG_PAIRED         = 0x1;
static constexpr uint16_t FLAG_PROPER_PAIR    = 0x2;
static constexpr uint16_t FLAG_UNMAPPED       = 0x4;
static constexpr uint16_t FLAG_REVERSE_STRAND = 0x10;
static constexpr uint16_t FLAG_DUPLICATE      = 0x400;

void RAMStats::Print() const
{
   auto pct = [](uint64_t a, uint64_t b) -> double {
      return b == 0 ? 0.0 : 100.0 * static_cast<double>(a) / static_cast<double>(b);
   };

   std::cout << "\n=== RAMTools Statistics ===\n";
   std::cout << std::left
             << std::setw(30) << "Total reads:"           << total_reads << "\n"
             << std::setw(30) << "Mapped reads:"          << mapped_reads
             << std::fixed << std::setprecision(2)
             << "  (" << pct(mapped_reads, total_reads) << "%)\n"
             << std::setw(30) << "Unmapped reads:"        << unmapped_reads
             << "  (" << pct(unmapped_reads, total_reads) << "%)\n"
             << std::setw(30) << "Duplicate reads:"       << duplicate_reads
             << "  (" << pct(duplicate_reads, total_reads) << "%)\n"
             << std::setw(30) << "Paired reads:"          << paired_reads
             << "  (" << pct(paired_reads, total_reads) << "%)\n"
             << std::setw(30) << "Properly paired reads:" << properly_paired_reads
             << "  (" << pct(properly_paired_reads, total_reads) << "%)\n"
             << std::setw(30) << "Forward strand:"        << forward_strand
             << "  (" << pct(forward_strand, mapped_reads) << "% of mapped)\n"
             << std::setw(30) << "Reverse strand:"        << reverse_strand
             << "  (" << pct(reverse_strand, mapped_reads) << "% of mapped)\n"
             << std::setw(30) << "Total bases:"           << total_bases << "\n"
             << std::setw(30) << "Mean read length:"
             << std::fixed << std::setprecision(2) << mean_read_length << "\n"
             << std::setw(30) << "Mean mapping quality:"
             << std::fixed << std::setprecision(2) << mean_mapping_quality << "\n";

   std::cout << "===========================\n\n";
}

/// Extract the original sequence length from an encoded seq field.
/// The field format is: [4-byte uint32_t length][packed 2-bit nucleotides]
/// We use memcpy to avoid undefined behaviour from unaligned reinterpret_cast.
static uint32_t DecodedSeqLength(const std::string &encoded_seq)
{
   if (encoded_seq.size() < 4)
      return 0;
   uint32_t length = 0;
   std::memcpy(&length, encoded_seq.data(), sizeof(length));
   return length;
}

RAMStatsResult ComputeStats(const char *filename)
{
   RAMStats stats;

   // Load reference name map stored alongside the RNTuple
   RAMNTupleRecord::InitializeRefs();
   RAMNTupleRecord::ReadAllRefs(filename);
   RAMNTupleRefs *rnameRefs = RAMNTupleRecord::GetRnameRefs();

   std::unique_ptr<ROOT::RNTupleReader> reader;
   try {
      reader = ROOT::RNTupleReader::Open("RAM", filename);
   } catch (const std::exception &e) {
      // Return error result so caller can distinguish bad file from empty file
      return RAMStatsResult{stats, false, e.what()};
   }

   if (!reader)
      return RAMStatsResult{stats, false, "RNTupleReader::Open returned nullptr"};

   auto viewFlag  = reader->GetView<uint16_t>("record.flag");
   auto viewMapQ  = reader->GetView<uint8_t>("record.mapq");
   auto viewSeq   = reader->GetView<std::string>("record.seq");
   auto viewRefId = reader->GetView<int32_t>("record.refid");

   uint64_t mapq_sum = 0;
   uint64_t len_sum  = 0;
   const uint64_t nEntries = reader->GetNEntries();

   for (uint64_t i = 0; i < nEntries; ++i) {
      stats.total_reads++;

      uint16_t flag  = viewFlag(i);
      uint8_t  mapq  = viewMapQ(i);
      const std::string &seq = viewSeq(i);
      int32_t  refid = viewRefId(i);

      if (flag & FLAG_UNMAPPED) {
         stats.unmapped_reads++;
         // Unmapped reads have no strand — do NOT count them in strand stats
      } else {
         stats.mapped_reads++;
         mapq_sum += mapq;
         // Strand is only meaningful for mapped reads
         if (flag & FLAG_REVERSE_STRAND)
            stats.reverse_strand++;
         else
            stats.forward_strand++;
      }

      if (flag & FLAG_DUPLICATE)   stats.duplicate_reads++;
      if (flag & FLAG_PAIRED)      stats.paired_reads++;
      if (flag & FLAG_PROPER_PAIR) stats.properly_paired_reads++;

      // Extract sequence length from encoded field (4-byte length prefix + packed nibbles)
      uint32_t rlen = DecodedSeqLength(seq);
      len_sum           += rlen;
      stats.total_bases += rlen;

      // Resolve chromosome name from integer refid via the ref name table
      if (refid >= 0 && rnameRefs &&
          refid < static_cast<int32_t>(rnameRefs->Size())) {
         const std::string &chrom = rnameRefs->GetRefName(refid);
         if (!chrom.empty() && chrom != "*")
            stats.reads_per_chromosome[chrom]++;
      }
   }

   if (stats.total_reads > 0)
      stats.mean_read_length = static_cast<double>(len_sum) / stats.total_reads;
   if (stats.mapped_reads > 0)
      stats.mean_mapping_quality = static_cast<double>(mapq_sum) / stats.mapped_reads;

   return RAMStatsResult{stats, true, ""};
}

} // namespace ramcore
