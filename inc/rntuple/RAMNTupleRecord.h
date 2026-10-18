//
// RAMNTupleRecord.h
// Header for RAM (ROOT Alignment/Map) format

#ifndef RAMNTupleRecord_h
#define RAMNTupleRecord_h

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RNTupleView.hxx>

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>
#include <cstddef>
#include <cstdint>

class RAMNTupleRefs;

/**
 * \class RAMNTupleRefs
 * \brief Reference-name ↔ ID mapping used by the RNTuple backend.
 *
 * Conceptually identical to `RAMRefs` but serialisable as an RNTuple field so
 * that it can be embedded directly inside the columnar file.
 */
class RAMNTupleRefs {
private:
   std::vector<std::string> m_refVec;
   std::unordered_map<std::string, int> m_index;
   mutable std::mutex m_mutex;

   void Rebuild();

public:
   RAMNTupleRefs();
   ~RAMNTupleRefs() = default;
   RAMNTupleRefs(const RAMNTupleRefs &) = delete;
   RAMNTupleRefs &operator=(const RAMNTupleRefs &) = delete;

   /// Returns the id of \p rname, adding it if new. Thread-safe.
   int GetRefId(const std::string &rname);
   /// Lookup-only; returns -1 if not found. Thread-safe.
   int FindRefId(const std::string &rname) const;
   /// The reference is invalidated by a concurrent GetRefId().
   const std::string &GetRefName(int rid) const;

   void Print() const;
   size_t Size() const;

   // For RNTuple serialization
   void Clear();
   void AddRef(const std::string &ref);
   /// A copy, taken under the lock. Thread-safe.
   std::vector<std::string> GetRefs() const;
   void SetRefs(const std::vector<std::string> &refs);
};

/// Running check of coordinate order: placed records ordered by (refid, pos),
/// unplaced ones (refid -1) after them.
struct RAMCoordinateOrder {
   int32_t last_refid = -1;
   int32_t last_pos = -1;
   bool seen_unplaced = false;
   bool sorted = true;

   void Note(int32_t refid, int32_t pos)
   {
      if (refid < 0) {
         seen_unplaced = true;
         return;
      }
      if (seen_unplaced || refid < last_refid || (refid == last_refid && pos < last_pos))
         sorted = false;
      last_refid = refid;
      last_pos = pos;
   }
};

/// Records per reference, counted as `samtools idxstats` counts them: a placed
/// record by whether FLAG 0x4 is set, a record with no reference as unplaced.
struct RAMRefCounts {
   std::vector<uint64_t> mapped;   ///< By refid.
   std::vector<uint64_t> unmapped; ///< By refid.
   uint64_t unplaced = 0;

   void Note(int32_t refid, uint16_t flag)
   {
      if (refid < 0) {
         unplaced++;
         return;
      }
      Resize(static_cast<size_t>(refid) + 1);
      ((flag & 0x4) ? unmapped : mapped)[static_cast<size_t>(refid)]++;
   }

   void Add(const RAMRefCounts &other)
   {
      Resize(other.mapped.size());
      for (size_t i = 0; i < other.mapped.size(); i++) {
         mapped[i] += other.mapped[i];
         unmapped[i] += other.unmapped[i];
      }
      unplaced += other.unplaced;
   }

   /// Grows both vectors to at least \p n references; never shrinks them.
   void Resize(size_t n)
   {
      if (n > mapped.size()) {
         mapped.resize(n);
         unmapped.resize(n);
      }
   }
};

/// What a RAM file's METADATA holds, one value per file.
struct RAMFileMetadata {
   std::vector<std::string> header_lines; ///< The SAM header, a line per entry, without the newline.
   std::vector<std::string> rname_refs;   ///< Reference names by refid.
   std::vector<std::string> rnext_refs;   ///< Mate reference names by refnext, "=" among them.
   /// Longest reference span of any record, so a region query knows how far
   /// before the region a read may start.
   uint32_t max_ref_span = 0;
   /// Region queries seek only in a sorted file.
   bool coordinate_sorted = true;
   /// Last row of each block of records; that row carries the block's names
   /// and qualities. See ramcore/RecordBlocks.h.
   std::vector<uint64_t> block_ends;
   RAMRefCounts counts;                   ///< One entry per reference.
   std::vector<std::string> tag_orders;   ///< The orders tag_order points into; see RAMTags.h.
   std::vector<std::string> tag_rg;       ///< RG values by tag_rg.
   std::vector<std::string> tag_pg;       ///< PG values by tag_pg.
};

/**
 * \class RAMNTupleRecord
 * \brief Alignment record stored in the ROOT Experimental RNTuple format.
 *
 * Uses standard C++ containers instead of raw C-style buffers and therefore
 * integrates naturally with the columnar storage back-end. Quality strings can
 * be stored as Phred+33, Illumina-binned, or dropped.
 *
 * A record holds reference ids, not names. The names are in its file's
 * METADATA: a converter fills them from its own RAMNTupleRefs tables, and a
 * reader looks them up in RAMFile.
 *
 * \sa RAMNTupleRefs, RAMFile
 */
class RAMNTupleRecord {
public:
   // Quality compression options
   enum EQualCompressionBits {
      kPhred33 = 1 << 14,         // Default Phred+33 quality score
      kIlluminaBinning = 1 << 15, // Illumina 8 bin compression
      kDrop = 1 << 16,            // Drop quality score
      kQualInBlock = 1 << 18,     // QUAL is in its block, not in `qual`
      kFqzcomp = 1 << 19          // In a quality policy: compress the quality blocks
                                  // with fqzcomp, not rANS (-fqzcomp); not kept in records
   };

   static constexpr const char *kQualBlockField = "qualblock";
   static constexpr const char *kNameBlockField = "nameblock";

   /// Version of the file layout, stored in METADATA. Until 1.0 a reader opens
   /// only files of its own version; bump it with every change to the layout.
   static constexpr uint32_t kFormatVersion = 8;

   // Alignment data fields
   std::string qname;             //! Query template NAME; not a column: the names are in kNameBlockField
   uint16_t flag;                 // Bitwise FLAG
   int32_t refid;                 // Reference sequence ID
   int32_t pos;                   // 0-based left most mapping POSition
   uint8_t mapq;                  // MAPing Quality
   std::vector<uint32_t> cigar;   // CIGAR operations
   int32_t refnext;               // Reference ID of the mate/next read
   // PNEXT and TLEN, stored so that they compress; SetMate() fills them and
   // GetPNEXT() and GetTLEN() give the SAM values back.
   int32_t pnext_delta = 0;       // 0-based PNEXT - pos when the mate is on this reference, else 0-based PNEXT
   int32_t tlen_delta = 0;        // TLEN - PredictTLEN()
   std::string seq;               // Segment sequence as uppercase IUPAC text, or "*"
   std::string qual;              // Quality scores (encoded); empty when kQualInBlock

   // Standard tags in columns of their own; RAMTagWriter fills them and
   // RAMTagReader writes them back. tag_order says which the record has.
   int32_t tag_nm = 0;
   int32_t tag_as = 0;
   int32_t tag_xs = 0;
   int32_t tag_mq = 0;
   uint32_t tag_rg = 0; // Id into METADATA's tag_rg
   uint32_t tag_pg = 0; // Id into METADATA's tag_pg
   std::string tag_md;
   std::string tag_mc;
   std::string tag_sa;
   std::string tag_xa;
   uint32_t tag_order = 0;
   std::vector<std::string> tags; // The other optional tags, as SAM text

   uint32_t compression_flags;

   RAMNTupleRecord();
   ~RAMNTupleRecord() = default;

   // Setters (SAM format, 1-based positions)
   void SetQNAME(const std::string &qname_) { qname = qname_; }
   void SetFLAG(uint16_t f) { flag = f; }
   void SetPOS(int32_t pos_) { pos = pos_ - 1; } // SAM is 1-based, we are 0-based
   void SetMAPQ(uint8_t mapq_) { mapq = mapq_; }
   void SetCIGAR(const std::string &cigar_str);
   /// Sets PNEXT (1-based) and TLEN, after SetFLAG(), SetPOS() and SetCIGAR():
   /// they are stored against pos, the strands and the CIGAR. \p sameRef is
   /// RAMMateOnSameRef() of the record's RNAME and RNEXT.
   void SetMate(int32_t pnext, int32_t tlen, bool sameRef);
   void SetSEQ(const std::string &seq_str);
   void SetQUAL(const std::string &qual_str);

   // Getters (SAM format, 1-based positions)
   /// Empty in a record read from a file: read the names with NameBlockReader.
   const std::string &GetQNAME() const { return qname; }
   uint16_t GetFLAG() const { return flag; }
   int32_t GetREFID() const { return refid; }
   int32_t GetPOS() const { return pos + 1; } // Convert back to 1-based for SAM
   uint8_t GetMAPQ() const { return mapq; }
   std::string GetCIGAR() const;
   int32_t GetREFNEXT() const { return refnext; }
   /// PNEXT (1-based) and TLEN as SetMate() was given them; \p sameRef as
   /// there, or RAMFile::MateOnSameRef().
   [[nodiscard]] int32_t GetPNEXT(bool sameRef) const;
   [[nodiscard]] int32_t GetTLEN(bool sameRef) const;
   /// The TLEN of a pair whose mate has the same reference span as this
   /// record: the distance between their 5' ends, as BWA and Picard compute
   /// it. 0 when either is unmapped or the mate is on another reference.
   [[nodiscard]] int32_t PredictTLEN(int32_t pnext0, bool sameRef) const;
   [[nodiscard]] const std::string &GetSEQ() const { return seq; }
   /// Empty in a record read from a file whose quality is in a block
   /// (kQualInBlock): read the qualities with QualityBlockReader.
   std::string GetQUAL() const;

   /// Number of bases; 0 for "*".
   int GetSEQLEN() const;

   size_t GetNCIGAROP() const { return cigar.size(); }
   int32_t GetCIGAROPLEN(size_t idx) const;
   int32_t GetCIGAROP(size_t idx) const;

   bool IsValid() const;
   void SetBit(uint32_t bit) { compression_flags |= bit; }
   bool TestBit(uint32_t bit) const { return compression_flags & bit; }

   /// Reference bases covered by this record's CIGAR (0 when it has none).
   uint32_t GetRefSpan() const;

   // File I/O
   /// Writes METADATA, with the format version, into \p file. A converter
   /// gathers \p meta while it writes the records and calls this once at the end.
   static void WriteMetadata(TFile &file, const RAMFileMetadata &meta);

   // RNTuple model creation
   static std::unique_ptr<ROOT::RNTupleModel> MakeModel();

   /// Replaces the quality mode.
   void SetCompressionMode(uint32_t flags) { compression_flags = flags; }
};

/// Whether a record's mate is on the record's own reference: RNEXT is "=" or
/// RNAME, and RNAME is not "*". Writers and readers decide it the same way.
inline bool RAMMateOnSameRef(std::string_view rname, std::string_view rnext)
{
   return rname != "*" && (rnext == "=" || rnext == rname);
}

// CIGAR operation codes (from BAM format)
#include "ramcore/CigarOps.h"

// Sequence and Quality utilities
namespace RAMNTupleUtils {
std::string NormalizeSequence(const std::string &seq);

std::string EncodeQuality(const std::string &qual, uint32_t compression_flags);
std::string DecodeQuality(const std::string &encoded_qual, uint32_t compression_flags);

/// Returns an empty vector and logs an error if the CIGAR is malformed.
std::vector<uint32_t> ParseCIGAR(const std::string &cigar_str);
std::string FormatCIGAR(const std::vector<uint32_t> &cigar_ops);

extern const uint8_t kIlluminaBinning[256];
} // namespace RAMNTupleUtils
#endif
