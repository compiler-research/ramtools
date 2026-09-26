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

/**
 * \class RAMNTupleRecord
 * \brief Alignment record stored in the ROOT Experimental RNTuple format.
 *
 * Uses standard C++ containers instead of raw C-style buffers and therefore
 * integrates naturally with the columnar storage back-end. Quality strings can
 * be stored as Phred+33, Illumina-binned, or dropped.
 *
 * Static managers (`fgRnameRefs`, `fgRnextRefs`) provide the shared
 * reference-name tables.
 *
 * \sa RAMNTupleRefs
 */
class RAMNTupleRecord {
public:
   // Quality compression options
   enum EQualCompressionBits {
      kPhred33 = 1 << 14,         // Default Phred+33 quality score
      kIlluminaBinning = 1 << 15, // Illumina 8 bin compression
      kDrop = 1 << 16,            // Drop quality score
      kQualInBlock = 1 << 18      // QUAL is in the fqzcomp block, not in `qual`
   };

   static constexpr const char *kQualBlockField = "qualblock";

   // Alignment data fields
   std::string qname;             // Query template NAME
   uint16_t flag;                 // Bitwise FLAG
   int32_t refid;                 // Reference sequence ID
   int32_t pos;                   // 0-based left most mapping POSition
   uint8_t mapq;                  // MAPing Quality
   std::vector<uint32_t> cigar;   // CIGAR operations
   int32_t refnext;               // Reference ID of the mate/next read
   int32_t pnext;                 // 0-based position of the mate/next read
   int32_t tlen;                  // Observed Template LENgth
   std::string seq;               // Segment sequence (encoded)
   std::string qual;              // Quality scores (encoded)
   std::vector<std::string> tags; // Optional SAM tags

   uint32_t compression_flags;

   // Static reference managers
   static std::unique_ptr<RAMNTupleRefs> fgRnameRefs;
   static std::unique_ptr<RAMNTupleRefs> fgRnextRefs;

   /// Longest reference span of any alignment in the file, so a region query
   /// knows how far before the region a read may start. 0 means unrecorded.
   static uint32_t fgMaxRefSpan;

   /// Order of the open file; region queries seek only when sorted. Files
   /// without the field are read as sorted.
   static RAMCoordinateOrder fgOrder;

   /// Last row of each quality block.
   static std::vector<uint64_t> fgQualBlockEnds;

public:
   RAMNTupleRecord();
   ~RAMNTupleRecord() = default;

   // Setters (SAM format, 1-based positions)
   void SetQNAME(const std::string &qname_) { qname = qname_; }
   void SetFLAG(uint16_t f) { flag = f; }
   void SetRNAME(const std::string &rname);
   void SetREFID(const std::string &rname) { SetRNAME(rname); }
   void SetPOS(int32_t pos_) { pos = pos_ - 1; } // SAM is 1-based, we are 0-based
   void SetMAPQ(uint8_t mapq_) { mapq = mapq_; }
   void SetCIGAR(const std::string &cigar_str);
   void SetRNEXT(const std::string &rnext);
   void SetREFNEXT(const std::string &rnext) { SetRNEXT(rnext); }
   void SetPNEXT(int32_t pnext_) { pnext = pnext_ - 1; } // SAM is 1-based, we are 0-based
   void SetTLEN(int32_t tlen_) { tlen = tlen_; }
   void SetSEQ(const std::string &seq_str);
   void SetQUAL(const std::string &qual_str);
   void AddTag(const std::string &tag) { tags.push_back(tag); }
   void SetOPT(const std::string &tag) { AddTag(tag); }
   void ClearTags() { tags.clear(); }
   void ResetNOPT() { ClearTags(); }

   // Getters (SAM format, 1-based positions)
   const std::string &GetQNAME() const { return qname; }
   uint16_t GetFLAG() const { return flag; }
   const std::string &GetRNAME() const;
   int32_t GetREFID() const { return refid; }
   int32_t GetPOS() const { return pos + 1; } // Convert back to 1-based for SAM
   uint8_t GetMAPQ() const { return mapq; }
   std::string GetCIGAR() const;
   const std::string &GetRNEXT() const;
   int32_t GetREFNEXT() const { return refnext; }
   int32_t GetPNEXT() const { return pnext + 1; } // Convert back to 1-based for SAM
   int32_t GetTLEN() const { return tlen; }
   std::string GetSEQ() const;
   std::string GetQUAL() const;
   const std::vector<std::string> &GetTags() const { return tags; }
   int GetNOPT() const { return static_cast<int>(tags.size()); }
   const std::string &GetOPT(int idx) const { return tags[idx]; }

   // Return sequence length without decoding (fast)
   int GetSEQLEN() const;

   size_t GetNCIGAROP() const { return cigar.size(); }
   int32_t GetCIGAROPLEN(size_t idx) const;
   int32_t GetCIGAROP(size_t idx) const;

   void Print(const char *option = "") const;
   bool IsValid() const;
   void SetBit(uint32_t bit) { compression_flags |= bit; }
   bool TestBit(uint32_t bit) const { return compression_flags & bit; }

   // Static managers
   static void InitializeRefs();
   static uint32_t GetMaxRefSpan() { return fgMaxRefSpan; }
   static void NoteRefSpan(uint32_t span)
   {
      if (span > fgMaxRefSpan)
         fgMaxRefSpan = span;
   }
   static bool IsCoordinateSorted() { return fgOrder.sorted; }
   static void SetCoordinateSorted(bool sorted) { fgOrder.sorted = sorted; }
   static void NotePlacement(int32_t refid, int32_t pos) { fgOrder.Note(refid, pos); }
   /// Reference bases covered by this record's CIGAR (0 when it has none).
   uint32_t GetRefSpan() const;
   static RAMNTupleRefs *GetRnameRefs() { return fgRnameRefs.get(); }
   static RAMNTupleRefs *GetRnextRefs() { return fgRnextRefs.get(); }
   static const std::vector<uint64_t> &GetQualBlockEnds() { return fgQualBlockEnds; }
   static void SetQualBlockEnds(std::vector<uint64_t> ends) { fgQualBlockEnds = std::move(ends); }

   // File I/O
   static std::unique_ptr<ROOT::RNTupleReader>
   OpenRAMFile(const std::string &filename, const std::string &ntupleName = "RAM");
   static void WriteAllRefs(TFile &file);
   static void ReadAllRefs(const std::string &filename = "");

   // RNTuple model creation
   static std::unique_ptr<ROOT::RNTupleModel> MakeModel();

   void SetCompressionMode(uint32_t flags) { compression_flags = flags; }

private:
   /// Creates the shared tables on first use without touching their contents.
   static void EnsureTables();
   static void WriteRefs(TFile &file, const RAMNTupleRefs *refs, const std::string &refname);
   static void ReadRefs(const std::string &filename, std::unique_ptr<RAMNTupleRefs> &refs, const std::string &refname);
   static void WriteRnameRefs(TFile &file) { WriteRefs(file, fgRnameRefs.get(), "RnameRefs"); }
   static void WriteRnextRefs(TFile &file) { WriteRefs(file, fgRnextRefs.get(), "RnextRefs"); }
   static void ReadRnameRefs(const std::string &filename = "") { ReadRefs(filename, fgRnameRefs, "RnameRefs"); }
   static void ReadRnextRefs(const std::string &filename = "") { ReadRefs(filename, fgRnextRefs, "RnextRefs"); }
};

// CIGAR operation codes (from BAM format)
#include "ramcore/CigarOps.h"

// Sequence and Quality utilities
namespace RAMNTupleUtils {
std::string EncodeSequence(const std::string &seq);
std::string DecodeSequence(const char *packed, size_t packed_size, size_t length);

std::string EncodeQuality(const std::string &qual, uint32_t compression_flags);
std::string DecodeQuality(const std::string &encoded_qual, uint32_t compression_flags);

/// Returns an empty vector and logs an error if the CIGAR is malformed.
std::vector<uint32_t> ParseCIGAR(const std::string &cigar_str);
std::string FormatCIGAR(const std::vector<uint32_t> &cigar_ops);

extern const uint8_t kIlluminaBinning[256];
} // namespace RAMNTupleUtils
#endif
