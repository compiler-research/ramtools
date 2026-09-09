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
#include <map>
#include <memory>
#include <cstdint>

class RAMNTupleRefs;
class RAMNTupleIndex;

/**
 * \class RAMNTupleRefs
 * \brief Reference-name ↔ ID mapping used by the RNTuple backend.
 *
 * Conceptually identical to `RAMRefs` but serialisable as an RNTuple field so
 * that it can be embedded directly inside the columnar file.
 */
class RAMNTupleRefs {
private:
   std::vector<std::string> fRefVec;
   mutable int fLastId;
   mutable std::string fLastName;

public:
   RAMNTupleRefs();
   ~RAMNTupleRefs() = default;

   int GetRefId(const std::string &rname);
   int FindRefId(const std::string &rname) const; ///< Lookup-only; returns -1 if not found.
   const std::string &GetRefName(int rid) const;

   void Print() const;
   size_t Size() const { return fRefVec.size(); }

   // For RNTuple serialization
   void Clear()
   {
      fRefVec.clear();
      fLastId = -1;
      fLastName.clear();
   }
   void AddRef(const std::string &ref) { fRefVec.push_back(ref); }
   const std::vector<std::string> &GetRefs() const { return fRefVec; }
   void SetRefs(const std::vector<std::string> &refs)
   {
      fRefVec = refs;
      fLastId = -1;
      fLastName.clear();
   }
};

/**
 * \class RAMNTupleIndex
 * \brief Sparse genomic index for fast region queries on RNTuple files.
 *
 * The index is stored as a plain `std::vector<IndexEntry>` for efficient
 * serialisation and complemented by a lazily-initialised `std::map` that
 * supports O(log n) look-ups by (refid,pos).
 */
class RAMNTupleIndex {
public:
   struct IndexEntry {
      int32_t refid;
      int32_t pos;
      int64_t entry;
   };

private:
   std::vector<IndexEntry> fIndex;
   mutable std::map<std::pair<int32_t, int32_t>, int64_t> fIndexMap;

   void RebuildMap() const;

public:
   RAMNTupleIndex() = default;
   ~RAMNTupleIndex() = default;

   void AddItem(int32_t refid, int32_t pos, int64_t row);
   int64_t GetRow(int32_t refid, int32_t pos) const;
   std::vector<int64_t> GetRowsInRange(int32_t refid, int32_t start, int32_t end) const;

   void Print() const;
   size_t Size() const { return fIndex.size(); }

   // For RNTuple serialization
   const std::vector<IndexEntry> &GetEntries() const { return fIndex; }
   void SetEntries(const std::vector<IndexEntry> &entries)
   {
      fIndex = entries;
      fIndexMap.clear();
   }
   void Clear()
   {
      fIndex.clear();
      fIndexMap.clear();
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
 * Static managers (`fgRnameRefs`, `fgRnextRefs`, `fgIndex`) provide shared
 * metadata for reference names and sparse region lookup.
 *
 * \sa RAMNTupleRefs, RAMNTupleIndex
 */
class RAMNTupleRecord {
public:
   // Quality compression options
   enum EQualCompressionBits {
      kPhred33 = 1 << 14,         // Default Phred+33 quality score
      kIlluminaBinning = 1 << 15, // Illumina 8 bin compression
      kDrop = 1 << 16             // Drop quality score
   };

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

   // Static reference and index managers
   static std::unique_ptr<RAMNTupleRefs> fgRnameRefs;
   static std::unique_ptr<RAMNTupleRefs> fgRnextRefs;
   static std::unique_ptr<RAMNTupleIndex> fgIndex;

   /// Longest reference span of any alignment in the file, so a region query
   /// knows how far before the region a read may start. 0 means unrecorded.
   static uint32_t fgMaxRefSpan;

   /// Whether the records are in coordinate order. A region query can only
   /// seek to an index entry and stop at the first record past the region when
   /// they are; on an unsorted file it reads every record. Files written before
   /// this field carry no answer and are read as sorted, as they always were.
   static bool fgCoordinateSorted;
   static int32_t fgLastPlacedRefId;
   static int32_t fgLastPlacedPos;

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
   static bool IsCoordinateSorted() { return fgCoordinateSorted; }
   static void SetCoordinateSorted(bool sorted) { fgCoordinateSorted = sorted; }
   /// Feeds one placed record to the running order check.
   static void NotePlacement(int32_t refid_, int32_t pos_)
   {
      if (refid_ < fgLastPlacedRefId || (refid_ == fgLastPlacedRefId && pos_ < fgLastPlacedPos))
         fgCoordinateSorted = false;
      fgLastPlacedRefId = refid_;
      fgLastPlacedPos = pos_;
   }
   /// Reference bases covered by this record's CIGAR (0 when it has none).
   uint32_t GetRefSpan() const;
   static RAMNTupleRefs *GetRnameRefs() { return fgRnameRefs.get(); }
   static RAMNTupleRefs *GetRnextRefs() { return fgRnextRefs.get(); }
   static RAMNTupleIndex *GetIndex() { return fgIndex.get(); }

   // File I/O
   static std::unique_ptr<ROOT::RNTupleReader>
   OpenRAMFile(const std::string &filename, const std::string &ntupleName = "RAM");
   static void WriteAllRefs(TFile &file);
   static void ReadAllRefs(const std::string &filename = "");
   static void WriteIndex(TFile &file);
   static void WriteIndex(TFile &file, const RAMNTupleIndex &index);
   static void ReadIndex(const std::string &filename = "");

   // RNTuple model creation
   static std::unique_ptr<ROOT::RNTupleModel> MakeModel();

   void SetCompressionMode(uint32_t flags) { compression_flags = flags; }

private:
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
