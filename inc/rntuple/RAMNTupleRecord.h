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
#include <cstddef>
#include <cstdint>

class RAMNTupleRefs;

/**
 * \class RAMNTupleRefs
 * \brief Table of reference names, stored in the METADATA ntuple.
 *
 * A record stores the index of its reference name in this table rather than
 * the name itself. Two tables exist per file: one for RNAME and one for RNEXT,
 * so that "=" can be an entry of its own.
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

   /// Id of \a rname, adding it to the table when it is new. "*" is -1. Thread-safe.
   int GetRefId(const std::string &rname);
   /// Id of \a rname without adding it; -1 when absent or "*". Thread-safe.
   int FindRefId(const std::string &rname) const;
   /// Name for an id; "*" for -1 and for ids outside the table. The reference
   /// is invalidated by a concurrent GetRefId().
   const std::string &GetRefName(int rid) const;

   void Print() const;  ///< prints the table on standard output
   size_t Size() const; ///< number of names

   /// Empties the table.
   void Clear();
   /// Appends a name without checking for duplicates; used when reading METADATA.
   void AddRef(const std::string &ref);
   /// The names in id order, copied under the lock. Thread-safe.
   std::vector<std::string> GetRefs() const;
   /// Replaces the whole table.
   void SetRefs(const std::vector<std::string> &refs);
};

/// Running check of coordinate order: placed records ordered by (refid, pos),
/// unplaced ones (refid -1) after them.
struct RAMCoordinateOrder {
   int32_t last_refid = -1;    ///< reference id of the last placed record
   int32_t last_pos = -1;      ///< position of the last placed record
   bool seen_unplaced = false; ///< an unplaced record has been seen
   bool sorted = true;         ///< every record so far was in order

   /// Checks the next record against the ones before it.
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
 * \brief One alignment record, the "record" field of the RAM ntuple.
 *
 * The setters take SAM values (1-based positions, text CIGAR, plain bases and
 * Phred+33 quality) and store the encoded form; the getters give the SAM
 * values back. RNTuple stores each data member as its own column, so a reader
 * that needs only positions and CIGARs reads only those.
 *
 * Reference names, the longest span and the sort flag are shared per process
 * through static members (`fgRnameRefs`, `fgRnextRefs`, `fgMaxRefSpan`,
 * `fgOrder`); InitializeRefs() resets the per-file parts before each
 * conversion or open.
 *
 * \sa RAMNTupleRefs, RAMCoordinateOrder, RAMNTupleUtils
 */
class RAMNTupleRecord {
public:
   /// How the QUAL column is stored; one bit is set in compression_flags.
   enum EQualCompressionBits {
      kPhred33 = 1 << 14,         ///< verbatim Phred+33 text
      kIlluminaBinning = 1 << 15, ///< one byte per base holding the Illumina-binned Phred value
      kDrop = 1 << 16             ///< nothing stored; reads back as "*"
   };

   std::string qname;             ///< read name, verbatim
   uint16_t flag;                 ///< SAM FLAG
   int32_t refid;                 ///< index into the RNAME table; -1 for "*"
   int32_t pos;                   ///< 0-based leftmost position; -1 when unplaced
   uint8_t mapq;                  ///< mapping quality
   std::vector<uint32_t> cigar;   ///< operations packed as (length << 4) | op, BAM style
   int32_t refnext;               ///< index into the RNEXT table; -1 for "*"
   int32_t pnext;                 ///< 0-based mate position; -1 when none
   int32_t tlen;                  ///< template length
   std::string seq;               ///< 4-byte LE length then bases packed two per byte; empty for "*"
   std::string qual;              ///< quality, encoded per compression_flags
   std::vector<std::string> tags; ///< optional fields verbatim, e.g. "NM:i:0"

   uint32_t compression_flags; ///< one of EQualCompressionBits

   // Static reference managers
   static std::unique_ptr<RAMNTupleRefs> fgRnameRefs;
   static std::unique_ptr<RAMNTupleRefs> fgRnextRefs;

   /// Longest reference span of any alignment in the file, so a region query
   /// knows how far before the region a read may start. 0 means unrecorded.
   static uint32_t fgMaxRefSpan;

   /// Order of the open file; region queries seek only when sorted. Files
   /// without the field are read as sorted.
   static RAMCoordinateOrder fgOrder;

public:
   RAMNTupleRecord();
   ~RAMNTupleRecord() = default;

   // Setters take SAM values: 1-based positions, text CIGAR, plain bases, Phred+33 quality.
   void SetQNAME(const std::string &qname_) { qname = qname_; } ///< read name
   void SetFLAG(uint16_t f) { flag = f; }                       ///< SAM FLAG
   /// Reference name, added to the RNAME table when new; "*" stores -1.
   void SetRNAME(const std::string &rname);
   void SetREFID(const std::string &rname) { SetRNAME(rname); }
   void SetPOS(int32_t pos_) { pos = pos_ - 1; } ///< 1-based position, stored 0-based; 0 means unplaced
   void SetMAPQ(uint8_t mapq_) { mapq = mapq_; } ///< mapping quality
   /// Text CIGAR, or "*". A malformed string is stored as no CIGAR; see RAMNTupleUtils::ParseCIGAR().
   void SetCIGAR(const std::string &cigar_str);
   /// Mate reference name, "=" or "*", added to the RNEXT table when new.
   void SetRNEXT(const std::string &rnext);
   void SetREFNEXT(const std::string &rnext) { SetRNEXT(rnext); }
   void SetPNEXT(int32_t pnext_) { pnext = pnext_ - 1; } ///< 1-based mate position, stored 0-based
   void SetTLEN(int32_t tlen_) { tlen = tlen_; }         ///< template length
   /// Bases, or "*"; stored packed, see RAMNTupleUtils::EncodeSequence().
   void SetSEQ(const std::string &seq_str);
   /// Phred+33 text, or "*"; stored according to compression_flags.
   void SetQUAL(const std::string &qual_str);
   void AddTag(const std::string &tag) { tags.push_back(tag); } ///< appends an optional field such as "NM:i:0"
   void SetOPT(const std::string &tag) { AddTag(tag); }
   void ClearTags() { tags.clear(); } ///< removes every optional field
   void ResetNOPT() { ClearTags(); }

   // Getters give SAM values back: 1-based positions, text CIGAR, decoded bases and quality.
   const std::string &GetQNAME() const { return qname; }            ///< read name
   uint16_t GetFLAG() const { return flag; }                        ///< SAM FLAG
   const std::string &GetRNAME() const;                             ///< reference name, "*" when unplaced
   int32_t GetREFID() const { return refid; }                       ///< id in the RNAME table, -1 for "*"
   int32_t GetPOS() const { return pos + 1; }                       ///< 1-based position, 0 when unplaced
   uint8_t GetMAPQ() const { return mapq; }                         ///< mapping quality
   std::string GetCIGAR() const;                                    ///< text CIGAR, "*" when there is none
   const std::string &GetRNEXT() const;                             ///< mate reference name, "=" or "*"
   int32_t GetREFNEXT() const { return refnext; }                   ///< id in the RNEXT table, -1 for "*"
   int32_t GetPNEXT() const { return pnext + 1; }                   ///< 1-based mate position, 0 when none
   int32_t GetTLEN() const { return tlen; }                         ///< template length
   std::string GetSEQ() const;                                      ///< bases, "*" when none were stored
   std::string GetQUAL() const;                                     ///< Phred+33 text; "*" when dropped or absent
   const std::vector<std::string> &GetTags() const { return tags; } ///< optional fields, verbatim
   int GetNOPT() const { return static_cast<int>(tags.size()); }
   const std::string &GetOPT(int idx) const { return tags[idx]; }

   /// Number of bases, read from the length prefix without decoding.
   int GetSEQLEN() const;

   size_t GetNCIGAROP() const { return cigar.size(); } ///< number of CIGAR operations
   int32_t GetCIGAROPLEN(size_t idx) const;            ///< length of operation \a idx
   int32_t GetCIGAROP(size_t idx) const;               ///< code of operation \a idx, see CigarOps.h

   /// Prints the record as one SAM line on standard output.
   void Print(const char *option = "") const;
   bool IsValid() const;
   void SetBit(uint32_t bit) { compression_flags |= bit; }              ///< sets a bit of compression_flags
   bool TestBit(uint32_t bit) const { return compression_flags & bit; } ///< tests a bit of compression_flags

   /// Creates the shared tables on first use and resets the per-file state:
   /// longest span and sort flag. Only the writers and OpenRAMFile() call it;
   /// the constructor must not.
   static void InitializeRefs();
   static uint32_t GetMaxRefSpan() { return fgMaxRefSpan; } ///< longest span in the open or written file
   /// Widens the longest span seen in the file being written.
   static void NoteRefSpan(uint32_t span)
   {
      if (span > fgMaxRefSpan)
         fgMaxRefSpan = span;
   }
   /// Whether the open or written file is in coordinate order.
   static bool IsCoordinateSorted() { return fgOrder.sorted; }
   /// Overrides the running order check, e.g. with a split file's own answer.
   static void SetCoordinateSorted(bool sorted) { fgOrder.sorted = sorted; }
   /// Feeds one record to the running order check.
   static void NotePlacement(int32_t refid, int32_t pos) { fgOrder.Note(refid, pos); }
   /// Reference bases covered by this record's CIGAR (0 when it has none).
   uint32_t GetRefSpan() const;
   static RAMNTupleRefs *GetRnameRefs() { return fgRnameRefs.get(); } ///< the RNAME table
   static RAMNTupleRefs *GetRnextRefs() { return fgRnextRefs.get(); } ///< the RNEXT table

   /// Opens the RAM ntuple of \a filename and loads its METADATA into the
   /// shared tables, so that name getters and region scans work.
   /// Returns null and logs an error when the file cannot be opened.
   static std::unique_ptr<ROOT::RNTupleReader>
   OpenRAMFile(const std::string &filename, const std::string &ntupleName = "RAM");
   /// Writes the METADATA ntuple: both name tables, the longest span and the sort flag.
   static void WriteAllRefs(TFile &file);
   /// Loads METADATA into the shared tables; fields missing from old files keep their defaults.
   static void ReadAllRefs(const std::string &filename = "");

   /// The RNTuple model of a RAM file: one field, "record", of this type.
   static std::unique_ptr<ROOT::RNTupleModel> MakeModel();

   /// Replaces compression_flags.
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

/// Encoders and decoders behind the RAMNTupleRecord setters and getters.
namespace RAMNTupleUtils {
/// Packs bases two per byte, 4 bits each in htslib's "=ACMGRSVTWYHKDBN" order,
/// behind a 4-byte little-endian length. Lower case is accepted; any other
/// byte becomes N. "*" becomes an empty string.
std::string EncodeSequence(const std::string &seq);
/// Inverse of EncodeSequence() for a payload without its length prefix.
/// Logs an error and returns an empty string when the payload is too short.
std::string DecodeSequence(const char *packed, size_t packed_size, size_t length);

/// Applies the policy in \a compression_flags: kPhred33 keeps the text,
/// kIlluminaBinning stores binned Phred values one byte per base ("*" becomes
/// empty), kDrop stores nothing.
std::string EncodeQuality(const std::string &qual, uint32_t compression_flags);
/// Inverse of EncodeQuality(); kDrop and an empty binned string give "*".
std::string DecodeQuality(const std::string &encoded_qual, uint32_t compression_flags);

/// Packs a text CIGAR as (length << 4) | op per operation. "*" and "" give an
/// empty vector; a malformed string logs an error and also gives an empty
/// vector, so callers that must reject it validate first (SamParser does).
std::vector<uint32_t> ParseCIGAR(const std::string &cigar_str);
/// Text CIGAR from packed operations; "*" for none.
std::string FormatCIGAR(const std::vector<uint32_t> &cigar_ops);

/// Phred value to its Illumina 8-level bin (0, 1, 6, 15, 22, 27, 33, 37, 40).
extern const uint8_t kIlluminaBinning[256];
} // namespace RAMNTupleUtils
#endif
