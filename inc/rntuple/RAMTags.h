#ifndef RNTUPLE_RAMTAGS_H
#define RNTUPLE_RAMTAGS_H

#include "rntuple/RAMNTupleRecord.h"

#include <array>
#include <cstdint>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

// Standard SAM tags have columns of their own in a record (tag_nm, tag_rg and
// the rest); every other tag stays SAM text in `tags`. A record's tag_order is
// an id into the file's list of orders, each the tags of a record in input
// order: a typed tag by its name, a text tag as "**". "**RGNM" is a text tag,
// then RG, then NM. The order also says which typed tags the record has.

/// The typed tags, in the order of RAMNTupleRecord's tag_ members.
inline constexpr std::array<const char *, 6> kRAMTypedTags = {"NM", "AS", "XS", "MQ", "RG", "PG"};
/// Marks a text tag in a tag order.
inline constexpr const char *kRAMTextTag = "**";

/// Strings a file's records share, such as RG values and tag orders, by id in
/// order of first use. Thread-safe, for the converters' workers.
class RAMStringIds {
public:
   uint32_t Id(const std::string &value);
   /// A copy, taken under the lock.
   [[nodiscard]] std::vector<std::string> Values() const;

private:
   std::vector<std::string> m_values;
   std::unordered_map<std::string, uint32_t> m_ids;
   mutable std::mutex m_mutex;
};

/// The tables a converter fills for one file, or for every file of -split.
struct RAMTagTables {
   RAMStringIds orders;
   RAMStringIds rg;
   RAMStringIds pg;

   /// Copies the tables into \p meta, for RAMNTupleRecord::WriteMetadata().
   void Store(RAMFileMetadata &meta) const;
};

/// Splits SAM tag text into a record's columns. One per thread: it keeps the
/// ids it has seen, so a worker takes the tables' lock only for a new value.
class RAMTagWriter {
public:
   explicit RAMTagWriter(RAMTagTables &tables) : m_tables(&tables) {}

   /// Sets every tag member of \p rec from \p tags ("TG:T:value" each). A
   /// typed tag goes into its column only if writing the column back gives the
   /// same text; otherwise, and for a second copy of a tag, it stays text.
   void Split(const std::vector<std::string> &tags, RAMNTupleRecord &rec);

private:
   RAMTagTables *m_tables;
   std::unordered_map<std::string, uint32_t> m_orders;
   std::unordered_map<std::string, uint32_t> m_rg;
   std::unordered_map<std::string, uint32_t> m_pg;
   std::string m_order; ///< The record's tag order, reused between records.
   std::string m_value; ///< An RG or PG value, reused between records.
};

/// Writes a record's tags back as SAM text, from a file's METADATA.
class RAMTagReader {
public:
   explicit RAMTagReader(const RAMFileMetadata &meta);

   /// Appends \p rec's tags to \p out, each after a tab, in input order.
   void Append(const RAMNTupleRecord &rec, std::string &out) const;

private:
   /// Per order, the index into kRAMTypedTags of each tag, or -1 for text.
   std::vector<std::vector<int8_t>> m_orders;
   std::vector<std::string> m_rg;
   std::vector<std::string> m_pg;
};

#endif // RNTUPLE_RAMTAGS_H
