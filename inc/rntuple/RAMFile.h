#ifndef RNTUPLE_RAMFILE_H
#define RNTUPLE_RAMFILE_H

#include "rntuple/RAMNTupleRecord.h"
#include "rntuple/RAMTags.h"

#include <ROOT/RNTupleReader.hxx>

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/// A RAM file open for reading: its records and what its METADATA says. Every
/// open file has its own, so several can be open in one process.
class RAMFile {
public:
   /// nullptr, with the reason printed, for a file that is not a RAM file or
   /// has another format version.
   static std::unique_ptr<RAMFile> Open(const std::string &path);

   RAMFile(std::unique_ptr<ROOT::RNTupleReader> reader, RAMFileMetadata meta);

   ROOT::RNTupleReader &GetReader() { return *m_reader; }
   [[nodiscard]] const RAMFileMetadata &GetMetadata() const { return m_meta; }

   /// The name of reference \p refid; "*" for -1 or an id the file does not have.
   [[nodiscard]] const std::string &GetRname(int32_t refid) const { return Name(m_meta.rname_refs, refid); }
   /// The name of mate reference \p refnext, "=" included; "*" as for GetRname().
   [[nodiscard]] const std::string &GetRnext(int32_t refnext) const { return Name(m_meta.rnext_refs, refnext); }
   /// The refid of \p rname; -1 when the file has no such reference.
   [[nodiscard]] int32_t FindRefId(const std::string &rname) const;
   /// Appends \p rec's optional tags to \p out as SAM text, each after a tab,
   /// in the order the input had them.
   void AppendTags(const RAMNTupleRecord &rec, std::string &out) const { m_tags.Append(rec, out); }

private:
   static const std::string &Name(const std::vector<std::string> &names, int32_t id);

   std::unique_ptr<ROOT::RNTupleReader> m_reader;
   RAMFileMetadata m_meta;
   RAMTagReader m_tags;
   std::unordered_map<std::string, int32_t> m_rnameIds;
};

#endif // RNTUPLE_RAMFILE_H
