#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <TError.h>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

// RNTupleReader reports a missing file or RNTuple only by throwing. Checking
// with TFile first costs more than it saves: opening a TFile starts ROOT's
// interpreter, about 0.35 s and 125 MB on every open.
std::unique_ptr<ROOT::RNTupleReader> TryOpen(const std::string &ntupleName, const std::string &path)
{
   try {
      return ROOT::RNTupleReader::Open(ntupleName, path);
   } catch (const std::exception &) {
      return nullptr;
   }
}

// ROOT's Error() is printf-style; the one vararg call is here.
void Refuse(const std::string &reason)
{
   ::Error("RAMFile::Open", "%s", reason.c_str()); // NOLINT(cppcoreguidelines-pro-type-vararg)
}

// The fields RAMNTupleRecord::WriteMetadata() writes; false, with the reason
// printed, when the file has no METADATA or another format version.
bool ReadMetadata(const std::string &path, RAMFileMetadata &meta)
{
   auto reader = TryOpen("METADATA", path);
   if (!reader || reader->GetNEntries() == 0) {
      Refuse(path + " has no METADATA; convert the input again with this RAMTools");
      return false;
   }
   if (reader->GetDescriptor().FindFieldId("format_version") == ROOT::kInvalidDescriptorId) {
      Refuse(path + " was written by RAMTools 0.1.0 or earlier; convert the input again with this RAMTools");
      return false;
   }
   const uint32_t version = reader->GetView<uint32_t>("format_version")(0);
   if (version != RAMNTupleRecord::kFormatVersion) {
      Refuse(path + " has RAM format version " + std::to_string(version) + " and this RAMTools reads version " +
             std::to_string(RAMNTupleRecord::kFormatVersion) + "; convert the input again with it");
      return false;
   }

   meta.header_lines = reader->GetView<std::vector<std::string>>("header_lines")(0);
   meta.rname_refs = reader->GetView<std::vector<std::string>>("rname_refs")(0);
   meta.rnext_refs = reader->GetView<std::vector<std::string>>("rnext_refs")(0);
   meta.max_ref_span = reader->GetView<uint32_t>("max_ref_span")(0);
   meta.coordinate_sorted = reader->GetView<bool>("coordinate_sorted")(0);
   meta.qual_block_ends = reader->GetView<std::vector<uint64_t>>("qual_block_ends")(0);
   meta.counts.mapped = reader->GetView<std::vector<uint64_t>>("ref_mapped")(0);
   meta.counts.unmapped = reader->GetView<std::vector<uint64_t>>("ref_unmapped")(0);
   meta.counts.unplaced = reader->GetView<uint64_t>("unplaced")(0);
   meta.tag_orders = reader->GetView<std::vector<std::string>>("tag_orders")(0);
   meta.tag_rg = reader->GetView<std::vector<std::string>>("tag_rg")(0);
   meta.tag_pg = reader->GetView<std::vector<std::string>>("tag_pg")(0);
   return true;
}

} // namespace

std::unique_ptr<RAMFile> RAMFile::Open(const std::string &path)
{
   auto reader = TryOpen("RAM", path);
   if (!reader) {
      Refuse(path + " has no RNTuple RAM");
      return nullptr;
   }
   RAMFileMetadata meta;
   if (!ReadMetadata(path, meta))
      return nullptr;
   return std::make_unique<RAMFile>(std::move(reader), std::move(meta));
}

RAMFile::RAMFile(std::unique_ptr<ROOT::RNTupleReader> reader, RAMFileMetadata meta)
   : m_reader(std::move(reader)), m_meta(std::move(meta)), m_tags(m_meta)
{
   for (std::size_t i = 0; i < m_meta.rname_refs.size(); i++)
      m_rnameIds.emplace(m_meta.rname_refs[i], static_cast<int32_t>(i));
}

int32_t RAMFile::FindRefId(const std::string &rname) const
{
   const auto it = m_rnameIds.find(rname);
   return it != m_rnameIds.end() ? it->second : -1;
}

const std::string &RAMFile::Name(const std::vector<std::string> &names, int32_t id)
{
   static const std::string kStar = "*";
   if (id < 0 || static_cast<std::size_t>(id) >= names.size())
      return kStar;
   return names[static_cast<std::size_t>(id)];
}
