#include "rntuple/RAMTags.h"
#include "rntuple/RAMNTupleRecord.h"

#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int kRG = 4;
constexpr int kPG = 5;

// The members kRAMTypedTags names, by its index; the integer ones come first.
int32_t *IntTag(RAMNTupleRecord &rec, int i)
{
   const std::array<int32_t *, 4> members = {&rec.tag_nm, &rec.tag_as, &rec.tag_xs, &rec.tag_mq};
   return members.at(static_cast<std::size_t>(i));
}

int32_t IntTag(const RAMNTupleRecord &rec, int i)
{
   const std::array<int32_t, 4> members = {rec.tag_nm, rec.tag_as, rec.tag_xs, rec.tag_mq};
   return members.at(static_cast<std::size_t>(i));
}

// The index into kRAMTypedTags of the tag named \p name, or -1.
int NameIndex(std::string_view name) noexcept
{
   int i = 0;
   for (const char *typed : kRAMTypedTags) {
      if (name == typed)
         return i;
      i++;
   }
   return -1;
}

// The index into kRAMTypedTags of a "TG:T:value" tag that has a column, or -1.
int TypedIndex(const std::string &tag)
{
   if (tag.size() < 5 || tag[2] != ':' || tag[4] != ':')
      return -1;
   const int i = NameIndex(std::string_view(tag.data(), 2));
   return i >= 0 && tag[3] == (i < kRG ? 'i' : 'Z') ? i : -1;
}

// The value of "TG:i:value", when writing it back as a number gives the same
// text: "007" or "+7" stay text, so the round trip is exact.
bool CanonicalInt(const std::string &tag, int32_t &out)
{
   const char *text = &tag[5];
   char *end = nullptr;
   errno = 0;
   const long long value = std::strtoll(text, &end, 10);
   if (end == text || *end != '\0' || errno == ERANGE || value < std::numeric_limits<int32_t>::min() ||
       value > std::numeric_limits<int32_t>::max() || std::to_string(value) != std::string_view(text))
      return false;
   out = static_cast<int32_t>(value);
   return true;
}

// The id of \p value in \p table, through \p seen, a thread's own copy of the
// ids it has used, so only a new value takes the table's lock.
uint32_t CachedId(RAMStringIds &table, std::unordered_map<std::string, uint32_t> &seen, const std::string &value)
{
   const auto it = seen.find(value);
   if (it != seen.end())
      return it->second;
   const uint32_t id = table.Id(value);
   seen.emplace(value, id);
   return id;
}

const std::string &Value(const std::vector<std::string> &values, uint32_t id)
{
   static const std::string kEmpty;
   return id < values.size() ? values[id] : kEmpty;
}

} // namespace

uint32_t RAMStringIds::Id(const std::string &value)
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   const auto it = m_ids.find(value);
   if (it != m_ids.end())
      return it->second;
   const auto id = static_cast<uint32_t>(m_values.size());
   m_values.push_back(value);
   m_ids.emplace(value, id);
   return id;
}

std::vector<std::string> RAMStringIds::Values() const
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   return m_values;
}

void RAMTagTables::Store(RAMFileMetadata &meta) const
{
   meta.tag_orders = orders.Values();
   meta.tag_rg = rg.Values();
   meta.tag_pg = pg.Values();
}

void RAMTagWriter::Split(const std::vector<std::string> &tags, RAMNTupleRecord &rec)
{
   rec.tag_nm = rec.tag_as = rec.tag_xs = rec.tag_mq = 0;
   rec.tag_rg = rec.tag_pg = 0;
   rec.tags.clear();
   m_order.clear();
   std::array<bool, kRAMTypedTags.size()> present{};
   for (const std::string &tag : tags) {
      const int i = TypedIndex(tag);
      bool typed = i >= 0 && !present.at(static_cast<std::size_t>(i));
      if (typed && i < kRG) {
         typed = CanonicalInt(tag, *IntTag(rec, i));
      } else if (typed) {
         m_value.assign(tag, 5);
         if (i == kRG)
            rec.tag_rg = CachedId(m_tables->rg, m_rg, m_value);
         else
            rec.tag_pg = CachedId(m_tables->pg, m_pg, m_value);
      }
      if (typed) {
         present.at(static_cast<std::size_t>(i)) = true;
         m_order += kRAMTypedTags.at(static_cast<std::size_t>(i));
      } else {
         rec.tags.push_back(tag);
         m_order += kRAMTextTag;
      }
   }
   rec.tag_order = CachedId(m_tables->orders, m_orders, m_order);
}

RAMTagReader::RAMTagReader(const RAMFileMetadata &meta) : m_rg(meta.tag_rg), m_pg(meta.tag_pg)
{
   for (const std::string &order : meta.tag_orders) {
      std::vector<int8_t> tokens;
      for (std::size_t at = 0; at + 1 < order.size(); at += 2)
         tokens.push_back(static_cast<int8_t>(NameIndex(std::string_view(&order[at], 2))));
      m_orders.push_back(std::move(tokens));
   }
}

void RAMTagReader::Append(const RAMNTupleRecord &rec, std::string &out) const
{
   std::size_t text = 0;
   if (rec.tag_order < m_orders.size()) {
      for (const int8_t i : m_orders[rec.tag_order]) {
         if (i < 0) {
            if (text < rec.tags.size()) {
               out += '\t';
               out += rec.tags[text++];
            }
            continue;
         }
         out += '\t';
         out += kRAMTypedTags.at(static_cast<std::size_t>(i));
         if (i < kRG) {
            out += ":i:";
            out += std::to_string(IntTag(rec, i));
         } else {
            out += ":Z:";
            out += Value(i == kRG ? m_rg : m_pg, i == kRG ? rec.tag_rg : rec.tag_pg);
         }
      }
   }
   // Text tags the order does not account for, as in a damaged file, still come out.
   for (; text < rec.tags.size(); text++) {
      out += '\t';
      out += rec.tags[text];
   }
}
