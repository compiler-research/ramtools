// ramdump: write a RAM file back out as SAM, the way `samtools view` does.
// `ramdump -h file.ram` should reproduce the SAM the file was built from, and
// `ramdump -c region` should agree with `samtools view -c region`.

#include "ramcore/RecordBlocks.h"
#include "ramcore/RAMNTupleView.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RFieldBase.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleReadOptions.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>
#include <Rtypes.h>
#include <TROOT.h>

#include <algorithm>
#include <cerrno>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

// Fields go into one buffer handed to fwrite in 1 MB blocks; formatting
// millions of records through std::ostream costs more than the query itself.
struct SamWriter {
   static constexpr size_t kFlushAt = 1U << 20;

   FILE *out;
   std::string buf{};

   explicit SamWriter(FILE *f) : out(f) { buf.reserve(kFlushAt + 4096); }
   ~SamWriter() { Flush(); }
   SamWriter(const SamWriter &) = delete;
   SamWriter(SamWriter &&) = delete;
   SamWriter &operator=(const SamWriter &) = delete;
   SamWriter &operator=(SamWriter &&) = delete;

   void Str(std::string_view s) { buf += s; }
   void Tab() { buf += '\t'; }
   void Int(long long v) { RAMNTupleUtils::AppendInt(buf, v); }

   void EndLine()
   {
      buf += '\n';
      if (buf.size() >= kFlushAt)
         Flush();
   }

   void Flush()
   {
      fwrite(buf.data(), sizeof(char), buf.size(), out);
      buf.clear();
   }

   /// Writes \p text after what is buffered.
   void Write(std::string_view text)
   {
      Flush();
      fwrite(text.data(), sizeof(char), text.size(), out);
   }
};

// SAM text kept in memory: a thread of a dump on several threads formats its
// records into one, to be written in order.
struct SamText {
   std::string buf;

   void Str(std::string_view s) { buf += s; }
   void Tab() { buf += '\t'; }
   void Int(long long v) { RAMNTupleUtils::AppendInt(buf, v); }
   void EndLine() { buf += '\n'; }
};

void WriteHeader(const RAMFile &file, SamWriter &out)
{
   for (const auto &line : file.GetMetadata().header_lines) {
      out.Str(line);
      out.EndLine();
   }
}

// A column of plain numbers, read with RNTuple's bulk reads into a vector: for
// such a field that copies whole pages, where a view call per value cost about
// a tenth of a full read. The vector grows with the runs asked for.
template <typename T>
class BulkColumn {
   ROOT::RNTupleView<T> m_view;
   ROOT::RFieldBase::RBulkValues m_bulk;
   std::vector<T> m_values;

public:
   BulkColumn(ROOT::RNTupleReader &reader, const char *name)
      : m_view(reader.GetView<T>(name)), m_bulk(m_view.CreateBulk())
   {
   }
   ~BulkColumn() = default;
   // The bulk reads into m_values, so the column must stay where it was built.
   BulkColumn(const BulkColumn &) = delete;
   BulkColumn(BulkColumn &&) = delete;
   BulkColumn &operator=(const BulkColumn &) = delete;
   BulkColumn &operator=(BulkColumn &&) = delete;

   void Read(ROOT::RNTupleLocalIndex first, std::size_t rows)
   {
      if (rows > m_values.size()) {
         m_values.resize(rows);
         m_bulk.AdoptBuffer(m_values.data(), rows);
      }
      m_bulk.ReadBulk(first, nullptr, rows);
   }
   T operator[](std::size_t i) const { return m_values[i]; }
};

// Reads records straight into a record's members. A view of the whole record
// (GetView<RAMNTupleRecord>) looks up ROOT's dictionary for the class, and that
// starts ROOT's interpreter: about 0.4 s and 150 MB before the first record. The
// members are plain types and need neither. The numbers are read in bulk, a run
// of rows of one cluster at a time; the strings and vectors one row at a time.
// A run starts at 1,024 rows and doubles up to 65,536 while the rows come in
// order, so that a small region reads no more pages than before.
class RecordReader {
   static constexpr std::size_t kFirstRun = 1024;
   static constexpr std::size_t kLongestRun = 65536;

   struct Cluster {
      ROOT::NTupleSize_t first;
      ROOT::NTupleSize_t rows;
      ROOT::DescriptorId_t id;
   };

   RAMNTupleRecord m_record;
   std::vector<Cluster> m_clusters; ///< By first row.
   ROOT::NTupleSize_t m_runFirst = 0;
   std::size_t m_runRows = 0;
   std::size_t m_run = kFirstRun / 2; ///< The length of the last run asked for.

   BulkColumn<uint16_t> m_flag;
   BulkColumn<int32_t> m_refid;
   BulkColumn<int32_t> m_pos;
   BulkColumn<uint8_t> m_mapq;
   BulkColumn<int32_t> m_refnext;
   BulkColumn<int32_t> m_pnextDelta;
   BulkColumn<int32_t> m_tlenDelta;
   BulkColumn<int32_t> m_tagNm;
   BulkColumn<int32_t> m_tagAs;
   BulkColumn<int32_t> m_tagXs;
   BulkColumn<int32_t> m_tagMq;
   BulkColumn<uint32_t> m_tagRg;
   BulkColumn<uint32_t> m_tagPg;
   BulkColumn<uint32_t> m_tagOrder;
   BulkColumn<uint32_t> m_compressionFlags;
   ROOT::RNTupleView<std::vector<uint32_t>> m_cigar;
   ROOT::RNTupleView<std::string> m_seq;
   ROOT::RNTupleView<std::string> m_qual;
   ROOT::RNTupleView<std::string> m_tagMd;
   ROOT::RNTupleView<std::string> m_tagMc;
   ROOT::RNTupleView<std::string> m_tagSa;
   ROOT::RNTupleView<std::string> m_tagXa;
   ROOT::RNTupleView<std::vector<std::string>> m_tags;

   // Reads the numbers of a run of rows from \p row on, within its cluster.
   void ReadRun(ROOT::NTupleSize_t row)
   {
      m_run = row == m_runFirst + m_runRows ? std::min(m_run * 2, kLongestRun) : kFirstRun;
      const auto next = std::upper_bound(m_clusters.begin(), m_clusters.end(), row,
                                         [](ROOT::NTupleSize_t r, const Cluster &c) { return r < c.first; });
      const Cluster &cluster = *std::prev(next);
      const ROOT::NTupleSize_t inCluster = row - cluster.first;
      const auto rows = static_cast<std::size_t>(std::min<ROOT::NTupleSize_t>(m_run, cluster.rows - inCluster));
      const ROOT::RNTupleLocalIndex first(cluster.id, inCluster);
      m_flag.Read(first, rows);
      m_refid.Read(first, rows);
      m_pos.Read(first, rows);
      m_mapq.Read(first, rows);
      m_refnext.Read(first, rows);
      m_pnextDelta.Read(first, rows);
      m_tlenDelta.Read(first, rows);
      m_tagNm.Read(first, rows);
      m_tagAs.Read(first, rows);
      m_tagXs.Read(first, rows);
      m_tagMq.Read(first, rows);
      m_tagRg.Read(first, rows);
      m_tagPg.Read(first, rows);
      m_tagOrder.Read(first, rows);
      m_compressionFlags.Read(first, rows);
      m_runFirst = row;
      m_runRows = rows;
   }

   static std::vector<Cluster> Clusters(ROOT::RNTupleReader &reader)
   {
      std::vector<Cluster> clusters;
      for (const auto &c : reader.GetDescriptor().GetClusterIterable())
         clusters.push_back({c.GetFirstEntryIndex(), c.GetNEntries(), c.GetId()});
      std::sort(clusters.begin(), clusters.end(), [](const Cluster &a, const Cluster &b) { return a.first < b.first; });
      return clusters;
   }

public:
   explicit RecordReader(ROOT::RNTupleReader &reader)
      : m_clusters(Clusters(reader)),
        m_flag(reader, /*name=*/"record.flag"),
        m_refid(reader, /*name=*/"record.refid"),
        m_pos(reader, /*name=*/"record.pos"),
        m_mapq(reader, /*name=*/"record.mapq"),
        m_refnext(reader, /*name=*/"record.refnext"),
        m_pnextDelta(reader, /*name=*/"record.pnext_delta"),
        m_tlenDelta(reader, /*name=*/"record.tlen_delta"),
        m_tagNm(reader, /*name=*/"record.tag_nm"),
        m_tagAs(reader, /*name=*/"record.tag_as"),
        m_tagXs(reader, /*name=*/"record.tag_xs"),
        m_tagMq(reader, /*name=*/"record.tag_mq"),
        m_tagRg(reader, /*name=*/"record.tag_rg"),
        m_tagPg(reader, /*name=*/"record.tag_pg"),
        m_tagOrder(reader, /*name=*/"record.tag_order"),
        m_compressionFlags(reader, /*name=*/"record.compression_flags"),
        m_cigar(reader.GetView("record.cigar", &m_record.cigar)),
        m_seq(reader.GetView("record.seq", &m_record.seq)),
        m_qual(reader.GetView("record.qual", &m_record.qual)),
        m_tagMd(reader.GetView("record.tag_md", &m_record.tag_md)),
        m_tagMc(reader.GetView("record.tag_mc", &m_record.tag_mc)),
        m_tagSa(reader.GetView("record.tag_sa", &m_record.tag_sa)),
        m_tagXa(reader.GetView("record.tag_xa", &m_record.tag_xa)),
        m_tags(reader.GetView("record.tags", &m_record.tags))
   {
   }
   ~RecordReader() = default;
   // The views write into m_record, so the reader must stay where it was built.
   RecordReader(const RecordReader &) = delete;
   RecordReader(RecordReader &&) = delete;
   RecordReader &operator=(const RecordReader &) = delete;
   RecordReader &operator=(RecordReader &&) = delete;

   const RAMNTupleRecord &operator()(ROOT::NTupleSize_t row)
   {
      if (row < m_runFirst || row >= m_runFirst + m_runRows)
         ReadRun(row);
      const auto i = static_cast<std::size_t>(row - m_runFirst);
      m_record.flag = m_flag[i];
      m_record.refid = m_refid[i];
      m_record.pos = m_pos[i];
      m_record.mapq = m_mapq[i];
      m_record.refnext = m_refnext[i];
      m_record.pnext_delta = m_pnextDelta[i];
      m_record.tlen_delta = m_tlenDelta[i];
      m_record.tag_nm = m_tagNm[i];
      m_record.tag_as = m_tagAs[i];
      m_record.tag_xs = m_tagXs[i];
      m_record.tag_mq = m_tagMq[i];
      m_record.tag_rg = m_tagRg[i];
      m_record.tag_pg = m_tagPg[i];
      m_record.tag_order = m_tagOrder[i];
      m_record.compression_flags = m_compressionFlags[i];
      m_cigar(row);
      m_seq(row);
      m_qual(row);
      m_tagMd(row);
      m_tagMc(row);
      m_tagSa(row);
      m_tagXa(row);
      m_tags(row);
      return m_record;
   }
};

template <typename Out>
void WriteRecord(const RAMFile &file, std::string_view name, const RAMNTupleRecord &rec, std::string_view qual,
                 Out &out)
{
   out.Str(name);
   out.Tab();
   out.Int(rec.GetFLAG());
   out.Tab();
   out.Str(file.GetRname(rec.GetREFID()));
   out.Tab();
   out.Int(rec.GetPOS());
   out.Tab();
   out.Int(rec.GetMAPQ());
   out.Tab();
   rec.AppendCIGAR(out.buf);
   out.Tab();
   out.Str(file.GetRnext(rec.GetREFNEXT()));
   out.Tab();
   const bool sameRef = file.MateOnSameRef(rec);
   out.Int(rec.GetPNEXT(sameRef));
   out.Tab();
   out.Int(rec.GetTLEN(sameRef));
   out.Tab();
   out.Str(rec.GetSEQ());
   out.Tab();
   out.Str(qual);
   file.AppendTags(rec, out.buf);
   out.EndLine();
}

// Writes every record of \p path on \p threads threads. A reader is not shared
// between threads, so each thread opens the file for itself and formats every
// threads-th block of records (those of a block of names and qualities, 10,000
// records) into a text of its own; the texts are written in order. A thread
// never takes two consecutive blocks, so the block readers decode nothing
// ahead, which other threads would need; and the file is opened without
// RNTuple's cluster cache, which would make every thread decompress the whole
// cluster of each block it takes. At most two blocks per thread wait to be
// written.
void DumpOnThreads(const std::string &path, const std::vector<uint64_t> &ends, int threads, uint16_t requireFlags,
                   uint16_t excludeFlags, SamWriter &writer)
{
   const std::size_t stride = std::max<std::size_t>(1, std::min(static_cast<std::size_t>(threads), ends.size()));
   const std::size_t window = 2 * stride;
   std::mutex mutex;
   std::condition_variable changed;
   std::map<std::size_t, std::string> texts; // Formatted, waiting to be written.
   std::size_t written = 0;                  // Blocks taken by the writer.
   bool failed = false;

   const auto work = [&](std::size_t first) {
      try {
         ROOT::RNTupleReadOptions options;
         options.SetClusterCache(ROOT::RNTupleReadOptions::EClusterCache::kOff);
         const auto file = RAMFile::Open(path, options);
         if (!file)
            throw std::runtime_error("cannot open " + path);
         RecordReader records(file->GetReader());
         QualityBlockReader quals(*file);
         NameBlockReader names(*file);
         for (std::size_t block = first; block < ends.size(); block += stride) {
            {
               std::unique_lock<std::mutex> lock(mutex);
               changed.wait(lock, [&] { return failed || block < written + window; });
               if (failed)
                  return;
            }
            SamText text;
            for (uint64_t row = block == 0 ? 0 : ends[block - 1] + 1; row <= ends[block]; row++) {
               const auto &rec = records(row);
               const uint16_t flag = rec.GetFLAG();
               if ((flag & requireFlags) != requireFlags || (flag & excludeFlags))
                  continue;
               WriteRecord(*file, names.Get(row), rec, quals.Get(rec, row), text);
            }
            {
               const std::lock_guard<std::mutex> lock(mutex);
               texts.emplace(block, std::move(text.buf));
            }
            changed.notify_all();
         }
      } catch (...) {
         {
            const std::lock_guard<std::mutex> lock(mutex);
            failed = true;
         }
         changed.notify_all();
         throw;
      }
   };

   std::vector<std::future<void>> workers;
   workers.reserve(stride);
   for (std::size_t i = 0; i < stride; i++)
      workers.push_back(std::async(std::launch::async, work, i));
   for (std::size_t block = 0; block < ends.size(); block++) {
      std::string text;
      {
         std::unique_lock<std::mutex> lock(mutex);
         changed.wait(lock, [&] { return failed || texts.count(block) != 0; });
         if (failed)
            break;
         const auto it = texts.find(block);
         text = std::move(it->second);
         texts.erase(it);
         written++;
      }
      changed.notify_all();
      writer.Write(text);
   }
   for (auto &w : workers)
      w.get();
}

void Usage()
{
   std::cerr << "Usage: ramdump [options] <in.ram> [region]\n"
             << "Options:\n"
             << "  -h        include the header in the output\n"
             << "  -H        print the header only\n"
             << "  -c        print only the number of matching records\n"
             << "  -f INT    only records with all of these FLAG bits set\n"
             << "  -F INT    only records with none of these FLAG bits set\n"
             << "  -o FILE   write to FILE instead of stdout\n"
             << "  -@ INT    with no region, read and format on INT threads (default 1)\n"
             << "\nRegion is rname[:start[-end]], 1-based and inclusive, as in samtools.\n";
}

// Accepts 0x900 as well as 2304, as samtools does.
bool ParseFlag(const std::string &text, uint16_t &out)
{
   constexpr int kAnyBase = 0;
   char *end = nullptr;
   errno = 0;
   const long value = std::strtol(text.c_str(), &end, kAnyBase);
   if (text.empty() || *end != '\0' || errno == ERANGE || value < 0 || value > UINT16_MAX)
      return false;
   out = static_cast<uint16_t>(value);
   return true;
}

constexpr int kMaxThreads = 1024;

// A decimal number of threads and nothing else, so "3x" is an error rather than 3.
bool ParseThreads(const std::string &text, int &out)
{
   constexpr int kDecimal = 10;
   char *end = nullptr;
   errno = 0;
   const long value = std::strtol(text.c_str(), &end, kDecimal);
   if (text.empty() || *end != '\0' || errno == ERANGE || value < 1 || value > kMaxThreads)
      return false;
   out = static_cast<int>(value);
   return true;
}

} // namespace

int main(int argc, char *argv[])
{
   // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
   const std::vector<std::string> args(argv + 1, argv + argc);

   bool withHeader = false;
   bool headerOnly = false;
   bool countOnly = false;
   uint16_t requireFlags = 0;
   uint16_t excludeFlags = 0;
   int threads = 1;
   std::string outPath{};
   std::string file{};
   std::string region{};

   for (size_t i = 0; i < args.size(); i++) {
      const std::string &arg = args[i];
      if (arg == "-h") {
         withHeader = true;
      } else if (arg == "-H") {
         headerOnly = true;
      } else if (arg == "-c") {
         countOnly = true;
      } else if (arg == "-f" || arg == "-F") {
         if (++i == args.size()) {
            std::cerr << "ramdump: " << arg << " needs a value\n";
            return 1;
         }
         uint16_t value = 0;
         if (!ParseFlag(args[i], value)) {
            std::cerr << "ramdump: invalid FLAG value '" << args[i] << "'\n";
            return 1;
         }
         (arg == "-f" ? requireFlags : excludeFlags) = value;
      } else if (arg == "-o") {
         if (++i == args.size()) {
            std::cerr << "ramdump: -o needs a file name\n";
            return 1;
         }
         outPath = args[i];
      } else if (arg == "-@") {
         if (++i == args.size()) {
            std::cerr << "ramdump: -@ needs a value\n";
            return 1;
         }
         if (!ParseThreads(args[i], threads)) {
            std::cerr << "ramdump: -@ takes a number of threads from 1 to " << kMaxThreads << "\n";
            return 1;
         }
      } else if (arg == "--help") {
         Usage();
         return 0;
      } else if (!arg.empty() && arg[0] == '-') {
         std::cerr << "ramdump: unknown option '" << arg << "'\n";
         return 1;
      } else if (file.empty()) {
         file = arg;
      } else {
         region = arg;
      }
   }

   if (file.empty()) {
      Usage();
      return 1;
   }

   std::unique_ptr<FILE, int (*)(FILE *)> owned(outPath.empty() ? nullptr : std::fopen(outPath.c_str(), "w"),
                                                std::fclose);
   if (!outPath.empty() && !owned) {
      std::cerr << "ramdump: cannot write " << outPath << "\n";
      return 1;
   }
   SamWriter writer(owned ? owned.get() : stdout);

   // RNTuple's cluster cache loads every page of the columns in use for a whole
   // cluster and the next one. A region needs a few pages of each column: on
   // chr20 a region of 1 kb reads 1 MB of the file without the cache and 7 MB
   // with it, and takes less time. A whole file is read faster with the cache
   // when it is not in memory yet.
   const bool onThreads = threads > 1 && region.empty() && !countOnly && !headerOnly;
   if (onThreads)
      ROOT::EnableThreadSafety();
   ROOT::RNTupleReadOptions options;
   if (!region.empty() && region != "*")
      options.SetClusterCache(ROOT::RNTupleReadOptions::EClusterCache::kOff);
   auto ram = RAMFile::Open(file, options);
   if (!ram) {
      std::cerr << "ramdump: cannot open " << file << "\n";
      return 1;
   }

   if (headerOnly) {
      WriteHeader(*ram, writer);
      return 0;
   }

   if (withHeader && !countOnly)
      WriteHeader(*ram, writer);

   if (onThreads) {
      DumpOnThreads(file, ram->GetMetadata().block_ends, threads, requireFlags, excludeFlags, writer);
      return 0;
   }

   // Counting with no FLAG filter never has to decode a record.
   if (countOnly && requireFlags == 0 && excludeFlags == 0) {
      writer.Int(ramntuplescan(*ram, region.c_str(), nullptr));
      writer.EndLine();
      return 0;
   }

   // The views are made at the first record of the region, after the search for
   // it. RNTuple loads a cluster with every column that has a view, so views made
   // before the search would load all columns of each cluster the search visits.
   // A count needs only the flags.
   std::optional<ROOT::RNTupleView<uint16_t>> flags;
   std::optional<RecordReader> records;
   std::optional<QualityBlockReader> quals;
   std::optional<NameBlockReader> names;
   Long64_t kept = 0;

   ramntuplescan(*ram, region.c_str(), [&](Long64_t row) {
      if (countOnly) {
         if (!flags)
            flags.emplace(ram->GetReader().GetView<uint16_t>("record.flag"));
         const uint16_t flag = (*flags)(static_cast<ROOT::NTupleSize_t>(row));
         if ((flag & requireFlags) == requireFlags && !(flag & excludeFlags))
            kept++;
         return;
      }
      if (!records) {
         records.emplace(ram->GetReader());
         quals.emplace(*ram);
         names.emplace(*ram);
      }
      const auto &rec = (*records)(row);
      const uint16_t flag = rec.GetFLAG();
      if ((flag & requireFlags) != requireFlags || (flag & excludeFlags))
         return;
      WriteRecord(*ram, names->Get(static_cast<ROOT::NTupleSize_t>(row)), rec,
                  quals->Get(rec, static_cast<ROOT::NTupleSize_t>(row)), writer);
   });

   if (countOnly) {
      writer.Int(kept);
      writer.EndLine();
   }

   return 0;
}
