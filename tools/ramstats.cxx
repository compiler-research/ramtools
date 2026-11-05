// ramstats: statistics of a RAM file, printed the way samtools prints them.
// `ramstats idxstats file.ram` prints what `samtools idxstats` prints for the
// same reads, from counts the converter stored, without reading a record.
// `ramstats flagstat file.ram` prints what `samtools flagstat` prints, reading
// only the flag, refid, refnext and mapq columns.

#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

void Usage()
{
   std::cerr << "Usage: ramstats idxstats <in.ram>\n"
             << "       ramstats flagstat [-O default|tsv|json] [--scan] <in.ram>\n"
             << "\nidxstats prints a line per reference, as samtools idxstats does: the name,\n"
             << "the length, the mapped and the unmapped records, then '*' and the records\n"
             << "with no reference.\n"
             << "\nflagstat prints the counts of the records by their FLAG bits, QC-passed and\n"
             << "QC-failed apart, in samtools flagstat's text, TSV or JSON. The converters\n"
             << "count them while writing; --scan counts them again from the columns.\n";
}

// LN of every @SQ line by SN, as written in the header.
std::map<std::string, std::string> SequenceLengths(const std::vector<std::string> &header)
{
   std::map<std::string, std::string> lengths;
   for (const auto &line : header) {
      if (line.rfind("@SQ\t", 0) != 0)
         continue;
      std::string name{};
      std::string length = "0";
      size_t start = 4;
      while (start < line.size()) {
         size_t end = line.find('\t', start);
         if (end == std::string::npos)
            end = line.size();
         if (line.compare(start, 3, "SN:") == 0)
            name = line.substr(start + 3, end - start - 3);
         else if (line.compare(start, 3, "LN:") == 0)
            length = line.substr(start + 3, end - start - 3);
         start = end + 1;
      }
      lengths.emplace(name, length);
   }
   return lengths;
}

uint64_t At(const std::vector<uint64_t> &v, size_t i)
{
   return i < v.size() ? v[i] : 0;
}

std::string Idxstats(const RAMFile &file)
{
   // References without an @SQ line, which samtools would not read, get length 0.
   const RAMFileMetadata &meta = file.GetMetadata();
   const auto lengths = SequenceLengths(meta.header_lines);
   const RAMRefCounts &counts = meta.counts;
   const std::vector<std::string> &names = meta.rname_refs;
   std::string out{};
   for (size_t i = 0; i < names.size(); i++) {
      const auto length = lengths.find(names[i]);
      out += names[i];
      out += '\t';
      out += length != lengths.end() ? length->second : "0";
      out += '\t';
      out += std::to_string(At(counts.mapped, i));
      out += '\t';
      out += std::to_string(At(counts.unmapped, i));
      out += '\n';
   }
   out += "*\t0\t0\t";
   out += std::to_string(counts.unplaced);
   out += '\n';
   return out;
}

// A line of the output, in samtools' order.
struct FlagLine {
   uint64_t RAMFlagStat::*count;
   const char *text;         ///< The label in the default and TSV outputs.
   const char *json;         ///< The key in the JSON output.
   uint64_t RAMFlagStat::*of; ///< What the percentage is of, or null for none.
};

const std::array<FlagLine, 16> kFlagLines = {{
   {&RAMFlagStat::total, "total (QC-passed reads + QC-failed reads)", "total", nullptr},
   {&RAMFlagStat::primary, "primary", "primary", nullptr},
   {&RAMFlagStat::secondary, "secondary", "secondary", nullptr},
   {&RAMFlagStat::supplementary, "supplementary", "supplementary", nullptr},
   {&RAMFlagStat::duplicates, "duplicates", "duplicates", nullptr},
   {&RAMFlagStat::primary_duplicates, "primary duplicates", "primary duplicates", nullptr},
   {&RAMFlagStat::mapped, "mapped", "mapped", &RAMFlagStat::total},
   {&RAMFlagStat::primary_mapped, "primary mapped", "primary mapped", &RAMFlagStat::primary},
   {&RAMFlagStat::paired, "paired in sequencing", "paired in sequencing", nullptr},
   {&RAMFlagStat::read1, "read1", "read1", nullptr},
   {&RAMFlagStat::read2, "read2", "read2", nullptr},
   {&RAMFlagStat::properly_paired, "properly paired", "properly paired", &RAMFlagStat::paired},
   {&RAMFlagStat::both_mapped, "with itself and mate mapped", "with itself and mate mapped", nullptr},
   {&RAMFlagStat::singletons, "singletons", "singletons", &RAMFlagStat::paired},
   {&RAMFlagStat::mate_other_chr, "with mate mapped to a different chr", "with mate mapped to a different chr",
    nullptr},
   {&RAMFlagStat::mate_other_chr_q5, "with mate mapped to a different chr (mapQ>=5)",
    "with mate mapped to a different chr (mapQ >= 5)", nullptr},
}};

// refnext indexes rnext_refs, which holds "=" besides names, and samtools
// compares the mate's reference id with the read's. This maps each refnext to
// a refid, kSameRef for "=", and -1 for a name the header lacks, as htslib does.
constexpr int32_t kSameRef = -2;

std::vector<int32_t> MateRefIds(const RAMFile &file)
{
   std::vector<int32_t> ids;
   for (const auto &name : file.GetMetadata().rnext_refs)
      ids.push_back(name == "=" ? kSameRef : file.FindRefId(name));
   return ids;
}

// The columns are read this many rows at a time.
constexpr ROOT::NTupleSize_t kChunk = 65536;

// Counts every record from the columns, as the converters count them into
// METADATA (flagstat --scan). The four columns are read a chunk at a time with
// RNTuple's bulk reads, which copy whole pages into the vectors: reading them
// one value per view call took four fifths of the time.
RAMFlagCounts CountFlags(RAMFile &file)
{
   ROOT::RNTupleReader &reader = file.GetReader();
   auto flagView = reader.GetView<uint16_t>("record.flag");
   auto refidView = reader.GetView<int32_t>("record.refid");
   auto refnextView = reader.GetView<int32_t>("record.refnext");
   auto mapqView = reader.GetView<uint8_t>("record.mapq");
   std::vector<uint16_t> flags(kChunk);
   std::vector<int32_t> refids(kChunk);
   std::vector<int32_t> refnexts(kChunk);
   std::vector<uint8_t> mapqs(kChunk);
   auto flagBulk = flagView.CreateBulk();
   auto refidBulk = refidView.CreateBulk();
   auto refnextBulk = refnextView.CreateBulk();
   auto mapqBulk = mapqView.CreateBulk();
   flagBulk.AdoptBuffer(flags.data(), kChunk);
   refidBulk.AdoptBuffer(refids.data(), kChunk);
   refnextBulk.AdoptBuffer(refnexts.data(), kChunk);
   mapqBulk.AdoptBuffer(mapqs.data(), kChunk);
   const std::vector<int32_t> mateRefIds = MateRefIds(file);

   // Bulk reads take an index inside a cluster. The reader's copy of the
   // descriptor may be replaced while reading, so the clusters are listed first.
   std::vector<std::pair<ROOT::DescriptorId_t, ROOT::NTupleSize_t>> clusters;
   for (const auto &cluster : reader.GetDescriptor().GetClusterIterable())
      clusters.emplace_back(cluster.GetId(), cluster.GetNEntries());

   RAMFlagCounts counts;
   for (const auto &[cluster, rows] : clusters) {
      for (ROOT::NTupleSize_t start = 0; start < rows; start += kChunk) {
         const std::size_t n = std::min(kChunk, rows - start);
         const ROOT::RNTupleLocalIndex first(cluster, start);
         flagBulk.ReadBulk(first, nullptr, n);
         refidBulk.ReadBulk(first, nullptr, n);
         refnextBulk.ReadBulk(first, nullptr, n);
         mapqBulk.ReadBulk(first, nullptr, n);
         for (std::size_t i = 0; i < n; i++) {
            const int32_t refid = refids[i];
            const int32_t refnext = refnexts[i];
            int32_t mate = -1;
            if (refnext >= 0 && static_cast<std::size_t>(refnext) < mateRefIds.size())
               mate = mateRefIds[static_cast<std::size_t>(refnext)] == kSameRef
                         ? refid
                         : mateRefIds[static_cast<std::size_t>(refnext)];
            counts.Note(flags[i], refid, mate, mapqs[i]);
         }
      }
   }
   return counts;
}

// samtools divides in float and prints two decimals; with nothing to divide
// by it prints N/A, or null in JSON.
std::string Percent(uint64_t n, uint64_t of, bool json)
{
   if (of == 0)
      return json ? "null" : "N/A";
   std::ostringstream s;
   s << std::fixed << std::setprecision(2)
     << static_cast<double>(static_cast<float>(n) / static_cast<float>(of)) * 100.0;
   if (!json)
      s << '%';
   return s.str();
}

std::string FlagstatDefault(const std::array<RAMFlagStat, 2> &c)
{
   std::string out{};
   for (const FlagLine &line : kFlagLines) {
      out += std::to_string(c[0].*line.count) + " + " + std::to_string(c[1].*line.count) + " ";
      if (line.count == &RAMFlagStat::total)
         out += "in ";
      out += line.text;
      if (line.of != nullptr)
         out += " (" + Percent(c[0].*line.count, c[0].*line.of, /*json=*/false) + " : " +
                Percent(c[1].*line.count, c[1].*line.of, /*json=*/false) + ")";
      out += '\n';
   }
   return out;
}

std::string FlagstatTsv(const std::array<RAMFlagStat, 2> &c)
{
   std::string out{};
   for (const FlagLine &line : kFlagLines) {
      out += std::to_string(c[0].*line.count) + "\t" + std::to_string(c[1].*line.count) + "\t" + line.text + "\n";
      if (line.of != nullptr)
         out += Percent(c[0].*line.count, c[0].*line.of, /*json=*/false) + "\t" +
                Percent(c[1].*line.count, c[1].*line.of, /*json=*/false) + "\t" + line.text + " %\n";
   }
   return out;
}

// Spaces before the newlines included, as samtools prints them.
std::string FlagstatJson(const std::array<RAMFlagStat, 2> &c)
{
   const std::array<const char *, 2> names = {"QC-passed reads", "QC-failed reads"};
   std::string out = "{\n";
   for (std::size_t w = 0; w < 2; w++) {
      out += std::string(" \"") + names.at(w) + "\": { \n";
      for (std::size_t i = 0; i < kFlagLines.size(); i++) {
         const FlagLine &line = kFlagLines.at(i);
         out += std::string("  \"") + line.json + "\": " + std::to_string(c.at(w).*line.count);
         if (line.of != nullptr)
            out += std::string(", \n  \"") + line.json +
                   " %\": " + Percent(c.at(w).*line.count, c.at(w).*line.of, /*json=*/true);
         out += i + 1 < kFlagLines.size() ? ", \n" : " \n";
      }
      out += w == 0 ? " },\n" : " }\n";
   }
   out += "}\n";
   return out;
}

// ramstats flagstat [-O default|tsv|json] [--scan] <in.ram>; the arguments after "flagstat".
int Flagstat(const std::vector<std::string> &args)
{
   std::string format = "default";
   std::string path{};
   bool scan = false;
   for (std::size_t i = 0; i < args.size(); i++) {
      if ((args[i] == "-O" || args[i] == "--output-fmt") && i + 1 < args.size()) {
         format = args[++i];
      } else if (args[i] == "--scan") {
         scan = true;
      } else if (path.empty() && !args[i].empty() && args[i][0] != '-') {
         path = args[i];
      } else {
         Usage();
         return 1;
      }
   }
   if (path.empty() ||
       (format != "default" && format != "tsv" && format != "TSV" && format != "json" && format != "JSON")) {
      Usage();
      return 1;
   }

   const auto file = RAMFile::Open(path);
   if (!file) {
      std::cerr << "ramstats: cannot open " << path << "\n";
      return 1;
   }
   const std::array<RAMFlagStat, 2> counts = scan ? CountFlags(*file).stat : file->GetMetadata().flags.stat;
   if (format == "json" || format == "JSON")
      std::cout << FlagstatJson(counts);
   else if (format == "tsv" || format == "TSV")
      std::cout << FlagstatTsv(counts);
   else
      std::cout << FlagstatDefault(counts);
   return 0;
}

int Run(const std::vector<std::string> &args)
{
   if (args.size() == 1 && args[0] == "--help") {
      Usage();
      return 0;
   }
   if (!args.empty() && args[0] == "flagstat")
      return Flagstat(std::vector<std::string>(args.begin() + 1, args.end()));
   if (args.size() != 2 || args[0] != "idxstats") {
      Usage();
      return 1;
   }

   const auto file = RAMFile::Open(args[1]);
   if (!file) {
      std::cerr << "ramstats: cannot open " << args[1] << "\n";
      return 1;
   }
   std::cout << Idxstats(*file);
   return 0;
}

} // namespace

int main(int argc, char *argv[])
{
   // A damaged file makes RNTuple throw while the records are read.
   try {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      return Run(std::vector<std::string>(argv + 1, argv + argc));
   } catch (const std::exception &e) {
      std::cerr << "ramstats: " << e.what() << "\n";
   }
   return 1;
}
