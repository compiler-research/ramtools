// ramstats: statistics of a RAM file, printed the way samtools prints them.
// `ramstats idxstats file.ram` prints what `samtools idxstats` prints for the
// same reads, from counts the converter stored, without reading a record.

#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

void Usage()
{
   std::cerr << "Usage: ramstats idxstats <in.ram>\n"
             << "\nidxstats prints a line per reference, as samtools idxstats does: the name,\n"
             << "the length, the mapped and the unmapped records, then '*' and the records\n"
             << "with no reference.\n";
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

} // namespace

int main(int argc, char *argv[])
{
   // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
   const std::vector<std::string> args(argv + 1, argv + argc);
   if (args.size() == 1 && args[0] == "--help") {
      Usage();
      return 0;
   }
   if (args.size() != 2 || args[0] != "idxstats") {
      Usage();
      return 1;
   }

   const auto file = RAMFile::Open(args[1]);
   if (!file) {
      std::cerr << "ramstats: cannot open " << args[1] << "\n";
      return 1;
   }

   // References without an @SQ line, which samtools would not read, get length 0.
   const RAMFileMetadata &meta = file->GetMetadata();
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
   std::cout << out;
   return 0;
}
