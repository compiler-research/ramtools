// ramdump: write a RAM file back out as SAM, the way `samtools view` does.
// `ramdump -h file.ram` should reproduce the SAM the file was built from, and
// `ramdump -c region` should agree with `samtools view -c region`.

#include "ramcore/QualityBlocks.h"
#include "ramcore/RAMNTupleView.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <Rtypes.h>
#include <TFile.h>
#include <TList.h>
#include <TNamed.h>

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
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

   void Str(const std::string &s) { buf += s; }
   void Tab() { buf += '\t'; }
   void Int(long long v) { buf += std::to_string(v); }

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
};

// The converters store each header line as a TNamed: the tag ("@SQ") as the
// name and everything after the first tab as the title.
void WriteHeader(const std::string &file, SamWriter &out)
{
   std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
   if (!f || f->IsZombie())
      return;

   auto *headers = f->Get<TList>("headers");
   if (!headers)
      return;

   TIter next(headers);
   while (auto *obj = next()) {
      auto *named = dynamic_cast<TNamed *>(obj);
      if (!named)
         continue;
      out.Str(named->GetName());
      const std::string content = named->GetTitle();
      if (!content.empty()) {
         out.Tab();
         out.Str(content);
      }
      out.EndLine();
   }
}

void WriteRecord(const RAMNTupleRecord &rec, const std::string &qual, SamWriter &out)
{
   out.Str(rec.GetQNAME());
   out.Tab();
   out.Int(rec.GetFLAG());
   out.Tab();
   out.Str(rec.GetRNAME());
   out.Tab();
   out.Int(rec.GetPOS());
   out.Tab();
   out.Int(rec.GetMAPQ());
   out.Tab();
   out.Str(rec.GetCIGAR());
   out.Tab();
   out.Str(rec.GetRNEXT());
   out.Tab();
   out.Int(rec.GetPNEXT());
   out.Tab();
   out.Int(rec.GetTLEN());
   out.Tab();
   out.Str(rec.GetSEQ());
   out.Tab();
   out.Str(qual);

   for (const auto &tag : rec.GetTags()) {
      out.Tab();
      out.Str(tag);
   }
   out.EndLine();
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

   if (headerOnly) {
      WriteHeader(file, writer);
      return 0;
   }

   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      std::cerr << "ramdump: cannot open " << file << "\n";
      return 1;
   }

   if (withHeader && !countOnly)
      WriteHeader(file, writer);

   // Counting with no FLAG filter never has to decode a record.
   if (countOnly && requireFlags == 0 && excludeFlags == 0) {
      writer.Int(ramntuplescan(*reader, region.c_str(), nullptr));
      writer.EndLine();
      return 0;
   }

   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   Long64_t kept = 0;

   ramntuplescan(*reader, region.c_str(), [&](Long64_t row) {
      const auto &rec = view(row);
      const uint16_t flag = rec.GetFLAG();
      if ((flag & requireFlags) != requireFlags || (flag & excludeFlags))
         return;
      kept++;
      if (!countOnly)
         WriteRecord(rec, quals.Get(rec, static_cast<ROOT::NTupleSize_t>(row)), writer);
   });

   if (countOnly) {
      writer.Int(kept);
      writer.EndLine();
   }

   return 0;
}
