#ifndef RAMCORE_SAMPARSER_H
#define RAMCORE_SAMPARSER_H

#include <string>
#include <vector>
#include <functional>
#include <cstdio>

namespace ramcore {

/**
 * \struct SamRecord
 * \brief Plain-old-data container representing one SAM alignment record.
 *
 * Mirrors the eleven mandatory SAM columns (`QNAME`, `FLAG`, `RNAME`, `POS`, `MAPQ`,
 * `CIGAR`, `RNEXT`, `PNEXT`, `TLEN`, `SEQ`, `QUAL`) and stores any auxiliary
 * fields verbatim in `optional_fields`.
 *
 * The struct is intentionally kept trivially-constructible so that it can be
 * allocated on the stack and cleared via `Clear()` without expensive heap
 * operations inside tight parsing loops.
 *
 * \sa SamParser
 */
struct SamRecord {
   std::string qname;
   int flag = 0;
   std::string rname;
   int pos = 0;
   int mapq = 0;
   std::string cigar;
   std::string rnext;
   int pnext = 0;
   int tlen = 0;
   std::string seq;
   std::string qual;
   std::vector<std::string> optional_fields;

   void Clear()
   {
      qname.clear();
      flag = 0;
      rname.clear();
      pos = 0;
      mapq = 0;
      cigar.clear();
      rnext.clear();
      pnext = 0;
      tlen = 0;
      seq.clear();
      qual.clear();
      optional_fields.clear();
   }
};

/**
 * \class SamParser
 * \brief Streaming parser for SAM text files.
 *
 * Reads a SAM file line-by-line, dispatching callbacks for header directives
 * and alignment records.  Parsing is zero-copy: the input buffer is split in
 * place and converted into a `SamRecord` which is then passed to the user
 * callback.
 *
 * Example usage:
 * \code
 * ramcore::SamParser parser;
 * parser.ParseFile("reads.sam",
 *     [](const std::string &tag, const std::string &content) {
 *         // handle @HD / @SQ / @PG ...
 *     },
 *     [](const ramcore::SamRecord &rec, size_t idx) {
 *         // process alignment #idx
 *     });
 * \endcode
 *
 * \sa SamRecord
 */
class SamParser {
public:
   using HeaderCallback = std::function<void(const std::string &tag, const std::string &content)>;
   using RecordCallback = std::function<void(const SamRecord &record, size_t record_number)>;

   bool ParseFile(const char *filename, HeaderCallback header_cb, RecordCallback record_cb);

   size_t GetLinesProcessed() const { return lines_processed_; }
   size_t GetRecordsProcessed() const { return records_processed_; }

private:
   size_t lines_processed_ = 0;
   size_t records_processed_ = 0;

   bool ParseLine(char *line, SamRecord &record);
   static constexpr int kMaxLineLength = 10240;
};

} // namespace ramcore

#endif // RAMCORE_SAMPARSER_H
