#ifndef ramcore_SAMPARSER_H
#define ramcore_SAMPARSER_H

#include <string>
#include <vector>
#include <functional>
#include <cstdio>

namespace ramcore {

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
    
    void Clear() {
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

class SamParser {
public:
    using HeaderCallback = std::function<void(const std::string& tag, const std::string& content)>;
    using RecordCallback = std::function<void(const SamRecord& record, size_t record_number)>;
    
    bool ParseFile(const char* filename, 
                   HeaderCallback header_cb, 
                   RecordCallback record_cb);
    
    size_t GetLinesProcessed() const { return lines_processed_; }
    size_t GetRecordsProcessed() const { return records_processed_; }
    
private:
    size_t lines_processed_ = 0;
    size_t records_processed_ = 0;
    
    bool ParseLine(char* line, SamRecord& record);
    static constexpr int kMaxLineLength = 10240;
};

void StripCRLF(char* str);

} 

#endif 

