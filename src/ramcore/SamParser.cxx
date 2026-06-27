#include "ramcore/SamParser.h"
#include <cerrno>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace ramcore {

namespace {

bool ParseInt(const char *value, int &out, const char *field_name, size_t line_number, int min_value, int max_value)
{
    if (!value || *value == '\0') {
        std::cerr << "[SamParser] Warning: empty value for field '" << field_name << "' at line " << line_number
                  << "; record skipped.\n";
        return false;
    }
    if (std::isspace(static_cast<unsigned char>(*value))) {
        std::cerr << "[SamParser] Warning: invalid integer value '" << value << "' for field '" << field_name
                  << "' at line " << line_number << "; record skipped.\n";
        return false;
    }

    errno = 0;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);

    if (end == value || *end != '\0') {
        std::cerr << "[SamParser] Warning: invalid integer value '" << value << "' for field '" << field_name
                  << "' at line " << line_number << "; record skipped.\n";
        return false;
    }

    if (errno == ERANGE || parsed < min_value || parsed > max_value) {
        std::cerr << "[SamParser] Warning: integer value '" << value << "' for field '" << field_name
                  << "' is out of range at line " << line_number << "; record skipped.\n";
        return false;
    }

    out = static_cast<int>(parsed);
    return true;
}

} // namespace

void StripCRLF(char* str) {
    size_t len = strlen(str);
    while (len > 0 && (str[len-1] == '\n' || str[len-1] == '\r')) {
        str[--len] = '\0';
    }
}

bool SamParser::ParseFile(const char* filename, 
                         HeaderCallback header_cb, 
                         RecordCallback record_cb) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        return false;
    }
    
    lines_processed_ = 0;
    records_processed_ = 0;
    
    char line[kMaxLineLength];
    SamRecord record;
    
    while (fgets(line, kMaxLineLength, fp)) {
        lines_processed_++;
        
        if (line[0] == '@') {
            char* tab = strchr(line, '\t');
            if (tab) {
                *tab = '\0';
                StripCRLF(tab + 1);
                if (header_cb) {
                    header_cb(line, tab + 1);
                }
            } else {
                StripCRLF(line);
                if (header_cb) {
                    header_cb(line, "");
                }
            }
            continue;
        }
        
        record.Clear();
        if (ParseLine(line, record)) {
            if (record_cb) {
                record_cb(record, records_processed_);
            }
            records_processed_++;
        }
    }
    
    fclose(fp);
    return true;
}

bool SamParser::ParseLine(char* line, SamRecord& record) {
    int field_num = 0;
    char* token = strtok(line, "\t");
    
    while (token) {
        switch (field_num) {
            case 0: record.qname = token; break;
            case 1:
                if (!ParseInt(token, record.flag, "flag", lines_processed_, 0, std::numeric_limits<uint16_t>::max()))
                    return false;
                break;
            case 2: record.rname = token; break;
            case 3:
                if (!ParseInt(token, record.pos, "pos", lines_processed_, 0, std::numeric_limits<int>::max()))
                    return false;
                break;
            case 4:
                if (!ParseInt(token, record.mapq, "mapq", lines_processed_, 0, std::numeric_limits<uint8_t>::max()))
                    return false;
                break;
            case 5: record.cigar = token; break;
            case 6: record.rnext = token; break;
            case 7:
                if (!ParseInt(token, record.pnext, "pnext", lines_processed_, 0, std::numeric_limits<int>::max()))
                    return false;
                break;
            case 8:
                if (!ParseInt(token, record.tlen, "tlen", lines_processed_, std::numeric_limits<int>::min() + 1,
                              std::numeric_limits<int>::max()))
                    return false;
                break;
            case 9: record.seq = token; break;
            case 10: 
                StripCRLF(token);
                record.qual = token; 
                break;
            default:
                StripCRLF(token);
                record.optional_fields.push_back(token);
                break;
        }
        
        field_num++;
        token = strtok(nullptr, "\t");
    }
    
    return field_num >= 11; 
}

} 

