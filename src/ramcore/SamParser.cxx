#include "ramcore/SamParser.h"
#include <cstring>
#include <cstdlib>

namespace ramcore {

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
            case 1: record.flag = atoi(token); break;
            case 2: record.rname = token; break;
            case 3: record.pos = atoi(token); break;
            case 4: record.mapq = atoi(token); break;
            case 5: record.cigar = token; break;
            case 6: record.rnext = token; break;
            case 7: record.pnext = atoi(token); break;
            case 8: record.tlen = atoi(token); break;
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

