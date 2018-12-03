%module BloomFilter
%include "std_string.i"
%include "stdint.i"

%{
#include "BloomFilter.hpp"
%}

class BloomFilter {
public:
        ~BloomFilter();
        BloomFilter(const char * filterFilePath);
        BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize);

        void insert(const char* kmer);
        void insert(const uint64_t *hVal);
        bool insert_make_change(const uint64_t *hVal);
        bool contains(const char* kmer);

        void storeFilter(const char * filterFilePath);
            
        size_t getPop();
        size_t getSize();
};
