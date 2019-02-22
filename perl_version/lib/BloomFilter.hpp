/*
 *
 * BloomFilter.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_

#include "iostream"
#include "fstream"
#include <cstring>
#include "nthash.hpp"

using namespace std;

inline unsigned popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
    & 0xf;
}

class BloomFilter {
public:
    
#pragma pack(push)
#pragma pack(1) //to maintain consistent values across platforms
    struct FileHeader {
        char magic[8];
        uint32_t hlen;
        uint64_t size;
        uint32_t nhash;
        uint32_t kmer;
        double dFPR;
        uint64_t nEntry;
        uint64_t tEntry;
    };
#pragma pack(pop)

    BloomFilter(const char * fPath) {
        FileHeader header;
        ifstream myFile(fPath, ios::in | ios::binary);
        myFile.seekg(0, ios::beg);
        myFile.read((char *)(&header), sizeof(struct FileHeader));
        m_size = header.size;
        m_hashNum = header.nhash;
        m_kmerSize = header.kmer;
        
        char magic[9];
        strncpy(magic, header.magic, 8);
        magic[8] = '\0';
        cerr << magic << "\tsize(bits): " << m_size << "\thash: " << m_hashNum << "\t k: " << m_kmerSize << endl;

        m_filter = new unsigned char [(m_size + 7)/8]();
        myFile.read((char *)m_filter, (m_size + 7)/8);
        myFile.close();
    }
    
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize):
        m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_filter = new unsigned char [(m_size + 7)/8]();
    }

    bool insert_make_change(const uint64_t *hVal) {
        bool change=false;
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            unsigned char old_byte = __sync_fetch_and_or(&m_filter[hLoc/8],(1<<(7-hLoc%8)));
            if((old_byte&(1<<(7-hLoc%8)))==0) change=true;
        }
        return change;
    }

    void insert(const uint64_t *hVal) {
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

    void insert(const char* kmer) {
	uint64_t hVal = NTC64(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
	    uint64_t mhVal=NTE64(hVal, m_kmerSize, i);
            size_t hLoc = mhVal % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

    bool contains(const char* kmer) const {
        uint64_t hVal = NTC64(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal=NTE64(hVal, m_kmerSize, i);
            size_t hLoc = mhVal % m_size;
            if ((m_filter[hLoc / 8] & (1 << (7 - hLoc % 8))) == 0)
                return false;
        }
        return true;
    }

    bool contains(const uint64_t *hVal) {
	    for (unsigned i=0; i< m_hashNum; i++) {
		    size_t hLoc = hVal[i] % m_size; 
		    if ((m_filter[hLoc / 8] & (1 << (7 - hLoc % 8))) == 0)
			    return false; 
	    }
	    return true; 
    }

    void storeFilter(const char * fPath) const {
        FileHeader header;
        strncpy(header.magic, "BlOOMFXX", 8);
        char magic[9];
        strncpy(magic, header.magic, 8);
        magic[8] = '\0';

        header.hlen = sizeof(struct FileHeader);
        header.size = m_size;
        header.nhash = m_hashNum;
        header.kmer = m_kmerSize;
        header.dFPR = pow(getPop()*1.0/m_size, m_hashNum);
        header.nEntry = 0;
        header.tEntry = 0;
        
        
        ofstream myFile(fPath, ios::out | ios::binary);
        myFile.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));
        myFile.write(reinterpret_cast<char*>(m_filter), (m_size + 7)/8);
        myFile.close();
    }

    void printBloomFilterDetails() {
	    std::cout << "BLOOM::\tsize: " << m_size << "\tnumber hash functions: " << m_hashNum << "\tkmer size: " << m_kmerSize 
		    << "\tFPR: " << pow(getPop()*1.0/m_size, m_hashNum) << std::endl; 
    }

    size_t getPop() const {
        size_t i, popBF=0;
        #pragma omp parallel for reduction(+:popBF)
        for(i=0; i<(m_size + 7)/8; i++)
            popBF = popBF + popCnt(m_filter[i]);
        return popBF;
    }
	
    size_t getSize() const {
        return m_size;
    }

    unsigned getHashNum() const {
	    return m_hashNum; 
    }

    unsigned getKmerSize() const {
	    return m_kmerSize; 
    }
     
    ~BloomFilter() {
        delete[] m_filter;
    }

private:
    BloomFilter(const BloomFilter& that);
    unsigned char * m_filter;
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
};

#endif /* BLOOMFILTER_H_ */
