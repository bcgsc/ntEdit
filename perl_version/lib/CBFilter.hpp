/*
 *
 * CBFilter.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#ifndef CBFILTER_H_
#define CBFILTER_H_
#include <string>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include "nthash.hpp"

using namespace std;


class CBFilter {
public:

    CBFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, unsigned repCap):
    m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_reCap(repCap) {
        m_filter = new unsigned char [m_size]();
    }

    bool insert_and_test(const uint64_t *hVal) {
        bool greaterFlag = true;
        unsigned minCount = 256;
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            if(m_filter[hLoc] < m_reCap) {
                if(m_filter[hLoc] < minCount)
                    minCount = m_filter[hLoc];
                greaterFlag = false;
            }
        }
        if(!greaterFlag) {
            for (unsigned i = 0; i < m_hashNum; i++) {
                size_t hLoc = hVal[i] % m_size;
                if(m_filter[hLoc] == minCount) {
                    #pragma omp atomic
                    ++m_filter[hLoc];
                }
            }
        }
        return greaterFlag;
    }

    ~CBFilter() {
        delete[] m_filter;
    }

private:
    CBFilter(const CBFilter& that); //to prevent copy construction
    unsigned char *m_filter;
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
    unsigned m_reCap;
};

#endif /* CBFILTER_H_ */
