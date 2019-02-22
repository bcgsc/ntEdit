#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "BloomFilter.hpp"

using namespace std;

namespace opt {
size_t j = 16;
size_t k = 20;
size_t h = 9;
}

int main(int argc, const char* argv[]) {
    BloomFilter myhBF(676938880, 9, opt::k);
    ifstream hitFile(argv[1]);
    string line;
    while(getline(hitFile, line))
        myhBF.insert(line.c_str());
    
    myhBF.storeFilter("nsolidBF_k20.bf");
    
    return 0;
}
