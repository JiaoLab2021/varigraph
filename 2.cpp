#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <bitset>
#include <sstream>

#include "include/haplotype_select.hpp"

using namespace std;

int main() {
    // vector<uint32_t> kmerVec = {1,2,3,2,3,2,0,0,0,3,4};

    std::string numbers = "203675 149442 149442 14479 14479 7995 7995 6948 6948 3567 3567 2369 2369 1992 1992 1792 1792 1435 1435 1119 1119 998 998 850 850 758 758 590 590 486 486 363 363 185 185 170 170 147 147 141 141 140 140 112 112 112 112 112 112 113 113 117 117 117 117 117 117 117 117 95 95 95 95 95 95 89 89 83 83 66 66 89 89 89 89 89 89 56 56 56 56 28 28 28 28 24 24 28 28 28 28 28 28 28 28 28 28 28 28 28 28";
    
    std::vector<uint32_t> kmerVec;
    std::istringstream iss(numbers);
    uint32_t num;

    while (iss >> num) {
        cerr << num << ' ';
        kmerVec.push_back(num);
    }
    cerr << endl;

    HaplotypeSelect* HaplotypeSelectClass = new HaplotypeSelect(kmerVec);
    HaplotypeSelectClass->calculate_sparsity();
    HaplotypeSelectClass->generate_sparse_frequency_vector();
    HaplotypeSelectClass->get_top_indices(15);
    for (auto it :  HaplotypeSelectClass->mTopHapVec)
    {
        cerr << it << endl;
    }
    
    delete HaplotypeSelectClass;
    
    bitset<8> bitsVector;
    bitsVector = 3;
    cerr << bitsVector[0] + bitsVector[2] << endl;
    return 0;
}