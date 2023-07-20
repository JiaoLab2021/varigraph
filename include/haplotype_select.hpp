#ifndef HAPLOTYPE_SELECT_HPP
#define HAPLOTYPE_SELECT_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <queue>
#include <utility>

using namespace std;

struct ComparePair {
    bool operator()(const std::pair<double, uint16_t>& a, const std::pair<double, uint16_t>& b) {
        return a.first > b.first; // descending order
    }
};

class HaplotypeSelect {
private:
    std::random_device rd;
    std::mt19937 prng;

    std::vector<uint32_t> kmerCountVec_;

    // Sparsity parameter
    double sparsity_;

    // sample frequency
    vector<double> freqVec_;

    std::vector<uint32_t> find_diff(const std::vector<uint32_t>& vec1, const std::vector<uint32_t>& vec2);

    void get_indices(std::vector<uint16_t>& kmerNonzeroIdx);

    static bool compare_pair(const std::pair<double, uint16_t>& a, const std::pair<double, uint16_t>& b);

public:
    std::vector<uint16_t> mTopHapVec;

    HaplotypeSelect(const std::vector<uint32_t>& kmerCountVec);

    void calculate_sparsity();

    void generate_sparse_frequency_vector();

    void get_top_indices(int numIndices);
};

#endif