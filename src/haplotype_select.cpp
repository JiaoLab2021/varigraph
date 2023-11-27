#include "../include/haplotype_select.hpp"

HaplotypeSelect::HaplotypeSelect(const std::vector<uint32_t>& kmerCountVec) : kmerCountVec_(kmerCountVec), prng(rd()) {}


std::vector<uint32_t> HaplotypeSelect::find_diff(const std::vector<uint32_t>& vec1, const std::vector<uint32_t>& vec2) {
    std::vector<uint32_t> diff;
    std::unordered_set<uint32_t> set2(vec2.begin(), vec2.end());

    for (const auto& element : vec1) {
        if (set2.find(element) == set2.end()) {
            diff.push_back(element);
        }
    }

    return diff;
}

void HaplotypeSelect::get_indices(std::vector<uint16_t>& kmerNonzeroIdx) {
    int preScore = -1;
    for (uint16_t i = 0; i < kmerNonzeroIdx.size(); i++) {
        if (kmerNonzeroIdx[i] != preScore) {
            preScore = kmerNonzeroIdx[i];
            kmerNonzeroIdx[i] = std::distance(kmerCountVec_.begin(), std::find(kmerCountVec_.begin(), kmerCountVec_.end(), kmerNonzeroIdx[i]));
        } else {
            preScore = kmerNonzeroIdx[i];
            kmerNonzeroIdx[i] = std::distance(kmerCountVec_.begin(), std::find(kmerCountVec_.begin() + kmerNonzeroIdx[i - 1] + 1, kmerCountVec_.end(), kmerNonzeroIdx[i]));
        }
    }
}

bool HaplotypeSelect::compare_pair(const std::pair<double, uint16_t>& a, const std::pair<double, uint16_t>& b) {
    return a.first > b.first;
}

void HaplotypeSelect::calculate_sparsity() {
    uint64_t nonZeroCount = 0;
    uint64_t vectorSize = kmerCountVec_.size();

    for (uint64_t i = 0; i < vectorSize; i++) {
        if (kmerCountVec_[i] != 0) {
            nonZeroCount++;
        }
    }

    sparsity_ = vectorSize > 0 && nonZeroCount > 0 ? static_cast<double>(nonZeroCount) / vectorSize : 1;
}

void HaplotypeSelect::generate_sparse_frequency_vector() {
    // select non-zero indices
    uint16_t hapNum = kmerCountVec_.size();
    uint16_t nonzeroNum = std::min(static_cast<uint16_t>(hapNum * sparsity_), static_cast<uint16_t>(std::count_if(kmerCountVec_.begin(), kmerCountVec_.end(), [](uint32_t num){return num != 0;})));
    std::vector<uint16_t> kmerNonzeroIdx;
    for (uint16_t i = 0; i < kmerCountVec_.size(); ++i) {
        if (kmerCountVec_[i] != 0) {
            kmerNonzeroIdx.push_back(i);
        }
    }
    
    // randomly assign frequencies to zero indices
    freqVec_.resize(hapNum, 0.0);

    // update frequency vector
    double freqDirichletSum = 0;
    for (const auto& it : kmerNonzeroIdx) {
        freqVec_[it] = std::gamma_distribution<double>(kmerCountVec_[it] + 1.0, 1)(prng);
        freqDirichletSum += freqVec_[it];
    }

    if (freqDirichletSum > 0) {
        for (auto& it : freqVec_) {
            it /= freqDirichletSum;
        }
    }
}


unordered_map<uint16_t, double> HaplotypeSelect::get_top_indices(int numIndices) {
    std::priority_queue<std::pair<double, uint16_t>, std::vector<std::pair<double, uint16_t> >, ComparePair> pq;  // descending order
    double sum = 0.0;  // For summation

    for (uint16_t i = 0; i < freqVec_.size(); i++) {
        pq.push(std::make_pair(freqVec_[i], i));
        sum += freqVec_[i];

        if (pq.size() > numIndices) {
            sum -= pq.top().first;
            pq.pop();
        }
    }

    unordered_map<uint16_t, double> hapIdxScoreMap;  // map<hapIdx, possibility>
    while (!pq.empty()) {
        mTopHapVec.push_back(pq.top().second);
        hapIdxScoreMap[pq.top().second] = pq.top().first / sum;  // normalization
        pq.pop();
    }

    return hapIdxScoreMap;
}