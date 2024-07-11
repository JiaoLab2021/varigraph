#ifndef FASTQ_KMER_CUH
#define FASTQ_KMER_CUH

#include "fastq_kmer.hpp"

#include "kmer.cuh"

using namespace std;

class FastqKmerKernel : public FastqKmer {
public:
    int buffer_ = 500;  // Buffer size

    // Constructor
    FastqKmerKernel(
        unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap, 
        const vector<string>& fastqFileNameVec, 
        const uint32_t& kmerLen, 
        const uint32_t & threads, 
        const int buffer
    ) : FastqKmer(GraphKmerHashHapStrMap, fastqFileNameVec, kmerLen, threads) {
        buffer_ = buffer;
    }

    // Destructor
    ~FastqKmerKernel() {}

    /**
     * @author zezhen du
     * @date 2024/04/30
     * @version v1.0.1
     * @brief building the kmer index of sequencing read
     * 
     * @return void
    **/
    void build_fastq_index_kernel();

    /**
     * @author zezhen du
     * @date 2024/04/30
     * @version v1.0
     * @brief building the kmer index of sequencing read
     * 
     * @param fastqFileName       sequencing read
     * 
     * @return void
    **/
    void fastq_file_open_kernel(const string & fastqFileName);
};

namespace fastq_kmer_kernel {
    /**
     * This function adds coverage to a map of k-mers.
     * 
     * @param kmerHashVec The hash vector of the k-mer.
     * @param kmerCovVec The coverage vector of the k-mer.
     * @param GraphKmerHashHapStrMap The map of k-mers to their coverage and frequency bit vectors.
     * @param startIndex The starting index of the k-mer.
     * @param endIndex The ending index of the k-mer.
    */
    void add_cov_to_map(
        const vector<uint64_t>& kmerHashVec, 
        const vector<uint8_t>& kmerCovVec,
        unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap,
        uint32_t startIndex,
        uint32_t endIndex
    );
}


#endif