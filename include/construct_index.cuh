#ifndef CONSTRUCT_INDEX_CUH
#define CONSTRUCT_INDEX_CUH

#include "construct_index.hpp"
#include "cuda_error_handling.hpp"

#include "kmer.cuh"
#include "kmer.hpp"
#include "counting_bloom_filter.cuh"


using namespace std;


// Define the ConstructIndexKernel class
class ConstructIndexKernel : public ConstructIndex {
public:
    BloomFilterKernel* mbfD = nullptr;  // Counting Bloom Filter on GPU

    int buffer_ = 100;  // Buffer size

    // Constructor
    ConstructIndexKernel(
        const string& refFileName, 
        const string& vcfFileName, 
        const string& inputGraphFileName, 
        const string& outputGraphFileName, 
        const bool& fastMode, 
        const bool& useUniqueKmers, 
        const uint32_t& kmerLen, 
        const uint32_t& vcfPloidy, 
        const bool& debug, 
        const uint32_t& threads, 
        const int buffer
    );

    // Destructor
    ~ConstructIndexKernel() {
        if (mbfD != nullptr) {
            delete mbfD;
            mbfD = nullptr;
        }
    }


    /**
     * @author zezhen du
     * @date 2024/04/29
     * @version v1.0.1
     * @brief Making Counting Bloom Filter
     * 
     * @return void
    **/
    void make_mbf_kernel();

    /**
     * @author zezhen du
     * @date 2024/04/28
     * @version v1.0
     * @brief building the k-mer index of graph on GPU
     * 
     * @return void
    **/
    void index_kernel();
};

/**
 * @brief graph index for k-mer on GPU
 * 
 * @date 2024/04/29
 * 
 * @param chromosome            mGraphMap output by construct  map<chr, map<start, nodeSrt> >
 * @param nodeIter              node iterator
 * @param startNodeMap          Chromosome all nodes
 * @param fastMode              fast mode
 * @param useUniqueKmers        use unique kmers for indexing
 * @param kmerLen               the length of kmer
 * @param bfD                   Counting Bloom Filter on GPU
 * @param vcfPloidy             ploidy of genotypes in VCF file
 * @param hapIdxQRmap           map<hapIdx, tuple<quotient, remainder> >
 * 
 * @return {nodeIter, tmpKmerHapBitMap, kmerHashFreMap}     kmer: map<kmerHash, vector<int8_t> >
**/
tuple<map<uint32_t, nodeSrt>::iterator, unordered_map<uint64_t, vector<int8_t> >, map<uint64_t, uint8_t> > index_run_kernel(
    string chromosome, 
    map<uint32_t, nodeSrt>::iterator nodeIter, 
    const map<uint32_t, nodeSrt>& startNodeMap, 
    const bool& fastMode, 
    const bool& useUniqueKmers,
    const uint32_t& kmerLen, 
    BloomFilterKernel* bfD, 
    const uint32_t& vcfPloidy, 
    const unordered_map<uint16_t, tuple<uint16_t, uint16_t> >& hapIdxQRmap
);

#endif