#ifndef CONSTRUCT_INDEX_HPP
#define CONSTRUCT_INDEX_HPP
#include <fstream>
#include <string>
#include <getopt.h>
#include <iostream>
#include <vector>
#include "zlib.h"
#include <map>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <malloc.h>
#include <iterator>
#include <sstream>

#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "ThreadPool.hpp"
#include "counting_bloom_filter.hpp"
#include "GzChunkReader.hpp"

#include "kseq.h"


// If the macro is not defined, define it
#ifndef DIVIDE_BY_8
    // Get the macro definition of the quotient divided by 8
    #define DIVIDE_BY_8(x) ((x) >> 3)
#endif

// If the macro is not defined, define it
#ifndef GET_LOW_3_BITS
    // Extract the macro definition of the lower 3 bits
    #define GET_LOW_3_BITS(x) ((x) & 0x7)
#endif


using namespace std;


// kmer coverage and frequency structure
struct kmerCovFre
{
    uint8_t c = 0;  // coverage: read
    uint8_t f = 0;  // frequency: genome or graph
};

// k-mer iterator and haplotype information
struct kmerHashHap
{
    uint64_t kmerHash;  // Hash value of k-mer
    vector<int8_t> BitVec;  // haplotype bitmap: 0000 0000, Each bits represents a haplotype, 0->False 1->True

    kmerHashHap() : kmerHash(0), BitVec() {}

    kmerHashHap(uint64_t hash, const vector<int8_t>& bitVector) : kmerHash(hash), BitVec(bitVector) {}

    void clear() {
        kmerHash = 0;
        vector<int8_t>().swap(BitVec);
    }
};


// Forward and backward algorithm scoring
struct HMMScore
{
    long double a = 0;  // store alpha scores of haplotype combination
    long double b = 0;  // Store betas of haplotype combination

    vector<uint16_t> hapVec;  // haplotype combination information
};


struct nodeSrt
{
    unordered_map<uint16_t, string> seqMap;  // Store sequence information of query: map<gt, seq>

    vector<uint16_t> hapGtVec;  // Genotype information for each haplotype at this locus: vector<GT>

    vector<kmerHashHap> kmerHashHapVec;  // Whether the storage haplotype contains the corresponding kmer: map<kmerHash, vector<int8_t> >:  0000 0000, Each bits represents a haplotype, 0->False 1->True

    vector<HMMScore> HMMScoreVec;  // store all alpha and betas scores of nodes

    tuple<long double, vector<uint16_t> > posteriorTup{0.0, vector<uint16_t>()};  // maximum posterior probability
};


class ConstructIndex
{
private:
    // input file names
    string refFileName_;
    string vcfFileName_;

    const string& inputMbfFileName_;  // Load Counting Bloom Filter index from file
    const string& outputMbfFileName_;  // Save Counting Bloom Filter index to file

    uint32_t kmerLen_;

    string prefix_;

    uint32_t ploidy_;

    uint32_t threads_;

    bool debug_;

    uint64_t genomeSize_ = 0;

public:
    map<string, map<uint32_t, nodeSrt> > mGraphMap;  // Store the sequence and k-mer information of the graph: map<chr, map<nodeStart, nodeSrt> >
    unordered_map<uint64_t, kmerCovFre> mGraphKmerCovFreMap;  // Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFre>
    map<uint16_t, string> mHapMap;  // Store haplotype information: map<hapIdx, hapName>
    uint16_t mHapNum;

    string mVcfHead;  // Store VCF file comment lines
    map<string, map<uint32_t, string> > mVcfInfoMap;  // Store VCF file information

    unordered_map<string, string> mFastaMap;  // map<chromosome, sequence>
    BloomFilter* mbf;  // The Counting Bloom filter of reference genome's k-mers informations

    ConstructIndex(
        const string& refFileName, 
        const string& vcfFileName, 
        const string& inputMbfFileName, 
        const string& outputMbfFileName, 
        const uint32_t& kmerLen, 
        const string& prefix, 
        const uint32_t& ploidy, 
        const uint32_t& threads, 
        const bool& debug
    );
    ~ConstructIndex();

    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
     * @brief Building the k-mers index of reference genome
     * 
     * @return void
	**/
    void build_fasta_index();

    /**
     * @author zezhen du
     * @date 2023/07/18
     * @version v1.0.1
     * @brief Making Counting Bloom Filter
     * 
     * @return void
    **/
    void make_mbf();

    /**
     * @author zezhen du
     * @date 2022/12/16
     * @version v1.0
	 * @brief Graph construction
     * 
     * 
     * @return void
	**/
    void construct();


    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0
	 * @brief building the k-mer index of graph
     * 
     * @return void
	**/
    void index();


    /**
     * @author zezhen du
     * @date 2023/07/14
     * @version v1.0
     * @brief Remove duplicate K-mers from graph genome
     * 
     * @return void
    **/
    void kmer_deduplication();

    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
	 * @brief Free memory
     * 
     * @return void
	**/
    void clear_memory();
};


namespace construct_index
{
    /**
     * @author zezhen du
     * @date 2023/06/11
     * @version v1.0
	 * @brief split the GT information
     * 
     * @param gtTxt     GT information
     * 
     * @return vector<string> gtVecTmp
	**/
    vector<string> gt_split(
        const string & gtTxt
    );


    /**
     * @author zezhen du
     * @date 2023/06/29
     * @version v1.0.1
	 * @brief Sets a specific bit of a given number to 1.
     * 
     * @param num         The value that needs to be modified.
     * @param bitIndex    The index of the bit to be set to 1.
     * 
     * @return void
	**/
    template <typename T>
    void set_bit_to_one(T& num, int bitIndex);

    /**
     * @author zezhen du
     * @date 2023/06/29
     * @version v1.0.1
	 * @brief Sets a specific bit of a given number to 0.
     * 
     * @param num         The value that needs to be modified.
     * @param bitIndex    The index of the bit to be set to 0.
     * 
     * @return void
	**/
    template <typename T>
    void set_bit_to_zero(T& num, int bitIndex);

    /**
     * @author zezhen du
     * @date 2023/06/29
     * @version v1.0.1
	 * @brief Retrieves the value of a specific bit in a given number.
     * 
     * @param num         The value from which to retrieve the bit.
     * @param bitIndex    The index of the bit to be queried.
     * 
     * @return int
	**/
    template <typename T>
    int get_bit(T num, int bitIndex);


    /**
     * @author zezhen du
     * @date 2023/07/03
     * @version v1.0.1
	 * @brief print bitmap.
     * 
     * @param num         The value from which to print.
     * 
     * @return int
	**/
    string print_bits(int8_t num);


    /**
     * @brief graph index for kmer (threads)
     * 
     * @date 2023/07/14
     * 
     * @param chromosome            mGraphMap output by construct，map<chr, map<start, nodeSrt> >
     * @param nodeIter              node iterator
     * @param startNodeMap          Chromosome all nodes
     * @param kmerLen               the length of kmer
     * @param bf                    Kmer frequency in the reference genome: Counting Bloom Filter
     * @param ploidy                Ploidy of vcf file
     * @param debug                 debug code
     * 
     * @return {nodeIter, tmpKmerHapBitMap}     kmer: map<kmerHash, vector<int8_t> >
    **/
    tuple<map<uint32_t, nodeSrt>::iterator, vector<kmerHashHap> > index_run(
        string chromosome, 
        map<uint32_t, nodeSrt>::iterator nodeIter, 
        const map<uint32_t, nodeSrt>& startNodeMap, 
        const uint32_t& kmerLen, 
        BloomFilter* bf, 
        const uint32_t& ploidy, 
        const bool& debug
    );


    /**
     * @author zezhen du
     * @date 2023/07/13
     * @version v1.0.1
     * @brief find the sequence information corresponding to the haplotype upstream and downstream of the node
     * 
     * @param haplotype       单倍型索引
     * @param seqLen          kmer的长度
     * @param nodeIter        startNodeMap中节点所在的迭代器
     * @param startNodeMap    染色体所有的节点信息startNodeMap，index函数构造
     * 
     * @return pair<upSeq, downSeq>
    **/
    pair<string, string> find_node_up_down_seq(
        const uint16_t & haplotype, 
        const uint32_t & seqLen,
        const map<uint32_t, nodeSrt>::iterator & nodeIter, 
        const map<uint32_t, nodeSrt> & startNodeMap
    );


    /**
     * @brief Remove duplicate K-mers from graph genome (threads)
     * 
     * @date 2023/07/13
     * 
     * @param kmerHashHapVec            Node k-mers information
     * @param GraphKmerCovFreMap        Total k-mers coverage and frequency information
     * 
     * @return tuple<vector<kmerHashHap>&, vector<kmerHashHap> >: tuple<ref(kmerHashHapVec), kmerHashHapVecTmp>
    **/
    tuple<vector<kmerHashHap>&, vector<kmerHashHap> > kmer_deduplication_run(
        vector<kmerHashHap>& kmerHashHapVec, 
        unordered_map<uint64_t, kmerCovFre>& GraphKmerCovFreMap
    );


    /**
     * @author zezhen du
     * @date 2023/07/13
     * @version v1.0.1
     * @brief Convert map to vector
     * 
     * @param inputMap            unordered_map<uint64_t, vector<int8_t> >
     * @param tmpkmerHashHapVec   Converted vector: vector<kmerHashHap>
     * 
     * @return void
    **/
    void convert_map2vec(const unordered_map<uint64_t, vector<int8_t> >& inputMap, vector<kmerHashHap>& tmpkmerHashHapVec);
}

#endif