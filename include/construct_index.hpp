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
#include <iomanip>

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


// k-mer iterator and haplotype information
struct kmerCovFreBitVec {
    uint8_t c = 0;  // coverage: read
    uint8_t f = 0;  // frequency: genome or graph

    vector<int8_t> BitVec;  // haplotype bitmap: 0000 0000, Each bits represents a haplotype, 0->False 1->True

    kmerCovFreBitVec() {}

    void clear() {
        c = 0;
        f = 0;
        vector<int8_t>().swap(BitVec);
    }

    // overloaded assignment operator
    kmerCovFreBitVec& operator=(const kmerCovFreBitVec& other) {
        if (this != &other) {
            c = other.c;
            f = other.f;
            BitVec = other.BitVec;
        }
        return *this;
    }
};


// Forward and backward algorithm scoring
struct HMMScore {
    long double a = 0;  // store alpha scores of haplotype combination
    long double b = 0;  // Store betas of haplotype combination

    vector<uint16_t> hapVec;  // haplotype combination information
};


// posterior
struct posteriorStr
{
    long double probability;  // Posterior probability, GPP
    vector<uint16_t> hapVec;  // haplotype information
    vector<uint64_t> kmerNumVec;  // The number of k-mers corresponding to the haplotype, NAK
    vector<float> kmerAveCovVec;  // The k-mer average depth corresponding to the haplotype, CAK

    uint8_t uniqueKmerNum;  // Number of node-unique k-mers, UK

    posteriorStr() : probability(0.0), uniqueKmerNum(0) {}
};


struct nodeSrt {
    vector<string> seqVec;  // Store sequence information of query: vector<seq>

    vector<uint16_t> hapGtVec;  // Genotype information for each haplotype at this locus: vector<GT>

    vector<uint64_t> kmerHashVec;  // All k-mer hashes in the node, vector<kmerHash>
    vector<unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator> GraphKmerHashHapStrMapIterVec;  // Iterator pointing to mGraphKmerHashHapStrMapP, vector<iter>
    
    vector<HMMScore> HMMScoreVec;  // store all alpha and betas scores of nodes

    posteriorStr posteriorInfo;  // maximum posterior probability
};


class ConstructIndex
{
private:
    // input file names
    string refFileName_;
    string vcfFileName_;

    const string& inputGraphFileName_;  // load the Genome Geaph from disk
    const string& outputGraphFileName_;  // save the Genome Geaph to disk

    bool fastMode_;

    uint32_t threads_;

public:
    map<string, map<uint32_t, nodeSrt> > mGraphMap;  // Store the sequence and k-mer information of the graph: map<chr, map<nodeStart, nodeSrt> >
    unordered_map<uint64_t, kmerCovFreBitVec> mGraphKmerHashHapStrMap;  // Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFreBitVec>

    map<uint16_t, string> mHapMap;  // Store haplotype information: map<hapIdx, hapName>
    uint16_t mHapNum;

    string mVcfHead;  // Store VCF file comment lines
    map<string, map<uint32_t, vector<string> > > mVcfInfoMap;  // Store VCF file information, map<chromosome, map<start, vector> >
    
    uint32_t mKmerLen;
    uint32_t mVcfPloidy;

    // number
    uint32_t mSnpNum = 0;
    uint32_t mIndelNum = 0;
    uint32_t mInsNum = 0;
    uint32_t mDelNum = 0;
    uint32_t mInvNum = 0;
    uint32_t mDupNum = 0;
    uint32_t mOtherNum = 0;

    unordered_map<string, string> mFastaSeqMap;  // map<chromosome, sequence>
    unordered_map<string, uint32_t> mFastaLenMap;  // map<chromosome, sequence length>
    uint64_t mGenomeSize = 0;  // Reference genome size

    BloomFilter* mbf;  // The Counting Bloom filter of reference genome's k-mers informations

    unordered_map<uint16_t, tuple<uint16_t, uint16_t> > mHapIdxQRmap;  // map<hapIdx, tuple<quotient, remainder> >

    uint64_t mGraphBaseNum = 0;  // The number of bases in the graph

    ConstructIndex(
        const string& refFileName, 
        const string& vcfFileName, 
        const string& inputGraphFileName, 
        const string& outputGraphFileName, 
        const bool& fastMode, 
        const uint32_t& kmerLen, 
        const uint32_t& vcfPloidy, 
        const bool& debug, 
        const uint32_t& threads
    );
    ~ConstructIndex();

    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
	 * @brief Free memory
     * 
     * @return void
	**/
    void clear_mbf();


    /**
     * @author zezhen du
     * @date 2023/07/21
     * @version v1.0.1
	 * @brief Free memory
     * 
     * @return void
	**/
    void clear_mGraphKmerCovFreMap();

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
     * @date 2024/01/04
     * @version v1.0
     * @brief make mHapIdxQRmap
     * 
     * @return void
    **/
    void make_QRmap();


    /**
     * @author zezhen du
     * @date 2023/09/12
     * @version v1.0
	 * @brief Build the index of vcf
     * 
     * @param chromosome
     * @param refStart
     * @param refLen
     * @param lineVec     vcf split list
     * @param qrySeqVec   the sequence Vector of ALT
     * @param gtIndex     Index where gt resides
     * 
     * @return void
	**/
    void vcf_construct(
        const string& chromosome, 
        const uint32_t& refStart, 
        const uint32_t& refLen, 
        const vector<string>& lineVec, 
        const vector<string>& qrySeqVec, 
        const int& gtIndex
    );


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
     * @date 2023/08/30
     * @version v1.0
     * @brief Merge k-mer information from Genome Graph into nodes.
     * 
     * @return void
    **/
    void graph2node();


    /**
     * @author zezhen du
     * @brief Save the graph index to file
     * 
     * @return void
    **/
    void save_index();

    /**
     * @author zezhen du
     * @brief load the graph index from file
     * 
     * @return void
    **/
    void load_index() ;
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
     * @brief graph index for k-mer (threads)
     * 
     * @date 2023/09/01
     * 
     * @param chromosome            mGraphMap output by construct??map<chr, map<start, nodeSrt> >
     * @param nodeIter              node iterator
     * @param startNodeMap          Chromosome all nodes
     * @param fastMode              fast mode
     * @param kmerLen               the length of kmer
     * @param bf                    Kmer frequency in the reference genome: Counting Bloom Filter
     * @param vcfPloidy             ploidy of genotypes in VCF file
     * @param hapIdxQRmap           map<hapIdx, tuple<quotient, remainder> >
     * 
     * @return {nodeIter, tmpKmerHapBitMap, kmerHashFreMap}     kmer: map<kmerHash, vector<int8_t> >
    **/
    tuple<map<uint32_t, nodeSrt>::iterator, unordered_map<uint64_t, vector<int8_t> >, map<uint64_t, uint8_t> > index_run(
        string chromosome, 
        map<uint32_t, nodeSrt>::iterator nodeIter, 
        const map<uint32_t, nodeSrt>& startNodeMap, 
        const bool& fastMode, 
        const uint32_t& kmerLen, 
        BloomFilter* bf, 
        const uint32_t& vcfPloidy, 
        const unordered_map<uint16_t, tuple<uint16_t, uint16_t> >& hapIdxQRmap
    );


    /**
     * @author zezhen du
     * @date 2023/08/13
     * @version v1.0.1
     * @brief find the sequence information corresponding to the haplotype upstream and downstream of the node
     * 
     * @param haplotype       haplotype index
     * @param altGt           genotype
     * @param altSeq          ALT sequence
     * @param seqLen          k-mer length - 1
     * @param nodeIter        startNodeMap iterator where the node is located
     * @param startNodeMap    All chromosome node information startNodeMap, index function construction
     * 
     * @return pair<upSeq, downSeq>
    **/
    pair<string, string> find_node_up_down_seq(
        const uint16_t& haplotype, 
        const uint16_t& altGt, 
        string& altSeq, 
        const uint32_t & seqLen,
        const map<uint32_t, nodeSrt>::iterator & nodeIter, 
        const map<uint32_t, nodeSrt> & startNodeMap
    );


    // sort by k-mer frequency in ascending order
    bool compare_frequency(
        const unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator& a, 
        const unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator& b
    );


    /**
     * @brief Merge k-mer information from Genome Graph into nodes. (threads)
     * 
     * @date 2023/08/30
     * 
     * @param kmerHashVec                     Node k-mers hash
     * @param GraphKmerHashHapStrMapIterVec   Iterator pointing to mGraphKmerHashHapStrMap, vector<iter>
     * @param GraphKmerHashHapStrMap          Total k-mers coverage and frequency information
     * 
     * @return 0
    **/
    int graph2node_run(
        vector<uint64_t>& kmerHashVec, 
        vector<unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator>& GraphKmerHashHapStrMapIterVec, 
        const unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap
    );
}

#endif