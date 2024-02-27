#ifndef FASTQ_KMER_HPP
#define FASTQ_KMER_HPP
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <regex>
#include <malloc.h>
#include <getopt.h>
#include <numeric>


#include "construct_index.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "ThreadPool.hpp"
#include "kmer.hpp"
#include "GzChunkReader.hpp"
#include "save.hpp"


using namespace std;


class FastqKmer
{
private:
    vector<string> fastqFileNameVec_;
    uint32_t kmerLen_;
    uint32_t threads_;

    std::mutex mtx;

    unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap_;  // // Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFreBitVecP>

public:
    uint64_t mReadBase = 0;  //Sequencing file size

    /**
     * @author zezhen du
     * @date 2023/07/17
     * @version v1.0
	 * @brief building the kmer index of sequencing read
     * 
     * @param GraphKmerCovFreMap     Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFreBitVecP>
     * @param fastqFileNameVec       the vector of sequencing read
     * @param kmerLen                the length of k-mer
     * @param threads                threads
     * 
     * @return void
	**/
    FastqKmer(
        unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap, 
        const vector<string>& fastqFileNameVec, 
        const uint32_t& kmerLen, 
        const uint32_t & threads
    );
    ~FastqKmer() {};

    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0
	 * @brief building the kmer index of sequencing read
     * 
     * @return void
	**/
    void build_fastq_index();

    /**
     * @author zezhen du
     * @date 2023/09/06
     * @version v1.0
	 * @brief building the kmer index of sequencing read
     * 
     * @param fastqFileName       sequencing read
     * 
     * @return void
	**/
    void fastq_file_open(
        const string & fastqFileName
    );


    /**
     * @author zezhen du
     * @date 2023/07/18
     * @version v1.0
	 * @brief Save the kmer statistical results of fastq
     * 
     * @param outputFileName
     * 
     * @return void
	**/
    void save_index(
        const string & outputFileName
    );


    /**
     * @author zezhen du
     * @date 2023/07/18
     * @version v1.0.1
     * @brief load the k-mer index from disk
     * 
     * @param inputFileName
     * 
     * @return void
    **/
    void load_index(
        const string & inputFileName
    );
};


namespace fastq_kmer
{
    /**
     * @author zezhen du
     * @date 2023/08/01
     * @version v1.0.1
     * @brief Component index multithreaded functions to reduce thread locks
     * 
     * @param sequenceVec                sequence vector
     * @param kmerLen                    The length of kmer
     * @param GraphKmerHashHapStrMap     Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFreBitVec>
     * 
     * @return hashVec               vector<uint64_t>
    **/
    vector<uint64_t> fastq_file_open_run(
        vector<string> sequenceVec, 
        uint32_t kmerLen, 
        const unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMapP
    );
}

#endif