#ifndef VARIGRAPH_HPP
#define VARIGRAPH_HPP
#include <iostream>
#include <vector>
#include "zlib.h"
#include <getopt.h>

#include "construct_index.hpp"
#include "fastq_kmer.hpp"
#include "genotype.hpp"

using namespace std;

class Varigraph
{
private:
    // input file name
    const string& refFileName_;  // refgenome files
    const vector<string>& fastqFileNameVec_;  // the vector of sequencing read files
    const string& vcfFileName_;  // The name of the input VCF file.
    const string& inputFastqKmerFileName_;  // load the k-mers index of the sequencing reads from disk
    const string& inputMbfFileName_;  // Load Counting Bloom Filter index from file

    // output file name
    const string& outputFastqKmerFileName_;  // save the k-mers index of the sequences reads to disk
    const string& outputMbfFileName_;  // Save Counting Bloom Filter index to file
    const string& outputFileName_;  // save the genotyping result to disk

    const uint32_t& kmerLen_;  // the length of k-mer
    const string& prefix_;  // the prefix_ of the sequencing reads
    const uint32_t& ploidy_;  // the ploidy_ of the sequencing reads

    const uint32_t& threads_;  // the thread number

    const bool& debug_;

    ConstructIndex* ConstructIndexClassPtr_;  // Record the index of graph and reference genome

    map<string, map<uint32_t, string> > vcfInfoMap_;  // Store VCF file information

    
public:
    Varigraph(
        const string& refFileName, 
        const vector<string>& fastqFileNameVec, 
        const string& vcfFileName, 
        const string& inputFastqKmerFileName, 
        const string& inputMbfFileName, 
        const string& outputFastqKmerFileName, 
        const string& outputMbfFileName, 
        const string& outputFileName, 
        const uint32_t& kmerLen, 
        const string& prefix, 
        const uint32_t& ploidy, 
        const uint32_t& threads, 
        const bool& debug
    );
    ~Varigraph();


    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
     * @brief build the kmer index of reference and construct graph
     * 
     * @return void
    **/
    void ref_idx_construct();


    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
     * @brief build the kmer index of files
     * 
     * @return void
    **/
    void kmer_read();


    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
     * @brief genotype
     * 
     * 
     * @return void
    **/
    void genotype();
};

#endif