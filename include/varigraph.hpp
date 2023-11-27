#ifndef VARIGRAPH_HPP
#define VARIGRAPH_HPP
#include <iostream>
#include <vector>
#include "zlib.h"
#include <iomanip>
#include <locale>
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
    const string& inputMbfFileName_;  // Load Counting Bloom Filter index from file
    const string& inputGraphFileName_;  // load the Genome Geaph from disk
    const string& inputFastqKmerFileName_;  // load the k-mers index of the sequencing reads from disk

    // output file name
    const string& outputMbfFileName_;  // Save Counting Bloom Filter index to file
    const string& outputGraphFileName_;  // save the Genome Geaph to disk
    const string& outputFastqKmerFileName_;  // save the k-mers index of the sequences reads to disk
    const string& outputFileName_;  // save the genotyping result to disk

    // fast mode
    const bool& fastMode_;

    const uint32_t& kmerLen_;  // the length of k-mer
    const string& sampleName_;  // the sampleName_ of the sequencing reads
    const string& genomeType_;  // specify the genotype of the reference genome (homozygous/heterozygous)
    const uint32_t& refPloidy_;  // genome ploidy
    const uint32_t& vcfPloidy_;  // ploidy of genotypes in VCF file

    const uint32_t& haploidNum_;  // the haploid number for genotyping 

    const bool& debug_;

    const uint32_t& threads_;  // the thread number

    ConstructIndex* ConstructIndexClassPtr_;  // Record the index of graph and reference genome

    map<string, map<uint32_t, string> > vcfInfoMap_;  // Store VCF file information

    float ReadDepth_ = 0.0f;  // sequencing data depth
public:
    Varigraph(
        const string& refFileName, 
        const vector<string>& fastqFileNameVec, 
        const string& vcfFileName, 
        const string& inputMbfFileName, 
        const string& inputGraphFileName, 
        const string& inputFastqKmerFileName, 
        const string& outputMbfFileName, 
        const string& outputGraphFileName, 
        const string& outputFastqKmerFileName, 
        const string& outputFileName, 
        const bool& fastMode, 
        const uint32_t& kmerLen, 
        const string& sampleName, 
        const string& genomeType, 
        const uint32_t& refPloidy, 
        const uint32_t& vcfPloidy, 
        const uint32_t& haploidNum, 
        const bool& debug, 
        const uint32_t& threads
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