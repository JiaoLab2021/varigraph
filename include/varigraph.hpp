#ifndef VARIGRAPH_HPP
#define VARIGRAPH_HPP
#include <iostream>
#include <vector>
#include "zlib.h"
#include <iomanip>
#include <locale>
#include <getopt.h>
#include <tuple>
#include <filesystem>

#include "construct_index.hpp"
#include "fastq_kmer.hpp"
#include "genotype.hpp"

using namespace std;

class Varigraph
{
protected:
    // input file name
    const string& refFileName_;  // refgenome files
    const string& vcfFileName_;  // The name of the input VCF file.

    // Sample configuration file
    const string& samplesConfigFileName_;  // The name of the input sample configuration file.
    vector<tuple<string, vector<string> > > sampleConfigTupleVec_;  // Store sample configuration information. vector<tuple<sampleName, vector<readPath> > >

    const string& inputGraphFileName_;  // load the Genome Geaph from disk

    // output file name
    const string& outputGraphFileName_;  // save the Genome Geaph to disk

    // fast mode
    const bool& fastMode_;

    uint32_t& kmerLen_;  // the length of k-mer
    const string& sampleType_;  // specify the genotype of the reference genome (hom/het)
    const uint32_t& samplePloidy_;  // genome ploidy
    uint32_t& vcfPloidy_;  // ploidy of genotypes in VCF file

    const uint32_t& haploidNum_;  // the haploid number for genotyping
    const uint32_t& chrLenThread_;  // Chromosome granularity
    const string& transitionProType_;  // transition probability type
    const bool& svGenotypeBool_;  // structural variation genotyping only

    const bool& debug_;

    const uint32_t& threads_;  // the thread number

    const float& minSupportingReads_;  // the minimum number of supporting reads for a variant

    const bool& useUniqueKmers_;  // use only unique k-mers for indexing
    const bool& useDepth_;  // use sequencing depth as the depth for homozygous k-mers

    ConstructIndex* ConstructIndexClassPtr_ = nullptr;  // Record the index of graph and reference genome

    map<string, map<uint32_t, string> > vcfInfoMap_;  // Store VCF file information

    float ReadDepth_ = 0.0f;  // sequencing data depth
    float hapKmerCoverage_ = 0.0f;  // haplotype k-mer coverage
public:
    Varigraph(
        const string& refFileName, 
        const string& vcfFileName, 
        const string& samplesConfigFileName, 
        const string& inputGraphFileName, 
        const string& outputGraphFileName, 
        const bool& fastMode, 
        uint32_t& kmerLen, 
        const string& sampleType, 
        const uint32_t& samplePloidy, 
        uint32_t& vcfPloidy, 
        const uint32_t& haploidNum, 
        const uint32_t& chrLenThread, 
        const string& transitionProType, 
        const bool& svGenotypeBool, 
        const bool& debug, 
        const uint32_t& threads, 
        const float& minSupportingReads, 
        const bool& useUniqueKmers, 
        const bool& useDepth
    );
    ~Varigraph() {
        if (ConstructIndexClassPtr_ != nullptr) {
            delete ConstructIndexClassPtr_;
            ConstructIndexClassPtr_ = nullptr;
        }
    }


    /**
     * @author zezhen du
     * @date 2024/01/04
     * @version v1.0.1
     * @brief Construct genome graph from reference genome and variants
     * 
     * @return void
    **/
    void construct();

    /**
     * @author zezhen du
     * @date 2024/01/04
     * @version v1.0.1
     * @brief Load genome graph
     * 
     * @return void
    **/
    void load();

    /**
     * @brief Parse sample configuration file (sampleName, readPath1, readPath2)
     * 
     * @return void
    */
    void parse_sample_config();

    /**
     * @brief fastq and genotype
     * 
     * @return void
    */
    void fastq_genotype();

    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
     * @brief build the kmer index of files
     * 
     * @param fastqFileNameVec    Sequencing data
     * 
     * @return void
    **/
    void kmer_read(vector<string> fastqFileNameVec);

    /**
     * @author zezhen du
     * @date 2024/01/19
     * @version v1.0.1
     * @brief calculate the average coverage of k-mers
     * 
     * @return void
    **/
    void cal_ave_cov_kmer();

    /**
     * @author zezhen du
     * @date 2024/01/22
     * @version v1.0.1
     * @brief Get homozygous k-mer
     * 
     * @return kmerCovFreMap   map<coverage, frequency>
    **/
    map<uint8_t, uint64_t> get_hom_kmer();

    /**
     * @author zezhen du
     * @date 2024/01/22
     * @version v1.0.1
     * @brief Looking for the k-mer peak
     * 
     * @param kmerCovFreMap   map<coverage, frequence>
     * 
     * @return maxCoverage, homCoverage
    **/
    tuple<uint8_t, uint8_t> get_hom_kmer_c(const map<uint8_t, uint64_t>& kmerCovFreMap);

    /**
     * @author zezhen du
     * @date 2024/01/22
     * @version v1.0.1
     * @brief calculate haplotype kmer coverage
     * 
     * @param homCoverage
     * 
     * @return void
    **/
    void cal_hap_kmer_cov(const uint8_t& homCoverage);

    /**
     * @author zezhen du
     * @date 2024/01/22
     * @version v1.0.1
     * @brief Print histogram
     * 
     * @param maxCoverage
     * @param homCoverage
     * @param kmerCovFreMap   map<coverage, frequence>
     * 
     * @return void
    **/
    void kmer_histogram(
        const uint8_t& maxCoverage, 
        const uint8_t& homCoverage, 
        const map<uint8_t, uint64_t>& kmerCovFreMap
    );


    /**
     * @author zezhen du
     * @date 2023/06/27
     * @version v1.0.1
     * @brief genotype
     * 
     * @param sampleName
     * 
     * @return void
    **/
    void genotype(string sampleName);
};

#endif