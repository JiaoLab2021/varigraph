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


/**
 * @author zezhen du
 * @date 2024/09/24
 * @brief Parameter configuration structure
 * 
 * @return void
**/
class VarigraphConfig {
public:
    string refFileName;               // r: refgenome files
    string vcfFileName;               // v: The name of the input VCF file.
    string samplesConfigFileName;     // s: The name of the input sample configuration file.
    string inputGraphFileName;        // 1: load the Genome Geaph from disk
    string outputGraphFileName;       // 1: save the Genome Geaph to disk
    bool fastMode;                    // 3: fast mode
    uint32_t kmerLen;                 // k: the length of k-mer
    string sampleType;                // g: specify the genotype of the reference genome (hom/het)
    uint32_t samplePloidy;            // 2: genome ploidy
    uint32_t vcfPloidy;               // 2: ploidy of genotypes in VCF file
    uint32_t haploidNum;              // n: the haploid number for genotyping
    uint32_t chrLenThread;            // 3: Chromosome granularity
    string transitionProType;         // m: transition probability type
    bool svGenotypeBool;              // 4: structural variation genotyping only
    bool debug;                       // D: Debug code
    uint32_t threads;                 // t: thread
    float minSupportingGQ;            // 5: minimum site quality (GQ) value for genotype
    bool useUniqueKmers;              // 4: use only unique k-mers for indexing
    bool useDepth;                    // 6: use sequencing depth as the depth for homozygous k-mers

    // Default constructor, setting default values
    VarigraphConfig() 
        : refFileName(""), 
          vcfFileName(""), 
          samplesConfigFileName(""), 
          inputGraphFileName("graph.bin"), 
          outputGraphFileName("graph.bin"),
          fastMode(false), 
          kmerLen(27), 
          sampleType("het"), 
          samplePloidy(2), 
          vcfPloidy(2), 
          haploidNum(15), 
          chrLenThread(1 * 1e6), 
          transitionProType("rec"), 
          svGenotypeBool(false), 
          debug(false), 
          threads(10), 
          minSupportingGQ(0.0), 
          useUniqueKmers(false), 
          useDepth(false) {}

    // Construct log output
    void logConstructionConfig() const {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Number of threads: " << threads << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "k-mer size: " << kmerLen << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Reference file path: " << refFileName << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Variants file path: " << vcfFileName << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Ploidy of genotypes in the VCF file: " << vcfPloidy << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Debug mode: " << (debug ? "Enabled" : "Disabled") << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Fast mode: " << (fastMode ? "Enabled" : "Disabled") << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Use only unique k-mers for indexing: " << (useUniqueKmers ? "Enabled" : "Disabled") << endl << endl << endl;
    }

    // Genotype log output
    void logGenotypeConfig() const {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Number of threads: " << threads << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Genome graph file: " << inputGraphFileName << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Sample configuration file: " << samplesConfigFileName << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Sample genome status: " << sampleType << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Sample ploidy: " << samplePloidy << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Number of haploids for genotyping: " << haploidNum << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Chromosome granularity: " << chrLenThread << " bp" << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Transition probability type: " << transitionProType << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Structural variation genotyping only: " << (svGenotypeBool ? "Enabled" : "Disabled") << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Minimum site quality (GQ): " << minSupportingGQ << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Use sequencing depth as the depth for homozygous k-mers: " << (useDepth ? "Enabled" : "Disabled") << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Debug mode: " << (debug ? "Enabled" : "Disabled") << endl << endl << endl;
    }
};


class Varigraph
{
protected:
    // input file name
    const string refFileName_;  // refgenome files
    const string vcfFileName_;  // The name of the input VCF file.

    // Sample configuration file
    const string samplesConfigFileName_;  // The name of the input sample configuration file.
    vector<tuple<string, vector<string> > > sampleConfigTupleVec_;  // Store sample configuration information. vector<tuple<sampleName, vector<readPath> > >

    const string& inputGraphFileName_;  // load the Genome Geaph from disk

    // output file name
    const string& outputGraphFileName_;  // save the Genome Geaph to disk

    // fast mode
    const bool fastMode_;

    uint32_t kmerLen_;  // the length of k-mer
    const string sampleType_;  // specify the genotype of the reference genome (hom/het)
    const uint32_t samplePloidy_;  // genome ploidy
    uint32_t vcfPloidy_;  // ploidy of genotypes in VCF file

    const uint32_t haploidNum_;  // the haploid number for genotyping
    const uint32_t chrLenThread_;  // Chromosome granularity
    const string transitionProType_;  // transition probability type
    const bool svGenotypeBool_;  // structural variation genotyping only

    const bool debug_;

    const uint32_t threads_;  // the thread number

    const float minSupportingGQ_;  // minimum site quality (GQ) value for genotype

    const bool useUniqueKmers_;  // use only unique k-mers for indexing
    const bool useDepth_;  // use sequencing depth as the depth for homozygous k-mers

    ConstructIndex* ConstructIndexClassPtr_ = nullptr;  // Record the index of graph and reference genome

    map<string, map<uint32_t, string> > vcfInfoMap_;  // Store VCF file information

    float ReadDepth_ = 0.0f;  // sequencing data depth
    float hapKmerCoverage_ = 0.0f;  // haplotype k-mer coverage
public:
    Varigraph(const VarigraphConfig& config) 
        : refFileName_(config.refFileName), 
          vcfFileName_(config.vcfFileName), 
          samplesConfigFileName_(config.samplesConfigFileName), 
          inputGraphFileName_(config.inputGraphFileName), 
          outputGraphFileName_(config.outputGraphFileName), 
          fastMode_(config.fastMode), 
          kmerLen_(config.kmerLen), 
          sampleType_(config.sampleType), 
          samplePloidy_(config.samplePloidy), 
          vcfPloidy_(config.vcfPloidy), 
          haploidNum_(config.haploidNum), 
          chrLenThread_(config.chrLenThread), 
          transitionProType_(config.transitionProType), 
          svGenotypeBool_(config.svGenotypeBool), 
          debug_(config.debug), 
          threads_(config.threads), 
          minSupportingGQ_(config.minSupportingGQ), 
          useUniqueKmers_(config.useUniqueKmers), 
          useDepth_(config.useDepth) {
        // cerr
        std::cerr.imbue(std::locale(""));  // Thousandth output
    }

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