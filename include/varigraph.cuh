#ifndef VARIGRAPH_CUH
#define VARIGRAPH_CUH

#include "construct_index.cuh"
#include "fastq_kmer.cuh"

#include "varigraph.hpp"

using namespace std;


/**
 * @author zezhen du
 * @date 2024/09/24
 * @brief Parameter configuration structure
 * 
 * @return void
**/
class VarigraphKernelConfig : public VarigraphConfig {
public:
    int gpu;      // GPU ID
    int buffer;   // Buffer size

    // Default constructor, setting default values
    VarigraphKernelConfig() 
        : VarigraphConfig(),
          gpu(0), 
          buffer(100) {}

    // Construct log output
    void logConstructionConfigKernel() const {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Number of threads: " << threads << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "k-mer size: " << kmerLen << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Reference file path: " << refFileName << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Variants file path: " << vcfFileName << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Ploidy of genotypes in the VCF file: " << vcfPloidy << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Selected GPU ID: " << gpu << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "GPU buffer size: " << buffer << " MB" << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Debug mode: " << (debug ? "Enabled" : "Disabled") << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Fast mode: " << (fastMode ? "Enabled" : "Disabled") << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Use only unique k-mers for indexing: " << (useUniqueKmers ? "Enabled" : "Disabled") << endl;
    }

    // Genotype log output
    void logGenotypeConfigKernel() const {
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

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Selected GPU ID: " << gpu << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "GPU buffer size: " << buffer << " MB" << endl;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Debug mode: " << (debug ? "Enabled" : "Disabled") << endl;
    }
};


class VarigraphKernel : public Varigraph {
public:
    ConstructIndexKernel* ConstructIndexKernelClassPtr_ = nullptr;  // Record the index of graph and reference genome
    int buffer_ = 100;  // Buffer size

    // Constructor
    VarigraphKernel(const VarigraphKernelConfig& config)
        : Varigraph(config) {
        buffer_ = config.buffer;
    }

    // Destructor
    ~VarigraphKernel() {
        if (ConstructIndexKernelClassPtr_ != nullptr) {
            delete ConstructIndexKernelClassPtr_;
            ConstructIndexKernelClassPtr_ = nullptr;
        }
    }

    /**
     * @author zezhen du
     * @date 2024/04/29
     * @version v1.0.1
     * @brief Construct genome graph from reference genome and variants
     * 
     * @return void
    **/
    void construct_kernel();

    /**
     * @brief fastq and genotype
     * 
     * @return void
    */
    void fastq_genotype_kernel();

    /**
     * @author zezhen du
     * @date 2024/04/30
     * @version v1.0.1
     * @brief build the kmer index of files
     * 
     * @param fastqFileNameVec    Sequencing data
     * 
     * @return void
    **/
    void kmer_read_kernel(vector<string> fastqFileNameVec);
};

#endif