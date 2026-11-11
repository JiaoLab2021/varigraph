#ifndef VARIGRAPH_CUH
#define VARIGRAPH_CUH

#include "construct_index.cuh"
#include "fastq_kmer.cuh"

#include "varigraph.hpp"

using namespace std;


class VarigraphKernel : public Varigraph {
public:
    ConstructIndexKernel* ConstructIndexKernelClassPtr_ = nullptr;  // Record the index of graph and reference genome
    int buffer_ = 100;  // Buffer size

    // Constructor
    VarigraphKernel(const VarigraphConfig& config, const GpuOptions& gpu)
        : Varigraph(config), buffer_(gpu.buffer) {}

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
