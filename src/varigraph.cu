// g++ -c graph.cpp -std=c++17 -O3 -march=native

#include "../include/varigraph.cuh"
#include "../include/logger.hpp"


/**
 * @author zezhen du
 * @date 2024/04/29
 * @version v1.0.1
 * @brief Construct genome graph from reference genome and variants
 *
 * @return void
**/
void VarigraphKernel::construct_kernel() {
    ConstructIndexKernelClassPtr_ = new ConstructIndexKernel(
        refFileName_,
        vcfFileName_,
        inputGraphFileName_,
        outputGraphFileName_,
        fastMode_,
        useUniqueKmers_,
        kmerLen_,
        vcfPloidy_,
        debug_,
        threads_,
        buffer_
    );

    // Building the k-mers index of reference genome
    ConstructIndexKernelClassPtr_->build_fasta_index();

    // Counting Bloom Filter
    ConstructIndexKernelClassPtr_->make_mbf_kernel();

    // Genome Graph construction
    ConstructIndexKernelClassPtr_->construct();

    // make mHapIdxQRmap
    ConstructIndexKernelClassPtr_->make_QRmap();

    // building the k-mer index of graph
    ConstructIndexKernelClassPtr_->index_kernel();

    // save Genome Graph to file
    ConstructIndexKernelClassPtr_->save_index();

    // Free memory (Counting Bloom Filter)
    ConstructIndexKernelClassPtr_->clear_mbf();

    // log
    stage_stream() << "Graph construction summary\n";
    log_stream() << "  - Bases: " << ConstructIndexKernelClassPtr_->mGraphBaseNum << '\n';
    log_stream() << "  - K-mers: " << ConstructIndexKernelClassPtr_->mGraphKmerHashHapStrMap.size() << '\n';
    log_stream() << "  - Haplotypes: " << ConstructIndexKernelClassPtr_->mHapMap.size() << '\n';
}

/**
 * @brief fastq and genotype
 *
 * @return void
*/
void VarigraphKernel::fastq_genotype_kernel() {
    // Merge k-mer information from Genome Graph into nodes.
    ConstructIndexClassPtr_->graph2node();

    // fastq and genotype
    for (size_t sampleIndex = 0; sampleIndex < sampleConfigTupleVec_.size(); ++sampleIndex) {
        const auto& [sampleName, fastqFileNameVec] = sampleConfigTupleVec_[sampleIndex];
        stage_stream() << "Processing sample: " << sampleName << '\n';

        // fastq
        kmer_read_kernel(fastqFileNameVec);

        // genotype
        genotype(sampleName);

        // Only another sample needs the graph coverage and posterior state reset.
        if (sampleIndex + 1 < sampleConfigTupleVec_.size()) {
            ConstructIndexClassPtr_->reset();
        }
        log_stream() << "  - Sample completed: " << sampleName << '\n';
    }
}

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
void VarigraphKernel::kmer_read_kernel(vector<string> fastqFileNameVec) {
    // Computing the frequency of variants and noisy k-mers
    FastqKmerKernel FastqKmerKernelClass(
        ConstructIndexClassPtr_->mGraphKmerHashHapStrMap,
        fastqFileNameVec,
        kmerLen_,
        threads_,
        buffer_
    );

    // Calculate coverage using sequencing files
    FastqKmerKernelClass.build_fastq_index_kernel();  // GPU

    // sequencing data depth
    ReadDepth_ = FastqKmerKernelClass.mReadBase / (float)ConstructIndexClassPtr_->mGenomeSize;

    // Calculate Average Coverage (variants k-mers)  2023/12/05
    cal_ave_cov_kmer();

    stage_stream() << "Sample coverage summary\n";
    diagnostic_stream() << fixed << setprecision(2);
    log_stream() << "  - Read bases: " << FastqKmerKernelClass.mReadBase / 1e9 << " Gb\n";
    log_stream() << "  - Read depth: " << ReadDepth_ << '\n';
    log_stream() << "  - Haplotype k-mer coverage: " << hapKmerCoverage_ << '\n';
    diagnostic_stream() << defaultfloat << setprecision(6);
}
