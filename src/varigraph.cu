// g++ -c graph.cpp -std=c++17 -O3 -march=native

#include "../include/varigraph.cuh"


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
    cerr << endl;
    cerr << "           - " << "Total number of bases in the Genome Graph: " << ConstructIndexKernelClassPtr_->mGraphBaseNum << endl;
    cerr << "           - " << "Total number of k-mers present in the Genome Graph: " << ConstructIndexKernelClassPtr_->mGraphKmerHashHapStrMap.size() << endl;
    cerr << "           - " << "Total number of haplotypes present in the Genome Graph: " << ConstructIndexKernelClassPtr_->mHapMap.size() << endl << endl << endl;
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
    for (const auto& [sampleName, fastqFileNameVec] : sampleConfigTupleVec_) {  // vector<tuple<sampleName, vector<readPath> > >
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Processing sample: " << sampleName << endl << endl;

        // fastq
        kmer_read_kernel(fastqFileNameVec);

        // genotype
        genotype(sampleName);

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Sample: " << sampleName << " has been processed." << endl << endl << endl;

        // Reset ConstructIndexClassPtr_
        ConstructIndexClassPtr_->reset();
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

    cerr << endl;
    cerr << fixed << setprecision(2);
    cerr << "           - " << "Size of the sequenced data: " << FastqKmerKernelClass.mReadBase / 1e9 << " Gb" << endl;
    cerr << "           - " << "Depth of the sequenced data: " << ReadDepth_ << endl;
    cerr << "           - " << "Coverage of haplotype k-mers: " << hapKmerCoverage_ << endl << endl << endl;
    cerr << defaultfloat << setprecision(6);
}