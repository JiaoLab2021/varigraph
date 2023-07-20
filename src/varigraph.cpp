// g++ -c graph.cpp -std=c++17 -O3 -march=native

#include "../include/varigraph.hpp"


Varigraph::Varigraph(
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
) : refFileName_(refFileName), fastqFileNameVec_(fastqFileNameVec), vcfFileName_(vcfFileName), inputFastqKmerFileName_(inputFastqKmerFileName), inputMbfFileName_(inputMbfFileName), 
    outputFastqKmerFileName_(outputFastqKmerFileName), outputMbfFileName_(outputMbfFileName),  outputFileName_(outputFileName), 
    kmerLen_(kmerLen), prefix_(prefix), ploidy_(ploidy), threads_(threads), debug_(debug) {}


Varigraph::~Varigraph()
{
    if (ConstructIndexClassPtr_ != nullptr) {
        delete ConstructIndexClassPtr_;
        ConstructIndexClassPtr_ = nullptr;
    }
};

/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief build the kmer index of reference and construct graph
 * 
 * @return void
**/
void Varigraph::ref_idx_construct()
{
    ConstructIndexClassPtr_ = new ConstructIndex(
        refFileName_, 
        vcfFileName_, 
        inputMbfFileName_, 
        outputMbfFileName_, 
        kmerLen_, 
        prefix_, 
        ploidy_, 
        threads_, 
        debug_
    );

    // Building the k-mers index of reference genome
    ConstructIndexClassPtr_->build_fasta_index();

    // Counting Bloom Filter
    ConstructIndexClassPtr_->make_mbf();

    // Graph construction
    ConstructIndexClassPtr_->construct();

    // building the k-mer index of graph
    ConstructIndexClassPtr_->index();
    cerr << endl;
    cerr << "           - " << "Number of k-mers in the graph: " << ConstructIndexClassPtr_->mGraphKmerCovFreMap.size() << endl;
    cerr << "           - " << "Number of haplotypes in the graph: " << ConstructIndexClassPtr_->mHapMap.size() << endl << endl << endl;

    // Free memory (Counting Bloom Filter)
    ConstructIndexClassPtr_->clear_memory();

    // K-mer deduplication in the graph genome
    ConstructIndexClassPtr_->kmer_deduplication();
    cerr << endl;
    cerr << "           - " << "Number of unique k-mers in the graph: " << ConstructIndexClassPtr_->mGraphKmerCovFreMap.size() << endl << endl << endl;
}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief build the kmer index of files
 * 
 * @return void
**/
void Varigraph::kmer_read()
{
    // Computing the frequency of variants and noisy k-mers
    FastqKmer FastqKmerClass(
        ConstructIndexClassPtr_->mGraphKmerCovFreMap, 
        fastqFileNameVec_, 
        kmerLen_, 
        threads_
    );
    
    // Calculate coverage using sequencing files
    if (inputFastqKmerFileName_.empty()) {
        FastqKmerClass.build_fastq_index();

        // save to file
        if (!outputFastqKmerFileName_.empty()) {
            FastqKmerClass.save_index(outputFastqKmerFileName_);
        }
    } else {  // load coverage from file
        FastqKmerClass.load_index(inputFastqKmerFileName_);
    }

    // Calculate Average Coverage (variants k-mers)
    uint64_t allKmerCoverage = std::accumulate(
        ConstructIndexClassPtr_->mGraphKmerCovFreMap.begin(), 
        ConstructIndexClassPtr_->mGraphKmerCovFreMap.end(), 
        0ull,
        [](int sum, const auto& pair)
        {
            return sum + pair.second.c;
        }
    );
    float aveKmerCoverage = (ConstructIndexClassPtr_->mGraphKmerCovFreMap.size() > 0) ? static_cast<float>(allKmerCoverage) / ConstructIndexClassPtr_->mGraphKmerCovFreMap.size() : 0.0f;
    
    cerr << endl;
    cerr << "           - " << "Average coverage of Variants k-mers: " << aveKmerCoverage << endl << endl << endl;  // print log
}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief genotype
 * 
 * 
 * @return void
**/
void Varigraph::genotype()
{
    // genotype
    GENOTYPE::genotype(
        ConstructIndexClassPtr_->mGraphKmerCovFreMap, 
        ConstructIndexClassPtr_->mGraphMap, 
        ConstructIndexClassPtr_->mHapMap, 
        ConstructIndexClassPtr_->mVcfHead, 
        ConstructIndexClassPtr_->mVcfInfoMap, 
        outputFileName_, 
        threads_, 
        debug_
    );
}