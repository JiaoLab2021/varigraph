// g++ -c graph.cpp -std=c++17 -O3 -march=native

#include "../include/varigraph.hpp"


Varigraph::Varigraph(
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
) : refFileName_(refFileName), fastqFileNameVec_(fastqFileNameVec), vcfFileName_(vcfFileName), inputMbfFileName_(inputMbfFileName), inputGraphFileName_(inputGraphFileName), inputFastqKmerFileName_(inputFastqKmerFileName), 
    outputMbfFileName_(outputMbfFileName), outputGraphFileName_(outputGraphFileName), outputFastqKmerFileName_(outputFastqKmerFileName), outputFileName_(outputFileName), 
    fastMode_(fastMode), kmerLen_(kmerLen), sampleName_(sampleName), genomeType_(genomeType), refPloidy_(refPloidy), vcfPloidy_(vcfPloidy), haploidNum_(haploidNum), debug_(debug), threads_(threads) {
        // cerr
        std::cerr.imbue(std::locale(""));  // Thousandth output
    }


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
        inputGraphFileName_, 
        outputMbfFileName_, 
        outputGraphFileName_, 
        fastMode_, 
        kmerLen_, 
        sampleName_, 
        vcfPloidy_, 
        debug_, 
        threads_
    );

    // Building the k-mers index of reference genome
    ConstructIndexClassPtr_->build_fasta_index();

    if (inputGraphFileName_.empty()) {
        // Counting Bloom Filter
        ConstructIndexClassPtr_->make_mbf();

        // Genome Graph construction
        ConstructIndexClassPtr_->construct();

        // building the k-mer index of graph
        ConstructIndexClassPtr_->index();

        // save Genome Graph to file
        if (!outputGraphFileName_.empty()) {
            ConstructIndexClassPtr_->save_index();
        }

        // Free memory (Counting Bloom Filter)
        ConstructIndexClassPtr_->clear_mbf();
        
    } else {
        // VCF index
        ConstructIndexClassPtr_->construct();

        // load Genome Graph from file
        ConstructIndexClassPtr_->load_index();
    }
    
    // log
    cerr << endl;
    cerr << "           - " << "Number of k-mers in the Genome Graph: " << ConstructIndexClassPtr_->mGraphKmerHashHapStrMap.size() << endl;
    cerr << "           - " << "Number of haplotypes in the Genome Graph: " << ConstructIndexClassPtr_->mHapMap.size() << endl << endl << endl;
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
        ConstructIndexClassPtr_->mGraphKmerHashHapStrMap, 
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

    // sequencing data depth
    ReadDepth_ = FastqKmerClass.mReadBase / (float)ConstructIndexClassPtr_->mGenomeSize;

    // Calculate Average Coverage (variants k-mers)
    uint64_t allKmerNum = 0;
    uint64_t allKmerCov = std::accumulate(
        ConstructIndexClassPtr_->mGraphKmerHashHapStrMap.begin(), 
        ConstructIndexClassPtr_->mGraphKmerHashHapStrMap.end(), 
        0ull,
        [&allKmerNum](int sum, const auto& pair) mutable {
            // if (pair.second.c > 0 && pair.second.f == 1) {
            if (pair.second.c > 1 && pair.second.f == 1) {  // 2023/09/27
                ++allKmerNum;
                return sum + pair.second.c;
            } else {
                return sum;
            }
        }
    );
    float aveKmerCoverage = (allKmerNum > 0) ? static_cast<float>(allKmerCov) / allKmerNum : 0.0f;

    cerr << endl;
    cerr << fixed << setprecision(2);
    cerr << "           - " << "Size of sequencing data: " << FastqKmerClass.mReadBase / 1e9 << " Gb" << endl;
    cerr << "           - " << "Sequencing data depth: " << ReadDepth_ << endl;
    cerr << "           - " << "Average coverage of k-mers: " << aveKmerCoverage << endl << endl << endl;
    cerr << defaultfloat << setprecision(6);

    // Merge k-mer information from Genome Graph into nodes.
    ConstructIndexClassPtr_->graph2node();
    // ConstructIndexClassPtr_->clear_mGraphKmerCovFreMap();
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
        ConstructIndexClassPtr_->mGraphMap, 
        ConstructIndexClassPtr_->mHapMap, 
        ConstructIndexClassPtr_->mVcfHead, 
        ConstructIndexClassPtr_->mVcfInfoMap, 
        genomeType_, 
        refPloidy_, 
        ReadDepth_, 
        outputFileName_, 
        kmerLen_, 
        haploidNum_, 
        threads_, 
        debug_
    );
}