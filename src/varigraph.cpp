// g++ -c graph.cpp -std=c++17 -O3 -march=native

#include "../include/varigraph.hpp"


Varigraph::Varigraph(
    const string& refFileName, 
    const string& vcfFileName, 
    const vector<string>& fastqFileNameVec, 
    const string& inputGraphFileName, 
    const string& outputGraphFileName, 
    const string& outputFileName, 
    const bool& fastMode, 
    uint32_t& kmerLen, 
    const string& sampleName, 
    const string& sampleType, 
    const uint32_t& samplePloidy, 
    uint32_t& vcfPloidy, 
    const uint32_t& haploidNum, 
    const uint32_t& chrLenThread, 
    const string& transitionProType, 
    const bool& svGenotypeBool, 
    const bool& debug, 
    const uint32_t& threads
) : refFileName_(refFileName), vcfFileName_(vcfFileName), fastqFileNameVec_(fastqFileNameVec), inputGraphFileName_(inputGraphFileName), 
    outputGraphFileName_(outputGraphFileName), outputFileName_(outputFileName), 
    fastMode_(fastMode), kmerLen_(kmerLen), sampleName_(sampleName), sampleType_(sampleType), samplePloidy_(samplePloidy), vcfPloidy_(vcfPloidy), haploidNum_(haploidNum), 
    chrLenThread_(chrLenThread), transitionProType_(transitionProType), svGenotypeBool_(svGenotypeBool), debug_(debug), threads_(threads) {
    // cerr
    std::cerr.imbue(std::locale(""));  // Thousandth output
}


Varigraph::~Varigraph() {
    if (ConstructIndexClassPtr_ != nullptr) {
        delete ConstructIndexClassPtr_;
        ConstructIndexClassPtr_ = nullptr;
    }
};


/**
 * @author zezhen du
 * @date 2024/01/04
 * @version v1.0.1
 * @brief Construct genome graph from reference genome and variants
 * 
 * @return void
**/
void Varigraph::construct()
{
    ConstructIndexClassPtr_ = new ConstructIndex(
        refFileName_, 
        vcfFileName_, 
        inputGraphFileName_, 
        outputGraphFileName_, 
        fastMode_, 
        kmerLen_, 
        vcfPloidy_, 
        debug_, 
        threads_
    );

    // Building the k-mers index of reference genome
    ConstructIndexClassPtr_->build_fasta_index();

    // Counting Bloom Filter
    ConstructIndexClassPtr_->make_mbf();

    // Genome Graph construction
    ConstructIndexClassPtr_->construct();

    // make mHapIdxQRmap
    ConstructIndexClassPtr_->make_QRmap();

    // building the k-mer index of graph
    ConstructIndexClassPtr_->index();

    // save Genome Graph to file
    ConstructIndexClassPtr_->save_index();

    // Free memory (Counting Bloom Filter)
    ConstructIndexClassPtr_->clear_mbf();

    // log
    cerr << endl;
    cerr << "           - " << "Total number of bases in the Genome Graph: " << ConstructIndexClassPtr_->mGraphBaseNum << endl;
    cerr << "           - " << "Total number of k-mers present in the Genome Graph: " << ConstructIndexClassPtr_->mGraphKmerHashHapStrMap.size() << endl;
    cerr << "           - " << "Total number of haplotypes present in the Genome Graph: " << ConstructIndexClassPtr_->mHapMap.size() << endl << endl << endl;
}


/**
 * @author zezhen du
 * @date 2024/01/04
 * @version v1.0.1
 * @brief Load genome graph
 * 
 * @return void
**/
void Varigraph::load()
{
    ConstructIndexClassPtr_ = new ConstructIndex(
        refFileName_, 
        vcfFileName_, 
        inputGraphFileName_, 
        outputGraphFileName_, 
        false, 
        kmerLen_, 
        vcfPloidy_, 
        debug_, 
        threads_
    );

    // load Genome Graph from file
    ConstructIndexClassPtr_->load_index();

    // make mHapIdxQRmap
    ConstructIndexClassPtr_->make_QRmap();

    // k-mer length
    kmerLen_ = ConstructIndexClassPtr_->mKmerLen;

    // vcf ploidy
    vcfPloidy_ = ConstructIndexClassPtr_->mVcfPloidy;
    
    // log
    cerr << endl;
    cerr << "           - " << "Total number of bases in the Genome Graph: " << ConstructIndexClassPtr_->mGraphBaseNum << endl;
    cerr << "           - " << "Total number of k-mers present in the Genome Graph: " << ConstructIndexClassPtr_->mGraphKmerHashHapStrMap.size() << endl;
    cerr << "           - " << "Total number of haplotypes present in the Genome Graph: " << ConstructIndexClassPtr_->mHapMap.size() << endl << endl << endl;
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
    FastqKmerClass.build_fastq_index();

    // sequencing data depth
    ReadDepth_ = FastqKmerClass.mReadBase / (float)ConstructIndexClassPtr_->mGenomeSize;

    // Calculate Average Coverage (variants k-mers)  2023/12/05
    cal_ave_cov_kmer();

    cerr << endl;
    cerr << fixed << setprecision(2);
    cerr << "           - " << "Size of the sequenced data: " << FastqKmerClass.mReadBase / 1e9 << " Gb" << endl;
    cerr << "           - " << "Depth of the sequenced data: " << ReadDepth_ << endl;
    cerr << "           - " << "Coverage of haplotype k-mers: " << hapKmerCoverage_ << endl << endl << endl;
    cerr << defaultfloat << setprecision(6);

    // Merge k-mer information from Genome Graph into nodes.
    ConstructIndexClassPtr_->graph2node();
}


/**
 * @author zezhen du
 * @date 2024/01/19
 * @version v1.0.1
 * @brief calculate the average coverage of k-mers
 * 
 * @return void
**/
void Varigraph::cal_ave_cov_kmer() {
    // homozygous k-mer (variants k-mers, map<coverage, frequence>)
    map<uint8_t, uint64_t> kmerCovFreMap = get_hom_kmer();

    // Get k-mers frequency
    uint8_t maxCoverage;  // Maximum coverage of k-mer
    uint8_t homCoverage;  // Homozygous k-mer coverage
    tie(maxCoverage, homCoverage) = get_hom_kmer_c(kmerCovFreMap);

    // haplotype k-mers coverage
    cal_hap_kmer_cov(homCoverage);

    // print hist lines
    kmer_histogram(
        maxCoverage, 
        homCoverage, 
        kmerCovFreMap
    );
}

/**
 * @author zezhen du
 * @date 2024/01/22
 * @version v1.0.1
 * @brief Get homozygous k-mer
 * 
 * @return kmerCovFreMap   map<coverage, frequency>
**/
map<uint8_t, uint64_t> Varigraph::get_hom_kmer() {
    // Calculate Average Coverage (variants k-mers)
    map<uint8_t, uint64_t> kmerCovFreMap;  // map<coverage, frequency>

    // Quotient and remainder for each haplotype in topHapVec
    const unordered_map<uint16_t, tuple<uint16_t, uint16_t> >& hapIdxQRmap = ConstructIndexClassPtr_->mHapIdxQRmap;  // map<hapIdx, tuple<quotient, remainder> >

    for (const auto& pair : ConstructIndexClassPtr_->mGraphKmerHashHapStrMap) {
        uint8_t c = pair.second.c;

        if (c == 0 || pair.second.f > 1) {  // If the coverage is 0 or the frequency is greater than 1, the k-mer is skipped
            continue;
        }

        int totalHomSampleCount = 0;  // The number of homozygous samples for this k-mer

        int index = 0;  // haplotype index
        int sampleCount = 0;  // the number of sample contained
        for (size_t i = 1; i < ConstructIndexClassPtr_->mHapNum; i++) {
            index++;

            // Determine whether the haplotype contains the k-mer
            if (construct_index::get_bit(pair.second.BitVec[get<0>(hapIdxQRmap.at(i))], get<1>(hapIdxQRmap.at(i))) > 0) {
                sampleCount++;
            }

            if (index == vcfPloidy_) {  // one sample
                index = 0;
                if (sampleCount == vcfPloidy_) {  // Only record homozygous k-mers
                    totalHomSampleCount++;
                    break;  // If it is a homozygous k-mer, because totalHomSampleCount is already equal to 1, jump out of the for loop directly.
                }
                sampleCount = 0;
            }
        }

        // Is the number of homozygous samples greater than 0?
        if (totalHomSampleCount > 0) {
            auto& emplacedValue = kmerCovFreMap.emplace(c, 0).first->second;
            emplacedValue++;
        }
    }
    return move(kmerCovFreMap);
}

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
tuple<uint8_t, uint8_t> Varigraph::get_hom_kmer_c(const map<uint8_t, uint64_t>& kmerCovFreMap) {
    // Get the k-mers with the highest frequency
    vector<uint8_t> coverageVec;
    vector<uint64_t> frequencyVec;
    int16_t index = -1;
    int16_t maxIndex = -1;
    uint8_t maxCoverage = 0;
    uint64_t maxFrequency = 0;
    uint8_t homCoverage = 0;  // Homozygous k-mer coverage
    for (const auto& [coverage, frequency] : kmerCovFreMap) {
        coverageVec.push_back(coverage);
        frequencyVec.push_back(frequency);
        index++;
        if (coverage > 1 && frequency >= maxFrequency && coverage < UINT8_MAX) {
            maxIndex = index;
            maxCoverage = coverage;
            maxFrequency = frequency;
            homCoverage = coverage;
        }
    }

    // Check if the data is correct
    if (maxIndex == -1) {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "Error: Failed to retrieve depth information of k-mers from the sequencing data. Please verify your data." << endl;
        exit(1);
    }
    
    // look for smaller peak on the right
    for (size_t i = maxIndex + 1; i < frequencyVec.size() - 1; i++) {
        if (coverageVec[i] > ReadDepth_) {  // the peak value must be smaller than the sequencing depth
            break;
        }
        
        // Determine whether it is a peak value
        if (frequencyVec[i] >= frequencyVec[i-1] && frequencyVec[i] >= frequencyVec[i+1]) {  // peak
            homCoverage = coverageVec[i];
        }
    }
    return {maxCoverage, homCoverage};
}

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
void Varigraph::cal_hap_kmer_cov(const uint8_t& homCoverage) {
    hapKmerCoverage_ = (homCoverage > 0 && samplePloidy_ > 0) ? static_cast<float>(homCoverage) / static_cast<float>(vcfPloidy_) : ReadDepth_ / static_cast<float>(vcfPloidy_);
}

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
void Varigraph::kmer_histogram(
    const uint8_t& maxCoverage, 
    const uint8_t& homCoverage, 
    const map<uint8_t, uint64_t>& kmerCovFreMap
) {
    uint64_t maxFrequency = kmerCovFreMap.at(maxCoverage);

    cerr << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "highest: count[" << +maxCoverage << "] = " << maxFrequency << endl;
    for (const auto& [coverage, frequence] : kmerCovFreMap) {
        int starNum = static_cast<int>(round(static_cast<double>(frequence) / maxFrequency * 100));  // The number of stars corresponding to k-mer
        if (starNum == 0) {  // If 0, skip
            continue;
        }
        // print hist log
        cerr << "[" << __func__ << "::" << getTime() << "] " << setw(3) << +coverage << ": ";
        for (int j = 0; j < min(starNum, 100); ++j) {
            cerr << '*';
        }
        if (starNum > 100) {
            cerr << ">";
        }
        cerr << " " << frequence << endl;
    }
    cerr << "[" << __func__ << "::" << getTime() << "] " << "peak_hom: " << +homCoverage << "; peak_hap: " << hapKmerCoverage_ << endl;
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
        ConstructIndexClassPtr_->mFastaLenMap, 
        ConstructIndexClassPtr_->mGraphMap, 
        ConstructIndexClassPtr_->mHapMap, 
        ConstructIndexClassPtr_->mVcfHead, 
        ConstructIndexClassPtr_->mVcfInfoMap, 
        ConstructIndexClassPtr_->mHapIdxQRmap, 
        sampleType_, 
        samplePloidy_, 
        hapKmerCoverage_, 
        sampleName_, 
        outputFileName_, 
        kmerLen_, 
        haploidNum_, 
        chrLenThread_, 
        transitionProType_, 
        svGenotypeBool_, 
        threads_, 
        debug_
    );
}