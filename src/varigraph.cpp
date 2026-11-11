// g++ -c graph.cpp -std=c++17 -O3 -march=native

#include "../include/varigraph.hpp"
#include "../include/logger.hpp"


/**
 * @author zezhen du
 * @date 2024/01/04
 * @version v1.0.1
 * @brief Construct genome graph from reference genome and variants
 *
 * @return void
**/
void Varigraph::construct() {
    ConstructIndexClassPtr_ = new ConstructIndex(
        refFileName_,
        vcfFileName_,
        inputGraphFileName_,
        outputGraphFileName_,
        fastMode_,
        useUniqueKmers_,
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
    stage_stream() << "Graph construction summary\n";
    log_stream() << "  - Bases: " << ConstructIndexClassPtr_->mGraphBaseNum << '\n';
    log_stream() << "  - K-mers: " << ConstructIndexClassPtr_->mGraphKmerHashHapStrMap.size() << '\n';
    log_stream() << "  - Haplotypes: " << ConstructIndexClassPtr_->mHapMap.size() << '\n';
}


/**
 * @author zezhen du
 * @date 2024/01/04
 * @version v1.0.1
 * @brief Load genome graph
 *
 * @return void
**/
void Varigraph::load() {
    ConstructIndexClassPtr_ = new ConstructIndex(
        refFileName_,
        vcfFileName_,
        inputGraphFileName_,
        outputGraphFileName_,
        fastMode_,
        useUniqueKmers_,
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
    stage_stream() << "Loaded graph summary\n";
    log_stream() << "  - Bases: " << ConstructIndexClassPtr_->mGraphBaseNum << '\n';
    log_stream() << "  - K-mers: " << ConstructIndexClassPtr_->mGraphKmerHashHapStrMap.size() << '\n';
    log_stream() << "  - Haplotypes: " << ConstructIndexClassPtr_->mHapMap.size() << '\n';
}


/**
 * @brief Parse sample configuration file (sampleName, readPath1, readPath2)
 *
 * @return void
*/
void Varigraph::parse_sample_config() {
    stage_stream() << "Reading sample configuration: " << samplesConfigFileName_ << '\n';
    // open file
    GzChunkReader GzChunkReaderClass(samplesConfigFileName_);

    // read line
    string line;
    while (GzChunkReaderClass.read_line(line)) {
        // skip empty lines
        if (line.empty()) continue;

        line = strip(line, '\n');  // remove trailing newline

        // split
        std::istringstream iss(line);
        vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{});

        // check
        if (lineVec.size() <= 1) {
            error_stream() << "Error: The samples configuration file is missing sequencing file information (" << line << ")." << endl;
            exit(1);
        }

        // sample name
        string& sampleName = lineVec[0];
        vector<string> fastqsVec;

        for (uint32_t i = 1; i < lineVec.size(); i++) {
            const auto& filePath = lineVec[i];
            if (filesystem::exists(filePath) && filesystem::file_size(filePath) > 0) {
                fastqsVec.push_back(filePath);
            } else {
                error_stream() << "Error: File '" << filePath << "' does not exist or is empty." << endl;
                exit(1);
            }
        }
        sampleConfigTupleVec_.push_back({sampleName, fastqsVec});
    }

    // Log
    log_stream() << "  - Samples: " << sampleConfigTupleVec_.size() << '\n';
}

/**
 * @brief fastq and genotype
 *
 * @return void
*/
void Varigraph::fastq_genotype() {
    // Merge k-mer information from Genome Graph into nodes.
    ConstructIndexClassPtr_->graph2node();

    // fastq and genotype
    for (size_t sampleIndex = 0; sampleIndex < sampleConfigTupleVec_.size(); ++sampleIndex) {
        const auto& [sampleName, fastqFileNameVec] = sampleConfigTupleVec_[sampleIndex];
        stage_stream() << "Processing sample: " << sampleName << '\n';

        // fastq
        kmer_read(fastqFileNameVec);

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
 * @date 2023/06/27
 * @version v1.0.1
 * @brief build the kmer index of files
 *
 * @param fastqFileNameVec    Sequencing data
 *
 * @return void
**/
void Varigraph::kmer_read(vector<string> fastqFileNameVec) {
    // Computing the frequency of variants and noisy k-mers
    FastqKmer FastqKmerClass(
        ConstructIndexClassPtr_->mGraphKmerHashHapStrMap,
        fastqFileNameVec,
        kmerLen_,
        threads_
    );

    // Calculate coverage using sequencing files
    FastqKmerClass.build_fastq_index();

    // sequencing data depth
    ReadDepth_ = FastqKmerClass.mReadBase / (float)ConstructIndexClassPtr_->mGenomeSize;

    // Calculate Average Coverage (variants k-mers)  2023/12/05
    cal_ave_cov_kmer();

    stage_stream() << "Sample coverage summary\n";
    diagnostic_stream() << fixed << setprecision(2);
    log_stream() << "  - Read bases: " << FastqKmerClass.mReadBase / 1e9 << " Gb\n";
    log_stream() << "  - Read depth: " << ReadDepth_ << '\n';
    log_stream() << "  - Haplotype k-mer coverage: " << hapKmerCoverage_ << '\n';
    diagnostic_stream() << defaultfloat << setprecision(6);
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

    // use sequencing depth as the depth for homozygous k-mers
    if (useDepth_) {
        homCoverage = ReadDepth_ * 0.8;
    }

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
        error_stream()
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
    hapKmerCoverage_ = (homCoverage > 0 && samplePloidy_ > 0) ? static_cast<float>(homCoverage) / static_cast<float>(samplePloidy_) : ReadDepth_ / static_cast<float>(samplePloidy_);
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

    stage_stream() << "K-mer coverage histogram\n";
    log_stream() << "  - Most frequent coverage: " << +maxCoverage << " (" << maxFrequency << " k-mers)\n";
    for (const auto& [coverage, frequence] : kmerCovFreMap) {
        int starNum = static_cast<int>(round(static_cast<double>(frequence) / maxFrequency * 100));  // The number of stars corresponding to k-mer
        if (starNum == 0) {  // If 0, skip
            continue;
        }
        // print hist log
        log_stream() << "  - " << setw(3) << +coverage << ": ";
        for (int j = 0; j < min(starNum, 100); ++j) {
            diagnostic_stream() << '*';
        }
        if (starNum > 100) {
            diagnostic_stream() << ">";
        }
        diagnostic_stream() << " " << frequence << endl;
    }
    log_stream() << "  - Homozygous peak: " << +homCoverage << "; haplotype coverage: " << hapKmerCoverage_ << '\n';
}


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
void Varigraph::genotype(string sampleName) {
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
        sampleName,
        kmerLen_,
        haploidNum_,
        chrLenThread_,
        transitionProType_,
        svGenotypeBool_,
        threads_,
        debug_,
        minSupportingGQ_
    );
}
