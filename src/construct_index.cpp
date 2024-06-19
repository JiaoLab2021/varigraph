// g++ -c construct_index.cpp -std=c++17 -lz -O3 -march=native

#include "../include/kmer.hpp"
#include "../include/construct_index.hpp"

using namespace std;

// global variable
bool debugConstruct = false;

std::mutex mtxCI;

// kseq.h
KSEQ_INIT(gzFile, gzread)


ConstructIndex::ConstructIndex(
    const string& refFileName, 
    const string& vcfFileName, 
    const string& inputGraphFileName, 
    const string& outputGraphFileName, 
    const bool& fastMode, 
    const bool& useUniqueKmers, 
    const uint32_t& kmerLen, 
    const uint32_t& vcfPloidy, 
    const bool& debug, 
    const uint32_t& threads
) : refFileName_(refFileName), vcfFileName_(vcfFileName), inputGraphFileName_(inputGraphFileName), outputGraphFileName_(outputGraphFileName), 
    fastMode_(fastMode), useUniqueKmers_(useUniqueKmers), mKmerLen(kmerLen), mVcfPloidy(vcfPloidy), threads_(threads) {

    mHapMap[0] = "reference";

    // debug
    debugConstruct = debug;
    if (debugConstruct) threads_ = 1;
}

ConstructIndex::~ConstructIndex() {
    if (mbf != nullptr) {
        delete mbf;  // Release the memory occupied by the Counting Bloom filter
        mbf = nullptr;
    }
}

/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief Free memory
 * 
 * @return void
**/
void ConstructIndex::clear_mbf() {
    unordered_map<string, string>().swap(mFastaSeqMap);
    if (mbf != nullptr) {
        delete mbf;
        mbf = nullptr;
    }
    malloc_trim(0);	// 0 is for heap memory
}

/**
 * @author zezhen du
 * @date 2023/07/21
 * @version v1.0.1
 * @brief Free memory
 * 
 * @return void
**/
void ConstructIndex::clear_mGraphKmerCovFreMap() {
    // Clear unordered_map
    unordered_map<uint64_t, kmerCovFreBitVec>().swap(mGraphKmerHashHapStrMap);

    malloc_trim(0);	// 0 is for heap memory
}

/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief Building the k-mers index of reference genome
 * 
 * @return void
**/
void ConstructIndex::build_fasta_index()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building refgenome index: " << refFileName_ << endl;

    // open fasta file
    gzFile gzfp = gzopen(refFileName_.c_str(), "rb");

    if(!gzfp) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << refFileName_ << "': No such file or directory." << endl;
        exit(1);
    }
    else {
        kseq_t *ks;
        ks = kseq_init(gzfp);

        // Define local variables to avoid frequent application and release of memory inside the loop
        string chromosome;
        string sequence;
        uint32_t sequenceLen;

        while( kseq_read(ks) >= 0 ) {
            // ks->name.s  name
            // ks->seq.s   sequence
            chromosome = ks->name.s;
            sequence = ks->seq.s;
            sequenceLen = ks->seq.l;

            // record the length of sequence
            mGenomeSize += ks->seq.l;

            // build fasta index
            mFastaSeqMap.emplace(chromosome, sequence);
            mFastaLenMap.emplace(chromosome, sequenceLen);

            // Check if chromosome length is greater than UINT32_MAX
            if (sequence.length() > UINT32_MAX) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'" << chromosome << "' length is greater than 4,294,967,295." << endl;
                exit(1);
            }
        }

        // free memory and close file
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    cerr << endl;
    cerr << fixed << setprecision(2);
    cerr << "           - " << "Size of reference genome: " << mGenomeSize / 1000.0 / 1000.0 << " Mb" << endl << endl << endl;
    cerr << defaultfloat << setprecision(6);

    malloc_trim(0);	// 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2023/07/18
 * @version v1.0.1
 * @brief Making Counting Bloom Filter
 * 
 * @return void
**/
void ConstructIndex::make_mbf() {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Initiating computation of k-mer frequencies in the reference genome on CPU ...\n";

    /* *************************************************** making or load *************************************************** */
    uint64_t bfSize = mGenomeSize - mKmerLen + 1;
    double errorRate = 0.01;
    mbf = new BloomFilter(bfSize, errorRate);

    // making
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Making Counting Bloom Filter with a false positive rate of " << errorRate << " ...\n";

    for (const auto& [chromosome, sequence] : mFastaSeqMap) {  // map<chromosome, sequence>
        // Constructing k-mer index
        kmerBit::kmer_sketch_bf(sequence, mKmerLen, mbf);

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Chromosome '" << chromosome << "' processed successfully ...\n";
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Counting Bloom Filter constructed successfully ..." << endl << endl;
    
    cerr << "           - " << "Counting Bloom Filter size: " << mbf->get_size() << endl;
    cerr << "           - " << "Hash functions count: " << mbf->get_num() << endl;
    cerr << fixed << setprecision(2);
    cerr << "           - " << "Counting Bloom Filter usage rate: " << mbf->get_cap() << endl << endl << endl;
    cerr << defaultfloat << setprecision(6);

    malloc_trim(0);	// 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0
 * @brief Graph construction
 * 
 * @return void
**/
void ConstructIndex::construct()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Constructing ...\n";  // print log

    mGraphBaseNum += mGenomeSize;  // The number of bases in the graph

    // Record the position of the previous node
    uint32_t tmpRefStart = 0;
    uint32_t tmpRefEnd = 0;
    string tmpChromosome;

    // open file
    GzChunkReader GzChunkReaderClass(vcfFileName_);

    // read line
    string line;
    while (GzChunkReaderClass.read_line(line)) {
        // skip empty lines
        if (line.empty()) {
            continue;
        }

        line = strip(line, '\n');  // remove trailing newline

        // First judge whether it is a comment line, if yes, skip
        if (line.find("##FORMAT") != string::npos) {
            continue;
        }
        if (line.find("#") != string::npos && line.find("#CHROM") == string::npos) {
            mVcfHead += line + "\n";  // Store VCF file comment line
            continue;
        }

        // Split the line into tokens based on whitespace
        std::istringstream iss(line);
        vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

        // Check if the file is correct, if not, jump out of the code
        if (lineVec.size() < 10) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error in '" << vcfFileName_  << "': Number of columns in the VCF file is less than 10. Current column count: " << lineVec.size() << endl;
            exit(1);
        }

        // store haplotype information
        if (line.find("#CHROM") != string::npos) {
            mVcfHead += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
            mVcfHead += "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality (phred-scaled 1 - max(GPP))\">\n";
            mVcfHead += "##FORMAT=<ID=GPP,Number=1,Type=String,Description=\"Genotype posterior probabilities\">\n";
            mVcfHead += "##FORMAT=<ID=NAK,Number=.,Type=Float,Description=\"Number of allele k-mers\">\n";
            mVcfHead += "##FORMAT=<ID=CAK,Number=.,Type=Float,Description=\"Coverage of allele k-mers\">\n";
            mVcfHead += "##FORMAT=<ID=UK,Number=1,Type=Integer,Description=\"Total number of unique kmers, capped at 255\">\n";

            mVcfHead += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";  // Store VCF file comment line

            uint16_t hapIdx = 1;  // haplotype index
            for (size_t i = 9; i < lineVec.size(); i++) {
                for (size_t j = 0; j < mVcfPloidy; j++) {
                    mHapMap[hapIdx] = lineVec[i];  // assignment
                
                    if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                        hapIdx++;  // Haplotype index +1
                    } else {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: The number of haplotypes exceeds the maximum limit of " << UINT16_MAX << "." << endl;
                        exit(1);
                    }
                }
            }

            // haplotype number
            mHapNum = mHapMap.size();
        } else {
            string chromosome = lineVec[0];  // 0
            uint32_t refStart = stoul(lineVec[1]);  // 1
            string refSeq = lineVec[3];  // 3
            uint32_t refLen = refSeq.size();  // the length of reference sequence
            uint32_t refEnd = refStart + refSeq.size() - 1;  // end position
            string qrySeq = lineVec[4];  // 4
            vector<string> qrySeqVec = split(qrySeq, ",");  // the vector of ALT

            // FORMAT
            vector<string> formatVec = split(strip(lineVec[8], '\n'), ":");
            vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

            int gtIndex = distance(formatVec.begin(), gtItera);
            if (gtIndex == formatVec.size()) {  // FORMAT: GT
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: Genotype (GT) information is missing in FORMAT: " << line << endl;
                exit(1);
            }

            // VCF index
            vcf_construct(
                chromosome, 
                refStart, 
                refLen, 
                lineVec, 
                qrySeqVec, 
                gtIndex
            );

            // Check if the chromosome is in the reference genome
            auto mFastaMapFindIter = mFastaSeqMap.find(chromosome);
            if (mFastaMapFindIter == mFastaSeqMap.end()) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: Chromosome '" << chromosome << "' not found in reference genome." << endl;
                exit(1);
            }

            // If it overlaps with the previous coordinate, skip it
            if (chromosome != tmpChromosome) {
                tmpRefStart = 0;
            }
            if (tmpRefStart == refStart) {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: Multiple variants detected, skipping this site -> " << chromosome << " " << refStart << endl;
                continue;
            } else if (tmpRefStart > refStart) {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: Variants are unsorted, skipping this site -> " << chromosome << " " << tmpRefStart << ">" << refStart << endl;
                continue;
            }
            
            // Check that the REF sequence is consistent
            string trueRefSeq = mFastaMapFindIter->second.substr(refStart-1, refSeq.length());
            if (trueRefSeq != refSeq) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Warning: Sequence discrepancy detected between reference genome and VCF. Replacing with sequence from reference genome -> "
                    << chromosome  << "\t" << refStart << endl;
                refSeq = trueRefSeq;
            }

            // genotype=0
            // If the first vcf position of the chromosome is not 1, the first node is constructed from the sequence at the front of the genome
            uint16_t genotype = 0;
            if (chromosome != tmpChromosome) {
                // Add the first part of the new chromosome
                string preRefSeq;
                uint32_t preRefStart;
                uint32_t preRefEnd;
                uint32_t preRefLen;
                
                // Add the second half of the old chromosome
                if (tmpRefEnd > 0 && tmpRefEnd < mFastaSeqMap[tmpChromosome].length()) {  // See if the last variant is at the end of the chromosome
                    preRefStart = tmpRefEnd + 1;  // Previous starting position
                    preRefEnd = mFastaSeqMap[tmpChromosome].length();  // The last END position, the end of the chromosome
                    preRefLen = preRefEnd - preRefStart + 1;  // The length of the previous sequence
                    preRefSeq = mFastaSeqMap[tmpChromosome].substr(preRefStart-1, preRefLen);  // The sequence information of the last part of the previous chromosome

                    // Variable Binding
                    nodeSrt& mGraphMapForChrForStart = mGraphMap[tmpChromosome][preRefStart];

                    mGraphMapForChrForStart.seqVec.push_back(preRefSeq);  // Add to graph, 0
                    mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                }

                // If the starting position of the first variant of the new chromosome is greater than 1, it means that ref information needs to be added.
                if (refStart > 1) {
                    // Add the first part of the new chromosome
                    preRefStart = 1;  // Initiation of the chromosome, 1
                    preRefEnd = refStart - 1;  // The first variant in the previous position
                    preRefLen = preRefEnd - preRefStart + 1;  // length
                    preRefSeq = mFastaMapFindIter->second.substr(0, preRefLen);  // sequence
                    
                    // Variable Binding
                    nodeSrt& mGraphMapForChrForStart = mGraphMap[chromosome][preRefStart];

                    mGraphMapForChrForStart.seqVec.push_back(preRefSeq);  // Add to graph, 0
                    mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                }
            } else {  // Otherwise the node is built from the sequence in the middle of the vcf
                string preRefSeq;
                uint32_t preRefStart;
                uint32_t preRefEnd;
                uint32_t preRefLen;

                preRefStart = tmpRefEnd + 1;
                preRefEnd = refStart - 1;
                preRefLen = preRefEnd - preRefStart + 1;

                // If it's adjacent, if the preRefLen is less than or equal to 0, skip it
                if (preRefStart <= preRefEnd) {
                    preRefSeq = mFastaMapFindIter->second.substr(preRefStart-1, preRefLen);

                    // Variable Binding
                    nodeSrt& mGraphMapForChrForStart = mGraphMap[chromosome][preRefStart];

                    mGraphMapForChrForStart.seqVec.push_back(preRefSeq);  // Add to graph, 0
                    mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                }
            }
            
            // Add a vcf node
            // Variable Binding
            nodeSrt& mGraphMapForChrForStart = mGraphMap[chromosome][refStart];

            // Add the ref node first
            mGraphMapForChrForStart.seqVec.push_back(refSeq);  // Add to graph, 0
            mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information

            // qrySeqVec
            mGraphMapForChrForStart.seqVec.insert(mGraphMapForChrForStart.seqVec.end(), qrySeqVec.begin(), qrySeqVec.end());
            // The number of bases in the graph
            mGraphBaseNum += accumulate(qrySeqVec.begin(), qrySeqVec.end(), 0, [](int sum, const string& str) {return sum + str.size();});
            // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
            if (mGraphMapForChrForStart.seqVec.size() > UINT16_MAX) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: The number of haplotypes exceeds the maximum limit of " << UINT16_MAX << "." << endl;
                exit(1);
            }

            uint16_t hapIdx = 1;  // Haplotype information index
            // Add genotyping information for each haplotype
            for (size_t i = 9; i < lineVec.size(); i++) {
                vector<string> gtVec = construct_index::gt_split(split(lineVec[i], ":")[gtIndex]);  // The corresponding typing list of the strains

                // Check if the vcfPloidy is consistent with the parameters
                if (gtVec.size() > mVcfPloidy) {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: The number of haplotypes at " << chromosome << "(" << refStart << ") exceeds the specified parameter. Excess haplotypes have been discarded." << endl;
                    gtVec.resize(mVcfPloidy);
                } else if (gtVec.size() < mVcfPloidy) {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: The number of haplotypes at " << chromosome << "(" << refStart << ") is less than the specified parameter. Filling the deficit with zeros." << endl;
                    while (gtVec.size() < mVcfPloidy) {
                        gtVec.push_back("0");
                    }
                }

                for (size_t j = 0; j < gtVec.size(); j++) {
                    if (gtVec[j] == ".") {
                        if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                            mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
                            hapIdx++;  // Haplotype index +1
                        } else {
                            cerr << "[" << __func__ << "::" << getTime() << "] "
                                << "Error: The number of haplotypes exceeds the maximum limit of " << UINT16_MAX << "." << endl;
                            exit(1);
                        }
                    } else {
                        mGraphMapForChrForStart.hapGtVec.push_back(stoul(gtVec[j]));  // haplotype information

                        if (hapIdx < UINT16_MAX) {  // Determine whether the variable is out of bounds, currently only supports 65535 haplotypes
                            hapIdx++;  // Haplotype index +1
                        } else {
                            cerr << "[" << __func__ << "::" << getTime() << "] "
                                << "Error: The number of haplotypes exceeds the maximum limit of " << UINT16_MAX << "." << endl;
                            exit(1);
                        }
                    }
                }
            }

            // update coordinates
            tmpRefStart = refStart;
            tmpRefEnd = refEnd;
            tmpChromosome = chromosome;
        }
    }

    // second half of the last chromosome
    if (tmpRefEnd < mFastaSeqMap[tmpChromosome].length()) {  // Determine whether the last variant is at the end of the chromosome
        uint32_t preRefStart = tmpRefEnd + 1;
        uint32_t preRefEnd = mFastaSeqMap[tmpChromosome].length();
        uint32_t preRefLen = preRefEnd - preRefStart + 1;

        string preRefSeq = mFastaSeqMap[tmpChromosome].substr(preRefStart-1, preRefLen);

        // Variable Binding
        nodeSrt& mGraphMapForChrForStart = mGraphMap[tmpChromosome][preRefStart];

        mGraphMapForChrForStart.seqVec.push_back(preRefSeq);  // Add to graph, 0
        mGraphMapForChrForStart.hapGtVec.push_back(0);  // haplotype information
    }

    // number of variants
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Parsed " << mSnpNum + mIndelNum + mInsNum + mDelNum + mInvNum + mDupNum + mOtherNum << " alternative alleles ...\n\n";  // print log
    cerr << "           - " << "SNP: " << mSnpNum << endl;
    cerr << "           - " << "InDels: " << mIndelNum << endl;
    cerr << "           - " << "Insertion: " << mInsNum << endl;
    cerr << "           - " << "Deletion: " << mDelNum << endl;
    cerr << "           - " << "Inversion: " << mInvNum << endl;
    cerr << "           - " << "Duplication: " << mDupNum << endl;
    cerr << "           - " << "Other: " << mOtherNum << endl << endl << endl;

    malloc_trim(0);	// 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2024/01/04
 * @version v1.0
 * @brief make mHapIdxQRmap
 * 
 * @return void
**/
void ConstructIndex::make_QRmap() {
    // Quotient and remainder for each haplotype in topHapVec
    for (uint16_t i = 0; i < mHapNum; i++) {
        mHapIdxQRmap[i] = {DIVIDE_BY_8(i), GET_LOW_3_BITS(i)};
    }
}


/**
 * @author zezhen du
 * @date 2023/09/12
 * @version v1.0
 * @brief Build the index of vcf
 * 
 * @param chromosome
 * @param refStart
 * @param refLen
 * @param lineVec     vcf split list
 * @param qrySeqVec   the sequence Vector of ALT
 * @param gtIndex     Index where gt resides
 * 
 * @return void
**/
void ConstructIndex::vcf_construct(
    const string& chromosome, 
    const uint32_t& refStart, 
    const uint32_t& refLen, 
    const vector<string>& lineVec, 
    const vector<string>& qrySeqVec, 
    const int& gtIndex
) {
    auto& emplacedValue = mVcfInfoMap[chromosome].emplace(refStart, vector<string>()).first->second;  // Store vcf file information for output

    for (const auto& qrySeq : qrySeqVec) {  // Loop over the qrySeqVec list
        // Record the amount of variation
        uint32_t qryLen = qrySeq.size();
        int32_t svLen = (int32_t)qryLen - (int32_t)refLen;
        double lengthRatio = qryLen / float(refLen);

        if (svLen == 0 && refLen == 1 && qryLen == 1) {
            mSnpNum++;
        } else if (svLen <= 49 && svLen >= -49 && refLen <= 49 && qryLen <= 49) {
            mIndelNum++;
        } else if (svLen >= -2 && svLen <= 2 && refLen > 49 && qryLen > 49) {
            mInvNum++;
        } else if (lengthRatio >= 1.8 && lengthRatio <= 2.2 && refLen > 49 && qryLen > 49) {
            mDupNum++;
        } else if (svLen < 0) {
            mDelNum++;
        } else if (svLen > 0) {
            mInsNum++;
        } else {
            mOtherNum++;
        }
    }

    // Add genotyping information for each haplotype
    for (size_t i = 0; i < lineVec.size(); i++) {
        if (i < 9) {
            emplacedValue.push_back(lineVec[i]);
            continue;
        }

        vector<string> gtVec = construct_index::gt_split(split(lineVec[i], ":")[gtIndex]);  // The corresponding typing list of the strains

        // If gt_split returns an empty list, it means the site is '.', an undefined line
        string gtTxt = "";
        
        if (gtVec.empty()) {
            for (size_t j = 0; j < mVcfPloidy; j++) {
                if (j == 0) {
                    gtTxt += "0";
                } else {
                    gtTxt += "|0";
                } 
            }
        } else {  // normal site
            // Check if the vcfPloidy is consistent with the parameters
            if (gtVec.size() >= mVcfPloidy) {
                for (size_t j = 0; j < mVcfPloidy; j++) {
                    if (j == 0) {
                        gtTxt += gtVec[j];
                    } else {
                        gtTxt += "|" + gtVec[j];
                    } 
                }
            } else if (gtVec.size() < mVcfPloidy) {
                gtTxt += join(gtVec, "|");
                for (size_t j = 0; j < (mVcfPloidy - gtVec.size()); j++) {
                    gtTxt += "|0";
                }
            }
        }

        // Add to the Vector
        emplacedValue.push_back(gtTxt);
    }
}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0
 * @brief building the k-mer index of graph
 * 
 * @return void
**/
void ConstructIndex::index() {
    // Log the start of the graph index construction on the CPU
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Initiating the construction of the graph index on CPU ...\n";

    // Save the results of multiple threads
    vector<future<tuple<map<uint32_t, nodeSrt>::iterator, unordered_map<uint64_t, vector<int8_t> >, map<uint64_t, uint8_t> > > > futureVec;  // nodeIter, KmerHapBitMap, kmerHashFreMap

    // Thread Pool
    ThreadPool pool(threads_);

    // init
    pool.init();

    // Total number of tasks
    uint64_t taskNum = 0;

    for(auto& [chromosome, startNodeMap] : mGraphMap) {  // map<chr, map<nodeStart, nodeSrt> >
        for(auto iter = startNodeMap.begin(); iter != startNodeMap.end(); iter++) {  // map<nodeStart, nodeSrt>
            if (iter->second.hapGtVec.size() == 1) {continue;}  // Skip the node if there is only one GT (0)

            // submit
            futureVec.push_back(
                pool.submit(
                    construct_index::index_run,
                    chromosome, 
                    iter, 
                    ref(startNodeMap), 
                    ref(fastMode_), 
                    ref(useUniqueKmers_), 
                    ref(mKmerLen), 
                    mbf, 
                    ref(mVcfPloidy), 
                    ref(mHapIdxQRmap)
                )
            );

            // number of tasks
            taskNum++;
        }
    }

    // Number of tasks processed
    uint64_t taskProcessNum = 0;

    // Save multithreaded results to graph
    for (auto& futureResult : futureVec) {  // vector<future<{nodeIter, KmerHapBitMap, kmerHashFreMap}> >
        auto [nodeIter, KmerHapBitMap, kmerHashFreMap] = move(futureResult.get());  // {nodeIter, KmerHapBitMap, kmerHashFreMap}

        if (KmerHapBitMap.empty()) continue;

        // Node k-mer information
        auto& kmerHashVec = nodeIter->second.kmerHashVec;

        // Bitmap corresponding to node k-mer
        for (const auto& [kmerHash, BitVec] : KmerHapBitMap) {  // unordered_map<uint64_t, vector<int8_t> > KmerHapBitMap
            kmerHashVec.push_back(kmerHash);

            auto emplaceResult = mGraphKmerHashHapStrMap.emplace(kmerHash, kmerCovFreBitVec{});
            auto& emplacedValue = emplaceResult.first->second;

            if (emplaceResult.second) {
                emplacedValue.BitVec = move(BitVec);
                // record frequency of the k-mers informations
                emplacedValue.f++;
            } else {
                auto& nodeBitVec = emplacedValue.BitVec;

                for (uint64_t i = 0; i < BitVec.size(); i++) {
                    nodeBitVec[i] |= BitVec[i];
                }
                
                if (emplacedValue.f < UINT8_MAX) {
                    emplacedValue.f++;
                }
            }
        }

        // Record the frequencies of k-mer with frequencies≥2 in bf to the graph
        for (const auto& [kmerHash, frequency] : kmerHashFreMap) {  // map<uint64_t, uint8_t>
            auto findIter = mGraphKmerHashHapStrMap.find(kmerHash);
            if (findIter == mGraphKmerHashHapStrMap.end()) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: The k-mer hash '" << kmerHash << "' is not found in the mGraphKmerHashHapStrMap." << endl;
                exit(1);
            }
            auto& f = findIter->second.f;
            if (f == 1) {
                f += frequency - 1;
            }
        }

        // tasks processed
        taskProcessNum++;

        // Print log every 5% of tasks
        if (taskProcessNum > 0 && taskProcessNum % (taskNum / 20) == 0) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << std::fixed << std::setprecision(0) << "Indexing progress: " << std::setw(3) << static_cast<double>(taskProcessNum) / taskNum * 100.0 << "%\n";
        }
    }

    // free memory
    vector<future<tuple<map<uint32_t, nodeSrt>::iterator, unordered_map<uint64_t, vector<int8_t> >, map<uint64_t, uint8_t> > > >().swap(futureVec);

    malloc_trim(0);	// 0 is for heap memory

    // Close the thread pool
    pool.shutdown();
}


/**
 * @author zezhen du
 * @date 2023/08/30
 * @version v1.0
 * @brief Merge k-mer information from Genome Graph into nodes.
 * 
 * @return void
**/
void ConstructIndex::graph2node() {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Merging k-mer information from Genome Graph into Nodes ...\n\n\n";  // print log

    // Save the results of multiple threads
    vector<future<int> > futureVec;  // vector<future<int> >

    // Thread Pool
    ThreadPool pool(threads_);

    // init
    pool.init();

    for (auto& [chromosome, startNodeMap] : mGraphMap) {  // map<chr, map<nodeStart, nodeSrt> >
        for (auto& [nodeStart, node] : startNodeMap) {  // map<nodeStart, nodeSrt>
        
            if (node.hapGtVec.size() == 1) {continue;}  // If only 0, skip the node

            // submit
            futureVec.push_back(
                pool.submit(
                    construct_index::graph2node_run,
                    ref(node.kmerHashVec), 
                    ref(node.GraphKmerHashHapStrMapIterVec), 
                    ref(mGraphKmerHashHapStrMap)
                )
            );
        }
    }

    // Save multithreaded results to graph
    for (auto& futureResult : futureVec) {  // vector<future<int> >
        auto result = move(futureResult.get());
    }

    // free memory
    vector<future<int> >().swap(futureVec);

    malloc_trim(0);  // 0 is for heap memory

    // Close the thread pool
    pool.shutdown();
}


/**
 * @author zezhen du
 * @brief Save the graph index to file
 * 
 * @return void
**/
void ConstructIndex::save_index() {
    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Genome Graph index saved to file: " << outputGraphFileName_ << endl;

    std::ofstream outFile(outputGraphFileName_, std::ios::binary);
    if (!outFile.is_open()) {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << outputGraphFileName_ << "': No such file or directory." << endl;
        exit(1);
    }

    /* -------------------------------------------------- mGraphBaseNum, mKmerLen and mVcfPloidy -------------------------------------------------- */
    // Write the total number of bases in the Genome Graph
    outFile.write(reinterpret_cast<const char*>(&mGraphBaseNum), sizeof(uint64_t));

    // Write the length of k-mer and mVcfPloidy
    outFile.write(reinterpret_cast<const char*>(&mKmerLen), sizeof(uint32_t));
    outFile.write(reinterpret_cast<const char*>(&mVcfPloidy), sizeof(uint32_t));

    /* -------------------------------------------------- mVcfHead and mVcfInfoMap  -------------------------------------------------- */
    // Write mVcfHead
    uint32_t vcfHeadLength = mVcfHead.length();
    outFile.write(reinterpret_cast<const char*>(&vcfHeadLength), sizeof(uint32_t));
    outFile.write(mVcfHead.data(), vcfHeadLength);

    // Write mVcfInfoMap
    uint32_t mVcfInfoMapSize = mVcfInfoMap.size();
    outFile.write(reinterpret_cast<const char*>(&mVcfInfoMapSize), sizeof(uint32_t));

    for (const auto& chrInfo : mVcfInfoMap) {
        const std::string& chrName = chrInfo.first;
        uint32_t chrNameLength = chrName.length();
        outFile.write(reinterpret_cast<const char*>(&chrNameLength), sizeof(uint32_t));
        outFile.write(chrName.data(), chrNameLength);

        // chromosome length
        uint32_t chrLen = mFastaLenMap.at(chrName);
        outFile.write(reinterpret_cast<const char*>(&chrLen), sizeof(uint32_t));

        const std::map<uint32_t, std::vector<std::string> >& startVecMap = chrInfo.second;
        uint32_t startVecMapSize = startVecMap.size();
        outFile.write(reinterpret_cast<const char*>(&startVecMapSize), sizeof(uint32_t));

        for (const auto& startVec : startVecMap) {
            uint32_t start = startVec.first;
            const std::vector<std::string>& infoVec = startVec.second;
            uint32_t infoVecSize = infoVec.size();
            outFile.write(reinterpret_cast<const char*>(&start), sizeof(uint32_t));
            outFile.write(reinterpret_cast<const char*>(&infoVecSize), sizeof(uint32_t));

            for (const auto& info : infoVec) {
                uint32_t infoLength = info.length();
                outFile.write(reinterpret_cast<const char*>(&infoLength), sizeof(uint32_t));
                outFile.write(info.data(), infoLength);
            }
        }
    }

    /* -------------------------------------------------- mHapNum and mHapMap  -------------------------------------------------- */
    outFile.write(reinterpret_cast<const char*>(&mHapNum), sizeof(uint16_t));

    for (const auto& hapEntry : mHapMap) {
        uint16_t hapIdx = hapEntry.first;
        const std::string& hapName = hapEntry.second;

        outFile.write(reinterpret_cast<const char*>(&hapIdx), sizeof(uint16_t));
        uint32_t hapNameLength = hapName.length();
        outFile.write(reinterpret_cast<const char*>(&hapNameLength), sizeof(uint32_t));
        outFile.write(reinterpret_cast<const char*>(hapName.data()), hapNameLength);
    }

    /* -------------------------------------------------- mGraphMap  -------------------------------------------------- */
    // Write the size of the outer map.
    uint32_t ChrNum = mGraphMap.size();
    outFile.write(reinterpret_cast<const char*>(&ChrNum), sizeof(uint32_t));

    for (const auto& chrNodes : mGraphMap) {  // map<chr, map<nodeStart, nodeSrt> >
        // Write the chromosome name and number of nodes in this chromosome.
        const std::string& chrName = chrNodes.first;
        uint32_t numNodes = chrNodes.second.size();

        uint32_t chrNameLength = chrName.length();
        outFile.write(reinterpret_cast<const char*>(&chrNameLength), sizeof(uint32_t));
        outFile.write(chrName.data(), chrNameLength);
        outFile.write(reinterpret_cast<const char*>(&numNodes), sizeof(uint32_t));

        for (const auto& nodeData : chrNodes.second) {  // map<nodeStart, nodeSrt>
            // Write the node start position and its data.
            uint32_t nodeStartPos = nodeData.first;
            const nodeSrt& nodeSrtData = nodeData.second;

            outFile.write(reinterpret_cast<const char*>(&nodeStartPos), sizeof(uint32_t));

            // Write the sequence information.
            uint32_t numSeqEntries = nodeSrtData.seqVec.size();
            outFile.write(reinterpret_cast<const char*>(&numSeqEntries), sizeof(uint32_t));

            for (const auto& seq : nodeSrtData.seqVec) {
                uint32_t seqLen = seq.length();
                outFile.write(reinterpret_cast<const char*>(&seqLen), sizeof(uint32_t));
                outFile.write(reinterpret_cast<const char*>(seq.data()), seqLen);
            }

            // Write the genotype information.
            uint32_t numGtEntries = nodeSrtData.hapGtVec.size();
            outFile.write(reinterpret_cast<const char*>(&numGtEntries), sizeof(uint32_t));
            outFile.write(reinterpret_cast<const char*>(&nodeSrtData.hapGtVec[0]), sizeof(uint16_t) * numGtEntries);

            // Write the k-mer hash information.
            uint32_t numKmerHashEntries = nodeSrtData.kmerHashVec.size();
            outFile.write(reinterpret_cast<const char*>(&numKmerHashEntries), sizeof(uint32_t));
            outFile.write(reinterpret_cast<const char*>(&nodeSrtData.kmerHashVec[0]), sizeof(uint64_t) * numKmerHashEntries);
        }
    }

    /* -------------------------------------------------- mGraphKmerHashHapStrMap  -------------------------------------------------- */
    // Read base
    uint64_t ReadBase = 0;
    outFile.write(reinterpret_cast<const char*>(&ReadBase), sizeof(uint64_t));

    for (const auto& kvp : mGraphKmerHashHapStrMap) {
        const uint64_t kmerHash = kvp.first;
        const kmerCovFreBitVec& kmerCovFreBitVecStr = kvp.second;
        const vector<int8_t>& BitVec = kmerCovFreBitVecStr.BitVec;
        uint64_t BitVecLen = BitVec.size();

        // Write kmerHash
        outFile.write(reinterpret_cast<const char*>(&kmerHash), sizeof(uint64_t));

        // Write kmerCov and kmerFre
        outFile.write(reinterpret_cast<const char*>(&kmerCovFreBitVecStr.c), sizeof(uint8_t));
        outFile.write(reinterpret_cast<const char*>(&kmerCovFreBitVecStr.f), sizeof(uint8_t));

        // Write BitVec
        outFile.write(reinterpret_cast<const char*>(&BitVecLen), sizeof(uint64_t));
        for (const int8_t& Bit : BitVec) {
            outFile.write(reinterpret_cast<const char*>(&Bit), sizeof(int8_t));
        }
    }

    // Close the file
    outFile.close();
}


/**
 * @author zezhen du
 * @brief load the graph index from file
 * 
 * @return void
**/
void ConstructIndex::load_index() {
    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Genome Graph index loaded from file: " << inputGraphFileName_ << endl;

    std::ifstream inFile(inputGraphFileName_, std::ios::binary);
    if (!inFile.is_open()) {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << inputGraphFileName_ << "': No such file or directory." << endl;
        exit(1);
    }

    /* -------------------------------------------------- mKmerLen and mVcfPloidy -------------------------------------------------- */
    // Load the total number of bases in the Genome Graph
    mGraphBaseNum = 0;
    inFile.read(reinterpret_cast<char*>(&mGraphBaseNum), sizeof(uint64_t));

    // load the length of k-mer and mVcfPloidy
    mKmerLen = 0;
    inFile.read(reinterpret_cast<char*>(&mKmerLen), sizeof(uint32_t));
    mVcfPloidy = 0;
    inFile.read(reinterpret_cast<char*>(&mVcfPloidy), sizeof(uint32_t));

    /* -------------------------------------------------- mVcfHead and mVcfInfoMap  -------------------------------------------------- */
    // Clear the existing data
    mVcfHead.clear();
    mFastaLenMap.clear();
    mVcfInfoMap.clear();

    // Read mVcfHead
    uint32_t vcfHeadLength;
    inFile.read(reinterpret_cast<char*>(&vcfHeadLength), sizeof(uint32_t));
    mVcfHead.resize(vcfHeadLength);
    inFile.read(&mVcfHead[0], vcfHeadLength);

    // Read mVcfInfoMap
    uint32_t mVcfInfoMapSize;
    inFile.read(reinterpret_cast<char*>(&mVcfInfoMapSize), sizeof(uint32_t));

    for (uint32_t i = 0; i < mVcfInfoMapSize; ++i) {
        uint32_t chrNameLength;
        inFile.read(reinterpret_cast<char*>(&chrNameLength), sizeof(uint32_t));
        std::string chrName;
        chrName.resize(chrNameLength);
        inFile.read(&chrName[0], chrNameLength);

        // chromosome length
        uint32_t chrLen;
        inFile.read(reinterpret_cast<char*>(&chrLen), sizeof(uint32_t));
        mFastaLenMap[chrName] = chrLen;
        mGenomeSize += chrLen;

        uint32_t startVecMapSize;
        inFile.read(reinterpret_cast<char*>(&startVecMapSize), sizeof(uint32_t));

        std::map<uint32_t, std::vector<std::string> > startVecMap;
        for (uint32_t j = 0; j < startVecMapSize; ++j) {
            uint32_t start;
            inFile.read(reinterpret_cast<char*>(&start), sizeof(uint32_t));

            uint32_t infoVecSize;
            inFile.read(reinterpret_cast<char*>(&infoVecSize), sizeof(uint32_t));

            std::vector<std::string> infoVec;
            for (uint32_t k = 0; k < infoVecSize; ++k) {
                uint32_t infoLength;
                inFile.read(reinterpret_cast<char*>(&infoLength), sizeof(uint32_t));
                std::string info;
                info.resize(infoLength);
                inFile.read(&info[0], infoLength);
                infoVec.push_back(std::move(info));
            }
            startVecMap[start] = infoVec;
        }
        mVcfInfoMap[chrName] = startVecMap;
    }

    /* -------------------------------------------------- mHapNum and mHapMap  -------------------------------------------------- */
    inFile.read(reinterpret_cast<char*>(&mHapNum), sizeof(uint16_t));

    // Clear the existing data
    mHapMap.clear();

    for (uint16_t i = 0; i < mHapNum; ++i) {
        uint16_t hapIdx;
        uint32_t hapNameLength;

        inFile.read(reinterpret_cast<char*>(&hapIdx), sizeof(uint16_t));
        inFile.read(reinterpret_cast<char*>(&hapNameLength), sizeof(uint32_t));

        std::string hapName;
        hapName.resize(hapNameLength);
        inFile.read(&hapName[0], hapNameLength);

        // Add the haplotype to the map.
        mHapMap.emplace(hapIdx, std::move(hapName));
    }

    /* -------------------------------------------------- mGraphMap  -------------------------------------------------- */
    // Clear the existing data
    mGraphMap.clear();

    // Read the size of the outer map.
    uint32_t chrNum;
    inFile.read(reinterpret_cast<char*>(&chrNum), sizeof(uint32_t));

    for (uint32_t i = 0; i < chrNum; ++i) {
        // Read the chromosome name and number of nodes in this chromosome.
        uint32_t chrNameLength;
        inFile.read(reinterpret_cast<char*>(&chrNameLength), sizeof(uint32_t));

        std::string chrName;
        chrName.resize(chrNameLength);
        inFile.read(&chrName[0], chrNameLength);

        uint32_t numNodes;
        inFile.read(reinterpret_cast<char*>(&numNodes), sizeof(uint32_t));

        for (uint32_t j = 0; j < numNodes; ++j) {
            // Read the node start position and its data.
            uint32_t nodeStartPos;
            nodeSrt nodeSrtData;

            inFile.read(reinterpret_cast<char*>(&nodeStartPos), sizeof(uint32_t));

            // Read the sequence information.
            uint32_t numSeqEntries;
            inFile.read(reinterpret_cast<char*>(&numSeqEntries), sizeof(uint32_t));

            for (uint32_t k = 0; k < numSeqEntries; ++k) {
                uint32_t seqLen;
                inFile.read(reinterpret_cast<char*>(&seqLen), sizeof(uint32_t));

                std::string seq;
                seq.resize(seqLen);
                inFile.read(&seq[0], seqLen);
                nodeSrtData.seqVec.push_back(seq);
            }

            // Read the genotype information.
            uint32_t numGtEntries;
            inFile.read(reinterpret_cast<char*>(&numGtEntries), sizeof(uint32_t));
            nodeSrtData.hapGtVec.resize(numGtEntries);
            inFile.read(reinterpret_cast<char*>(&nodeSrtData.hapGtVec[0]), sizeof(uint16_t) * numGtEntries);

            // Read the k-mer hash information.
            uint32_t numKmerHashEntries;
            inFile.read(reinterpret_cast<char*>(&numKmerHashEntries), sizeof(uint32_t));
            nodeSrtData.kmerHashVec.resize(numKmerHashEntries);
            inFile.read(reinterpret_cast<char*>(&nodeSrtData.kmerHashVec[0]), sizeof(uint64_t) * numKmerHashEntries);

            mGraphMap[chrName][nodeStartPos] = move(nodeSrtData);
        }
    }

    /* -------------------------------------------------- mGraphKmerHashHapStrMap  -------------------------------------------------- */
    // Clear the existing data in mGraphKmerHashHapStrMap
    mGraphKmerHashHapStrMap.clear();

    uint64_t kmerHash;
    kmerCovFreBitVec kmerCovFreBitVecStr;
    uint8_t c;
    uint8_t f;
    uint64_t BitVecLen;

    // Read base
    uint64_t ReadBase = 0;
    // Read base
    inFile.read(reinterpret_cast<char*>(&ReadBase), sizeof(uint64_t));

    // Read kmerHash, kmerCov, and kmerFre from the file
    while (inFile.read(reinterpret_cast<char*>(&kmerHash), sizeof(uint64_t))) {

        // Read kmerCov and kmerFre
        inFile.read(reinterpret_cast<char*>(&c), sizeof(uint8_t));
        inFile.read(reinterpret_cast<char*>(&f), sizeof(uint8_t));

        // Read BitVec length
        inFile.read(reinterpret_cast<char*>(&BitVecLen), sizeof(uint64_t));

        // Read BitVec
        vector<int8_t> BitVec(BitVecLen);
        for (auto& Bit : BitVec) {
            inFile.read(reinterpret_cast<char*>(&Bit), sizeof(int8_t));
        }

        // Store data in the map
        kmerCovFreBitVecStr.c = c;
        kmerCovFreBitVecStr.f = f;
        kmerCovFreBitVecStr.BitVec = BitVec;
        mGraphKmerHashHapStrMap[kmerHash] = kmerCovFreBitVecStr;
    }

    // Close the file
    inFile.close();
}


/**
 * @brief graph index for k-mer (threads)
 * 
 * @date 2024/04/16
 * 
 * @param chromosome            mGraphMap output by construct map<chr, map<start, nodeSrt> >
 * @param nodeIter              node iterator
 * @param startNodeMap          Chromosome all nodes
 * @param fastMode              fast mode
 * @param useUniqueKmers        use unique k-mers for indexing
 * @param kmerLen               the length of kmer
 * @param bf                    Kmer frequency in the reference genome: Counting Bloom Filter
 * @param vcfPloidy             ploidy of genotypes in VCF file
 * @param hapIdxQRmap           map<hapIdx, tuple<quotient, remainder> >
 * 
 * @return {nodeIter, tmpKmerHapBitMap, kmerHashFreMap}     kmer: map<kmerHash, vector<int8_t> >
**/
tuple<map<uint32_t, nodeSrt>::iterator, unordered_map<uint64_t, vector<int8_t> >, map<uint64_t, uint8_t> > construct_index::index_run(
    string chromosome, 
    map<uint32_t, nodeSrt>::iterator nodeIter, 
    const map<uint32_t, nodeSrt>& startNodeMap, 
    const bool& fastMode, 
    const bool& useUniqueKmers,
    const uint32_t& kmerLen, 
    BloomFilter* bf, 
    const uint32_t& vcfPloidy, 
    const unordered_map<uint16_t, tuple<uint16_t, uint16_t> >& hapIdxQRmap
) {
    // Minimum k-mer frequency
    uint8_t MIN_KMER_FRE = 1;
    unordered_map<uint64_t, pair<vector<int8_t>, uint8_t > > KmerHapBitFrePairMap;  // map<kmerHash, pair<vector<int8_t>, frequence > >

    unordered_map<uint64_t, vector<int8_t> > KmerHapBitMap;  // kmer: map<kmerHash, vector<int8_t> >:  0000 0000, Each bits represents a haplotype, 0->False 1->True

    map<uint64_t, uint8_t> kmerHashFreMap;  // map<kmerHash, frequence≥2>

    const auto& seqVec = nodeIter->second.seqVec;  // vector<seq>
    const auto& hapGtVec = nodeIter->second.hapGtVec;  // vector<gt>

    uint16_t haplotype = 0;  // Index of the haplotype
    
    for (const auto& gt : hapGtVec) {  // Iterate over the genotypes
        // fast mode
        // calculate whether the genotype of the corresponding sample is empty or zero. If it is, skip all its haplotypes
        if (fastMode && haplotype > 0 && gt == (uint16_t)0) {
            uint16_t groupIdx = (haplotype - 1) / vcfPloidy;
            uint16_t hapIdxL = groupIdx * vcfPloidy + 1;
            uint16_t hapIdxR = (groupIdx + 1) * vcfPloidy;

            uint16_t gtSum = std::accumulate(
                hapGtVec.begin() + hapIdxL,
                hapGtVec.begin() + hapIdxR + 1,
                0
            );
            
            if (gtSum == 0) {
                // Haplotype index increment
                ++haplotype;
                continue;
            }
        }
        
        // Check if the array is out of range
        if (gt >= seqVec.size()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Error: The node '" << chromosome << "-" << nodeIter->first << "' lacks sequence information for haplotype " << gt << "." << endl;
            exit(1);
        }

        // Get ALT sequence
        string seqTmp = seqVec[gt];

        // Sequence with upstream and downstream 1 k-mer sequence
        pair<string, string> upDownSeq = construct_index::find_node_up_down_seq(
            haplotype, 
            gt, 
            seqTmp, 
            kmerLen - 1, 
            nodeIter, 
            startNodeMap
        );
        if (debugConstruct) {
            cerr << "Node Start:" << nodeIter->first << ", Haplotype:" << haplotype << ", GT:" << +gt << ", Upstream:" << upDownSeq.first << ", Current:" << seqTmp << ", Downstream:" << upDownSeq.second << endl << endl;
        }
        
        seqTmp = upDownSeq.first + seqTmp + upDownSeq.second;  // Add upstream and downstream sequences

        // k-mer indexing
        map<uint8_t, unordered_set<uint64_t> > freKmerHashSetMap;  // map<frequency, unordered_set<kmerHash> >
        kmerBit::kmer_sketch_construct(seqTmp, kmerLen, freKmerHashSetMap, bf);

        // Variable to store the quotient and remainder
        const uint16_t& quotient = get<0>(hapIdxQRmap.at(haplotype));
        const uint16_t& remainder = get<1>(hapIdxQRmap.at(haplotype));

        // Record all k-mers to KmerHapBitFrePairMap
        for (const auto& [frequency, kmerHashSet] : freKmerHashSetMap) {
            for (const auto& kmerHash : kmerHashSet) {
                auto emplacedValue = KmerHapBitFrePairMap.emplace(kmerHash, make_pair(vector<int8_t>(DIVIDE_BY_8(hapGtVec.size()) + 1, 0), frequency)).first;
                emplacedValue->second.second = frequency;
                int8_t& hapBitTmp = emplacedValue->second.first[quotient];
                construct_index::set_bit_to_one(hapBitTmp, remainder);  // Set the corresponding haplotype to 1

                // If the k-mer is present in the MBF, but haplotype 0 of that node does not contain it, set the last bit of hapBit to 1, which means that the frequency in the genome is 1.
                if (gt != (uint16_t)0 && bf->find(kmerHash) && construct_index::get_bit(emplacedValue->second.first[0], 0) == 0) {
                    // set the last bit of hapBit to 1
                    construct_index::set_bit_to_one(emplacedValue->second.first.back(), 7);
                }
            }

            // Update MIN_KMER_FRE if the frequency is greater than or equal to 2
            if (frequency >= 2 && !useUniqueKmers) {
                MIN_KMER_FRE = min(MIN_KMER_FRE, frequency);
            }
        }

        // clear memory (freKmerHashSetMap)
        map<uint8_t, unordered_set<uint64_t> >().swap(freKmerHashSetMap);

        // Haplotype index increment
        ++haplotype;
    }

    // Transfer k-mers with frequency <= MIN_KMER_FRE to KmerHapBitMap
    for (auto it = KmerHapBitFrePairMap.begin(); it != KmerHapBitFrePairMap.end(); ) {
        if (it->second.second <= MIN_KMER_FRE) {
            KmerHapBitMap.emplace(it->first, move(it->second.first));
            // Record the frequency of k-mer with frequencies≥2
            if (it->second.second >= 2) {
                kmerHashFreMap.emplace(it->first, it->second.second);
            }
            it = KmerHapBitFrePairMap.erase(it);
        } else {
            ++it;
        }
    }

    return {nodeIter, move(KmerHapBitMap), move(kmerHashFreMap)};
}


/**
 * @author zezhen du
 * @date 2024/02/20
 * @version v1.0.2
 * @brief find the sequence information corresponding to the haplotype upstream and downstream of the node
 * 
 * @param haplotype       haplotype index
 * @param altGt           genotype
 * @param altSeq          ALT sequence
 * @param seqLen          k-mer length - 1
 * @param nodeIter        startNodeMap iterator where the node is located
 * @param startNodeMap    All chromosome node information startNodeMap, index function construction
 * 
 * @return pair<upSeq, downSeq>
**/
pair<string, string> construct_index::find_node_up_down_seq(
    const uint16_t & haplotype, 
    const uint16_t& altGt, 
    string& altSeq, 
    const uint32_t & seqLen,
    const map<uint32_t, nodeSrt>::iterator & nodeIter, 
    const map<uint32_t, nodeSrt> & startNodeMap
) {
    // start and end position
    uint32_t altStart = nodeIter->first;
    uint32_t altEnd = static_cast<uint32_t>(altStart + nodeIter->second.seqVec[0].size() - 1);
    uint32_t altLen = altSeq.size();

    // Upstream/Downstream sequence
    string upSeq;
    string downSeq;

    // Record the previous ALT information
    vector<uint32_t> preQryLenVec = {altLen};
    vector<uint16_t> preGtVec = {altGt};
    vector<uint32_t> preNodeStartVec = {altStart};
    vector<uint32_t> preNodeEndVec = {altEnd};
    uint32_t preQryLen = altLen;
    uint16_t preGt = altGt;

    // Temporary iterator
    auto iterTmp = nodeIter;

    // Add the end sequence of the previous node
    while (upSeq.size() < seqLen && iterTmp != startNodeMap.begin()) {
        // iterator advance
        --iterTmp;

        const uint32_t& nodeStartTmp = iterTmp->first;  // start position
        const nodeSrt& nodeTmp = iterTmp->second;  // node information
        uint32_t nodeEndTmp = nodeStartTmp + nodeTmp.seqVec.at(0).size() - 1;  // end position

        uint16_t gtTmp = (haplotype < nodeTmp.hapGtVec.size()) ? nodeTmp.hapGtVec[haplotype] : 0;

        // Check if the array is out of bounds
        if (gtTmp >= nodeTmp.seqVec.size()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Error: The node '" << altStart << "' lacks sequence information for haplotype " << gtTmp << "." << endl;
            exit(1);
        }
        // Get ALT sequence
        string seqTmp = nodeTmp.seqVec[gtTmp];

        /* The current version of the nested graph only has coordinates for "ref." Therefore, when nodes intersect and the genotype is 0, the sequence are truncated. */
        /*
         1.
            |_______________|
         |______|

         2.
                |_______________|
         |______|
        */
        while (preNodeStartVec.size() > 0 && nodeEndTmp >= preNodeStartVec.back() && !seqTmp.empty()) {
            if (gtTmp == 0) {  // If genotype 0, truncate the sequence
                seqTmp = seqTmp.substr(0, preNodeStartVec.back() - nodeStartTmp);
                break;
            } else if (preGtVec.back() == 0 && !upSeq.empty()) {
                /* 
                 Check where the coordinates overlap to see if there is a sequence of the genotype, and if so, replace it
                 Node1: 63124  DEL3  CT  C  0/0
                 Node2: 63125  SNP134  T  A  1/1 
                 In the given coordinate situation, at the previous node 
                 the corresponding genotype was 0, so the sequence corresponded 
                 to CT. However, at the next node 63125 actually corresponds to 
                 the sequence A. The following code solves this problem.
                */

                // query sequence
                uint32_t preQryLenTmp = min(nodeEndTmp - preNodeStartVec.back() + 1, preQryLenVec.back());
                upSeq = upSeq.substr(preQryLenTmp, upSeq.size() - preQryLenTmp);

                // delete from vector
                preQryLenVec.pop_back();
                preGtVec.pop_back();
                preNodeStartVec.pop_back();
                preNodeEndVec.pop_back();

                continue;
            }
            break;
        }

        if (seqTmp.empty()) {continue;};
        
        // Update coordinates
        preNodeStartVec.push_back(nodeStartTmp);
        preNodeEndVec.push_back(nodeEndTmp);

        if (debugConstruct) {
            cerr << "UP - Start:" << nodeStartTmp << ", GT:" << +gtTmp << ", sequence:" << seqTmp << endl;
        }

        int64_t remainingLen = seqLen - upSeq.size();
        if (seqTmp.size() >= remainingLen) {
            upSeq.insert(0, seqTmp.substr(seqTmp.size() - remainingLen, remainingLen));
            preQryLen = remainingLen;
            preQryLenVec.push_back(preQryLen);
        } else {
            upSeq.insert(0, seqTmp);
            preQryLen = seqTmp.size();
            preQryLenVec.push_back(preQryLen);
        }

        preGt = gtTmp;
        preGtVec.push_back(preGt);
    }

    // Reset the iterator and coordinates
    iterTmp = nodeIter;
    preQryLen = altLen;
    preGt = altGt;
    preQryLenVec = {altLen};
    preGtVec = {altGt};
    preNodeStartVec = {altStart};
    preNodeEndVec = {altEnd};

    // Iterate over the next node's start sequence
    while ((downSeq.size() < seqLen) && (++iterTmp != startNodeMap.end())) {
        const uint32_t& nodeStartTmp = iterTmp->first;  // start position
        const nodeSrt& nodeTmp = iterTmp->second;  // node information
        uint32_t nodeLenTmp = nodeTmp.seqVec[0].size();  // node length
        uint32_t nodeEndTmp = nodeStartTmp + nodeLenTmp - 1;  // end position

        uint16_t gtTmp = (haplotype < nodeTmp.hapGtVec.size()) ? nodeTmp.hapGtVec[haplotype] : 0;  // the genotype in this node

        // Check if the array is out of bounds
        if (gtTmp >= nodeTmp.seqVec.size()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Error: The node '" << altStart << "' lacks sequence information for haplotype " << gtTmp << "." << endl;
            exit(1);
        }
        // Get ALT sequence
        string seqTmp = nodeTmp.seqVec[gtTmp];

        /*
         If the original node represents a deletion or insertion, and the genotype of the haplotype is either 0 (deletion) or 1 (insertion), 
         but the haplotype also contains a SNP on the deletion or insertion, then the deletion corresponding to the haplotype should be replaced. 
         The corresponding position is the SNP. For example:

         chr1    17008   Deletion        TTTTTTT       T       0/1
         chr1    17009   SNP     T       A       1/1
         1.
            |_______________|
                |______|
         2.
            |_______________|
                     |______|
         3.
            |_______________|
                |_______________|
         4.
            |_______________|
                            |_______________|

         In the given coordinate situation, the sequence of haplotype 1 corresponds to TTTTTTT. 
         However, in the subsequent node, 17009 actually corresponds to 'A', hence the sequence should be 'TATTTTT'. 
         The following code addresses this issue.
        */
        // 1 and 2. The end of the next node is less than the start of the current node (SNP).
        if (altGt == 0 && gtTmp != 0 && nodeEndTmp <= altEnd && seqTmp.size() == 1 && nodeLenTmp == 1) {
            if (debugConstruct) {
                cerr << "REPLACE - " << "altSeq.size:" << altSeq.size() << ", altStart:" << altStart << ", altEnd:" << altEnd 
                    << ", nodeStartTmp:" << nodeStartTmp << ", nodeEndTmp:" << nodeEndTmp 
                    << ", replaceStartTmp:" << nodeStartTmp - altStart << ", replaceLenTmp:" << nodeLenTmp << ", seqTmp:" << seqTmp << endl;
            }

            altSeq.replace(nodeStartTmp - altStart, nodeLenTmp, seqTmp);
        }
        // 3 and 4. The end of the next node is greater than the end of the current node.

        // Check if nodes intersect
        if (nodeEndTmp <= altEnd) continue;

        /* The current version of the nested graph only has coordinates for "ref." Therefore, when nodes intersect and the genotype is 0, the sequence are truncated. */
        // Check if nodes intersect
        /*
         |_______________|
               |______|
        */
        while (preNodeEndVec.size() > 0 && nodeEndTmp <= preNodeEndVec.back() && !seqTmp.empty()) {  // If the end position of the next node is less than the current node, proceed to the next iteration directly.
            if (gtTmp == 0) {
                seqTmp = "";
                break;
            } else if (preGt == 0 && !downSeq.empty()) {
                /* 
                 Check where the coordinates overlap to see if there is a sequence of the genotype, and if so, replace it
                 Node1: 63124  DEL3  CT  C  0/0
                 Node2: 63125  SNP134  T  A  1/1 
                 In the given coordinate situation, at the previous node 
                 the corresponding genotype was 0, so the sequence corresponded 
                 to CT. However, at the next node 63125 actually corresponds to 
                 the sequence A. The following code solves this problem.
                */

                uint32_t preQryLenTmp = min(preNodeEndVec.back() - nodeStartTmp + 1, preQryLenVec.back());
                downSeq = downSeq.substr(0, downSeq.size() - preQryLenTmp);

                // delete from vector
                preQryLenVec.pop_back();
                preGtVec.pop_back();
                preNodeStartVec.pop_back();
                preNodeEndVec.pop_back();

                continue;
            }
            break;
        }

        /*
         1.
         |_______________|
               |_______________|

         2.
         |_______________|
                         |_______________|
        */
        while (preNodeEndVec.size() > 0 && nodeStartTmp <= preNodeEndVec.back() && !seqTmp.empty()) {  // The end of this node is greater than the start of the next node.
            if (gtTmp == 0) {  // When referring to the sequence, update the sequence.
                seqTmp = seqTmp.substr(preNodeEndVec.back() - nodeStartTmp + 1, nodeEndTmp - preNodeEndVec.back());  // Truncate the sequence of ref.
                break;
            } else if (preGt == 0 && !downSeq.empty()) {
                /* 
                 This section of the code checks for overlapping coordinates to determine if there is a sequence corresponding to the genotype. 
                 If such a sequence exists, it replaces it. Here's an example to illustrate this:

                 Consider two nodes:
                 Node1: At position 63124, there's a deletion (DEL3) changing 'CT' to 'C'. The genotype is 0/0.
                 Node2: At position 63125, there's a single nucleotide polymorphism (SNP134) changing 'T' to 'A'. The genotype is 1/1.

                 In this scenario, the genotype at Node1 is 0, so the sequence corresponds to 'CT'. 
                 However, at Node2 (position 63125), the sequence actually corresponds to 'A', not 'T'. 

                 The code below addresses this discrepancy by replacing the sequence at the overlapping coordinates, ensuring the sequence accurately reflects the genotype at each position.
                */

                // query sequence
                uint32_t preQryLenTmp = min(preNodeEndVec.back() - nodeStartTmp + 1, preQryLenVec.back());
                downSeq = downSeq.substr(0, downSeq.size() - preQryLenTmp);

                // delete from vector
                preQryLenVec.pop_back();
                preGtVec.pop_back();
                preNodeStartVec.pop_back();
                preNodeEndVec.pop_back();

                continue;
            }
            break;
        }
       
       if (seqTmp.empty()) {continue;};

        // update coordinates
        preNodeStartVec.push_back(nodeStartTmp);
        preNodeEndVec.push_back(nodeEndTmp);

        if (debugConstruct) {
            cerr << "DOWN - Start:" << nodeStartTmp << ", GT:" << +gtTmp << ", sequence:" << seqTmp << endl;
        }

        int64_t remainingLen = seqLen - downSeq.size();
        if (seqTmp.size() >= remainingLen) {
            downSeq.append(seqTmp, 0, remainingLen);
            preQryLen = remainingLen;
            preQryLenVec.push_back(preQryLen);
        } else {
            downSeq.append(seqTmp);
            preQryLen = seqTmp.size();
            preQryLenVec.push_back(preQryLen);
        }

        preGt = gtTmp;
        preGtVec.push_back(preGt);
    }

    return make_pair(upSeq, downSeq);
}


// sort by k-mer frequency in ascending order
bool construct_index::compare_frequency(
    const unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator& a, 
    const unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator& b
) {
    return a->second.f < b->second.f;
}


/**
 * @brief Merge k-mer information from Genome Graph into nodes. (threads)
 * 
 * @date 2023/08/30
 * 
 * @param kmerHashVec                     Node k-mers hash
 * @param GraphKmerHashHapStrMapIterVec   Iterator pointing to mGraphKmerHashHapStrMap, vector<iter>
 * @param GraphKmerHashHapStrMap          Total k-mers coverage and frequency information
 * 
 * @return 0
**/
int construct_index::graph2node_run(
    vector<uint64_t>& kmerHashVec, 
    vector<unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator>& GraphKmerHashHapStrMapIterVec, 
    const unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap
) {
    vector<unordered_map<uint64_t, kmerCovFreBitVec>::const_iterator> GraphKmerHashHapStrMapIterVecTmp;  // Iterator pointing to mGraphKmerHashHapStrMap, vector<iter>

    // If there are few node-specific k-mers, they will not be filtered (SNPs, Indels).
    // uint16_t frequenceNum = (kmerHashVec.size() > 100) ? 1 : UINT16_MAX;
    uint16_t frequenceNum = UINT16_MAX;

    for (const auto& kmerHash: kmerHashVec) {
        auto findIter = GraphKmerHashHapStrMap.find(kmerHash);  // The frequency of the k-mer in the graph
        if (findIter == GraphKmerHashHapStrMap.end() || findIter->second.f > frequenceNum) {
            continue;  // delete from node
        } else {
            GraphKmerHashHapStrMapIterVecTmp.push_back(findIter);  // save to node
        }
    }

    // Keep the first 128 k-mers
    if (GraphKmerHashHapStrMapIterVecTmp.size() > 128) {
        std::sort(GraphKmerHashHapStrMapIterVecTmp.begin(), GraphKmerHashHapStrMapIterVecTmp.end(), construct_index::compare_frequency);  // sort k-mers based on frequency
        GraphKmerHashHapStrMapIterVecTmp.resize(128); // keep only the first 128 k-mers
    }

    std::lock_guard<std::mutex> mtx_locker(mtxCI);
    vector<uint64_t>().swap(kmerHashVec); // free memory
    GraphKmerHashHapStrMapIterVec = move(GraphKmerHashHapStrMapIterVecTmp);  // save k-mer informations

    return 0;
}


/**
 * @author zezhen du
 * @date 2023/06/11
 * @version v1.0
 * @brief split the GT information
 * 
 * @param gtTxt     GT information
 * 
 * @return vector<string> gtVecTmp
**/
vector<string> construct_index::gt_split(const string & gtTxt)
{
    vector<string> gtVecTmp;

    if (gtTxt == ".") {
        return gtVecTmp;
    }
    
    if (gtTxt.find("/") != string::npos) {
        gtVecTmp = split(gtTxt, "/");
    } else if (gtTxt.find("|") != string::npos) {
        gtVecTmp = split(gtTxt, "|");
    } else {
        try {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: sample has only one genotype, attempting to correct to diploid -> " << gtTxt << endl;
            stoul(gtTxt);
            gtVecTmp.push_back(gtTxt);
        } catch (const std::invalid_argument&) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: GT is not separated by '/' or '|' -> " << gtTxt << endl;
            exit(1);
        } catch (const std::out_of_range&) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: GT is not separated by '/' or '|' -> " << gtTxt << endl;
            exit(1);
        }
    }

    return gtVecTmp;
}


/**
 * @author zezhen du
 * @date 2023/06/29
 * @version v1.0.1
 * @brief Sets a specific bit of a given number to 1.
 * 
 * @param num         The value that needs to be modified.
 * @param bitIndex    The index of the bit to be set to 1.
 * 
 * @return void
**/
template <typename T>
void construct_index::set_bit_to_one(T& num, int bitIndex) {
    T mask = static_cast<T>(1) << bitIndex;
    num |= mask;
}

/**
 * @author zezhen du
 * @date 2023/06/29
 * @version v1.0.1
 * @brief Sets a specific bit of a given number to 0.
 * 
 * @param num         The value that needs to be modified.
 * @param bitIndex    The index of the bit to be set to 0.
 * 
 * @return void
**/
template <typename T>
void construct_index::set_bit_to_zero(T& num, int bitIndex) {
    T mask = ~(static_cast<T>(1) << bitIndex);
    num &= mask;
}

/**
 * @author zezhen du
 * @date 2023/06/29
 * @version v1.0.1
 * @brief Retrieves the value of a specific bit in a given number.
 * 
 * @param num         The value from which to retrieve the bit.
 * @param bitIndex    The index of the bit to be queried.
 * 
 * @return int
**/
template <typename T>
int construct_index::get_bit(T num, int bitIndex) {
    T mask = static_cast<T>(1) << bitIndex;
    return ((num & mask) ? 1 : 0);
}

/**
 * @author zezhen du
 * @date 2023/07/03
 * @version v1.0.1
 * @brief print bitmap.
 * 
 * @param num         The value from which to print.
 * 
 * @return int
**/
string construct_index::print_bits(int8_t num) {
    stringstream ss;
    for (int i = 7; i >= 0; i--) {
        ss << ((num >> i) & 1);
    }
    return ss.str();
}


// Explicitly instantiate the template
// int8_t
template void construct_index::set_bit_to_one<int8_t>(int8_t&, int);
template void construct_index::set_bit_to_zero<int8_t>(int8_t&, int);
template int construct_index::get_bit<int8_t>(int8_t, int);

// int32_t
template void construct_index::set_bit_to_one<int32_t>(int32_t&, int);
template void construct_index::set_bit_to_zero<int32_t>(int32_t&, int);
template int construct_index::get_bit<int32_t>(int32_t, int);

// int64_t
template void construct_index::set_bit_to_one<int64_t>(int64_t&, int);
template void construct_index::set_bit_to_zero<int64_t>(int64_t&, int);
template int construct_index::get_bit<int64_t>(int64_t, int);