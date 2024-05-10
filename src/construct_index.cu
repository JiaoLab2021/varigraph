// g++ -c construct_index.cpp -std=c++17 -lz -O3 -march=native

#include "../include/kmer.cuh"
#include "../include/construct_index.cuh"

using namespace std;

// global variable
bool debugConstructKernel = false;

// Constructor
ConstructIndexKernel::ConstructIndexKernel(
    const string& refFileName, 
    const string& vcfFileName, 
    const string& inputGraphFileName, 
    const string& outputGraphFileName, 
    const bool& fastMode, 
    const uint32_t& kmerLen, 
    const uint32_t& vcfPloidy, 
    const bool& debug, 
    const uint32_t& threads, 
    const int buffer
) : ConstructIndex(refFileName, vcfFileName, inputGraphFileName, outputGraphFileName, fastMode, kmerLen, vcfPloidy, debug, threads) {
    buffer_ = buffer;  // GPU buffer size
    // debug
    debugConstructKernel = debug;
}


/**
 * @author zezhen du
 * @date 2024/04/29
 * @version v1.0.1
 * @brief Making Counting Bloom Filter
 * 
 * @return void
**/
void ConstructIndexKernel::make_mbf_kernel() {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Initiating computation of k-mer frequencies in the reference genome on GPU ...\n";

    /* *************************************************** making or load *************************************************** */
    uint64_t bfSize = mGenomeSize - mKmerLen + 1;
    double errorRate = 0.01;
    mbfD = new BloomFilterKernel(bfSize, errorRate);

    // making
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Making Counting Bloom Filter with a false positive rate of " << errorRate << " ...\n";

    const uint32_t MAX_CHUNK_SIZE = buffer_ * 1024 * 1024;  // GPU buffer size
    const uint32_t MAX_KMER_NUM = MAX_CHUNK_SIZE - mKmerLen + 1;  // sequence.size() 1 mKmerLen + 1

    // Calculate the size of the grid and thread blocks
    int blockSize = 64;
    uint32_t numBlocks = (MAX_KMER_NUM / blockSize) + 1;

    // sequence
    char* strD;  // device memory
    gpuErrchk(cudaMallocManaged((void**)&strD, sizeof(char) * MAX_CHUNK_SIZE));  // Allocate unified memory to store the sequence

    // k-mers
    uint64_t* kmersD;  // device memory
    gpuErrchk(cudaMallocManaged(&kmersD, sizeof(uint64_t) * MAX_KMER_NUM));  // Allocate unified memory to store the hash results

    for (const auto& [chromosome, sequence] : mFastaSeqMap) {  // map<chromosome, sequence>
        const char* sequencePtr = sequence.c_str();
        uint32_t sequenceSize = sequence.size();

        for (uint32_t i = 0; i < sequenceSize; i += MAX_CHUNK_SIZE) {
            uint32_t chunkSize = std::min(MAX_CHUNK_SIZE, sequenceSize - i);
            uint32_t chunkKmerNum = chunkSize - mKmerLen + 1;

            gpuErrchk(cudaMemcpy(strD, sequencePtr + i, sizeof(char) * chunkSize, cudaMemcpyHostToDevice));  // Copy the sequence to the device

            // ************************************************************* k-mer sketch ************************************************************* //
            // Launch the GPU kernel
            gpuErrchk(cudaMemset(kmersD, -1, sizeof(uint64_t) * MAX_KMER_NUM));  // Initialize the memory to UINT64_MAX
            kmer_sketch_kernel<<<numBlocks, blockSize>>>(strD, chunkSize, mKmerLen, kmersD);
            gpuErrchk(cudaDeviceSynchronize());  // Wait for the GPU to finish operations

            // ************************************************************* Counting Bloom Filter ************************************************************* //
            // Add the k-mers to the counting bloom filter
            mbfD->add_kernel(kmersD, chunkKmerNum);
        }

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Chromosome '" << chromosome << "' processed successfully ...\n";
    }
    
    // copies the content of the Counting Bloom Filter from the device (GPU) to the host (CPU).
    mbfD->copyFilterDToHost();
    gpuErrchk(cudaDeviceSynchronize());  // Wait for the GPU to finish operations

    // Free unified memory
    gpuErrchk(cudaFree(strD));
    gpuErrchk(cudaFree(kmersD));

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Counting Bloom Filter constructed successfully ..." << endl << endl;

    cerr << "           - " << "Counting Bloom Filter size: " << mbfD->get_size() << endl;
    cerr << "           - " << "Hash functions count: " << mbfD->get_num() << endl;
    cerr << fixed << setprecision(2);
    cerr << "           - " << "Counting Bloom Filter usage rate: " << mbfD->get_cap() << endl << endl << endl;
    cerr << defaultfloat << setprecision(6);

    malloc_trim(0);	// 0 is for heap memory
}

/**
 * @author zezhen du
 * @date 2024/04/28
 * @version v1.0
 * @brief building the k-mer index of graph on GPU
 * 
 * @return void
**/
void ConstructIndexKernel::index_kernel() {
    // Log the start of the graph index construction on the GPU
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
                    index_run_kernel,
                    chromosome, 
                    iter, 
                    ref(startNodeMap), 
                    ref(fastMode_), 
                    ref(mKmerLen), 
                    mbfD, 
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

        if (KmerHapBitMap.empty()) {
            continue;
        }

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
 * @brief graph index for k-mer on GPU
 * 
 * @date 2024/04/29
 * 
 * @param chromosome            mGraphMap output by construct��map<chr, map<start, nodeSrt> >
 * @param nodeIter              node iterator
 * @param startNodeMap          Chromosome all nodes
 * @param fastMode              fast mode
 * @param kmerLen               the length of kmer
 * @param bfD                   Counting Bloom Filter on GPU
 * @param vcfPloidy             ploidy of genotypes in VCF file
 * @param hapIdxQRmap           map<hapIdx, tuple<quotient, remainder> >
 * 
 * @return {nodeIter, tmpKmerHapBitMap, kmerHashFreMap}     kmer: map<kmerHash, vector<int8_t> >
**/
tuple<map<uint32_t, nodeSrt>::iterator, unordered_map<uint64_t, vector<int8_t> >, map<uint64_t, uint8_t> > index_run_kernel(
    string chromosome, 
    map<uint32_t, nodeSrt>::iterator nodeIter, 
    const map<uint32_t, nodeSrt>& startNodeMap, 
    const bool& fastMode, 
    const uint32_t& kmerLen, 
    BloomFilterKernel* bfD, 
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
        if (debugConstructKernel) {
            cerr << "Node Start:" << nodeIter->first << ", Haplotype:" << haplotype << ", GT:" << +gt << ", Upstream:" << upDownSeq.first << ", Current:" << seqTmp << ", Downstream:" << upDownSeq.second << endl << endl;
        }

        seqTmp = upDownSeq.first + seqTmp + upDownSeq.second;  // Add upstream and downstream sequences

        // k-mer indexing
        map<uint8_t, unordered_set<uint64_t> > freKmerHashSetMap;  // map<frequency, unordered_set<kmerHash> >
        kmer_sketch_construct_kernel(seqTmp, kmerLen, freKmerHashSetMap, bfD);

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
                if (gt != (uint16_t)0 && bfD->find(kmerHash) && construct_index::get_bit(emplacedValue->second.first[0], 0) == 0) {
                    // set the last bit of hapBit to 1
                    construct_index::set_bit_to_one(emplacedValue->second.first.back(), 7);
                }
            }

            // Update MIN_KMER_FRE if the frequency is greater than or equal to 2
            if (frequency >= 2) {
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