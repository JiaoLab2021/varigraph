// g++ -c fastq_kmer.cpp -o fastq_kmer.o -std=c++17 -O3 -march=native

#include "../include/fastq_kmer.hpp"

using namespace std;


// kseq.h
KSEQ_INIT(gzFile, gzread)


/**
 * @author zezhen du
 * @date 2023/07/17
 * @version v1.0
 * @brief building the kmer index of sequencing read
 * 
 * @param GraphKmerCovFreMap     Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFreBitVecP>
 * @param fastqFileNameVec       the vector of sequencing read
 * @param kmerLen                the length of k-mer
 * @param threads                threads
 * 
 * @return void
**/
FastqKmer::FastqKmer(
    unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap, 
    const vector<string>& fastqFileNameVec, 
    const uint32_t& kmerLen, 
    const uint32_t & threads
) : GraphKmerHashHapStrMap_(GraphKmerHashHapStrMap), fastqFileNameVec_(fastqFileNameVec), kmerLen_(kmerLen), threads_(threads) {}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief building the kmer index of sequencing read
 * 
 * @return void
**/
void FastqKmer::build_fastq_index() {
    if (fastqFileNameVec_.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -f\n";
        exit(1);
    }

    for (auto fastqFileName : fastqFileNameVec_) {
        FastqKmer::fastq_file_open(fastqFileName);
    }

    malloc_trim(0); // 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2023/09/06
 * @version v1.0
 * @brief building the kmer index of sequencing read
 * 
 * @param fastqFileName       sequencing read
 * 
 * @return void
**/
void FastqKmer::fastq_file_open(
    const string & fastqFileName
) {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Collecting kmers from read on CPU: " << fastqFileName << endl;  // print log

    // Save the result of multithreading
    vector<future<vector<uint64_t> > > hashVecVec;  // vector<kmerVec>

    // open read file
    gzFile gzfp = gzopen(fastqFileName.c_str(), "rb");

    // Thread pool
    ThreadPool pool(threads_);
    const int MAX_THREADS_NUM = threads_*100;  // maximum number of queues

    // Initializes the thread pool
    pool.init();

    // Save a temporary list of sequences for multithreaded submission tasks
    vector<string> sequenceVecTmp;

    // open file
    if(!gzfp) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << fastqFileName 
                << "': No such file or directory." 
                << endl;
        exit(1);
    } else {
        kseq_t *ks;
        ks = kseq_init(gzfp);
    
        while( kseq_read(ks) >= 0 ) {
            // ks->name.s records the name
            // ks->seq.s records the sequence
            // Uppercase
            string sequence = ks->seq.s;
            transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

            // read size
            mReadBase += ks->seq.l;

            // Store sequences for multithreaded submission
            sequenceVecTmp.push_back(sequence);

            // If the value is greater than the threshold, submit the task
            if (sequenceVecTmp.size() >= MAX_THREADS_NUM) {
                // submit
                hashVecVec.push_back(
                    pool.submit(
                        fastq_kmer::fastq_file_open_run, 
                        sequenceVecTmp, 
                        kmerLen_, 
                        ref(GraphKmerHashHapStrMap_)
                    )
                );
                // clear
                sequenceVecTmp.clear();
            }

            // If the maximum number of queues is reached
            if (hashVecVec.size() >= MAX_THREADS_NUM) {
                // Save multithreaded results into hash table
                for (size_t i = 0; i < hashVecVec.size(); i++) {  // vector<kmerVec>
                    auto hashVecTmp = move(hashVecVec[i].get());

                    // Record coverage of variant k-mers
                    for (const auto& hash : hashVecTmp) {  // vector<uint64_t>
                        auto& kmerCov = GraphKmerHashHapStrMap_[hash].c;
                        // if (kmerCov < UINT8_MAX - 1) {  // Prevent variable out of bounds
                        if (kmerCov < UINT8_MAX) {  // Prevent variable out of bounds
                            kmerCov++;
                        }
                    }
                }
                hashVecVec.clear();
            }
        }

        // free memory, close file
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    // submit
    if (sequenceVecTmp.size() > 0)
    {
        hashVecVec.push_back(
            pool.submit(
                fastq_kmer::fastq_file_open_run, 
                sequenceVecTmp, 
                kmerLen_, 
                ref(GraphKmerHashHapStrMap_)
            )
        );
        // free memory
        vector<string>().swap(sequenceVecTmp);
    }
    
    // last save
    if (hashVecVec.size() > 0){
        // Save multithreaded results into hash table
        for (size_t i = 0; i < hashVecVec.size(); i++) {  // vector<kmerVec>
            auto hashVecTmp = move(hashVecVec[i].get());

            // Record coverage of variant k-mers
            for (const auto& hash : hashVecTmp) {  // vector<uint64_t>
                auto& kmerCov = GraphKmerHashHapStrMap_[hash].c;
                // if (kmerCov < UINT8_MAX - 1) {  // Prevent variable out of bounds
                if (kmerCov < UINT8_MAX) {  // Prevent variable out of bounds
                    kmerCov++;
                }
            }
        }
        vector<future<vector<uint64_t> > >().swap(hashVecVec);
    }

    // Close the thread pool
    pool.shutdown();

    // free memory
    malloc_trim(0);	 // 0 is for heap memory
}


/**
 * @author zezhen du
 * @date 2023/07/18
 * @version v1.0.1
 * @brief Save the k-mer statistical results of fastq
 * 
 * @param outputFileName
 * 
 * @return void
**/
void FastqKmer::save_index(const string & outputFileName) {
    // log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Reads index saved to file: " << outputFileName << endl;

    std::ofstream outputFile(outputFileName, std::ios::binary);
    if (outputFile) {
        // Read base
        outputFile.write(reinterpret_cast<const char*>(&mReadBase), sizeof(uint64_t));

        for (const auto& kvp : GraphKmerHashHapStrMap_) {
            const uint64_t kmerHash = kvp.first;
            const kmerCovFreBitVec& kmerCovFreBitVecStr = kvp.second;
            const vector<int8_t>& BitVec = kmerCovFreBitVecStr.BitVec;
            uint64_t BitVecLen = BitVec.size();

            // Write kmerHash
            outputFile.write(reinterpret_cast<const char*>(&kmerHash), sizeof(uint64_t));

            // Write kmerCov and kmerFre
            outputFile.write(reinterpret_cast<const char*>(&kmerCovFreBitVecStr.c), sizeof(uint8_t));
            outputFile.write(reinterpret_cast<const char*>(&kmerCovFreBitVecStr.f), sizeof(uint8_t));

            // Write BitVec
            outputFile.write(reinterpret_cast<const char*>(&BitVecLen), sizeof(uint64_t));
            for (const int8_t& Bit : BitVec) {
                outputFile.write(reinterpret_cast<const char*>(&Bit), sizeof(int8_t));
            }
        }

        outputFile.close();
    } else {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << outputFileName << "': No such file or directory." << endl;
        exit(1);
    }
}


/**
 * @author zezhen du
 * @date 2023/07/18
 * @version v1.0.1
 * @brief load the k-mer index from disk
 * 
 * @param inputFileName
 * 
 * @return void
**/
void FastqKmer::load_index(
    const string & inputFileName
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Reads index loaded from file: " << inputFileName << endl;  // print log

    std::ifstream inputFile(inputFileName, std::ios::binary);
    if (inputFile) {
        // Clear the existing data in GraphKmerHashHapStrMap_
        GraphKmerHashHapStrMap_.clear();

        uint64_t kmerHash;
        kmerCovFreBitVec kmerCovFreBitVecStr;
        uint8_t c;
        uint8_t f;
        uint64_t BitVecLen;

        // Read base
        inputFile.read(reinterpret_cast<char*>(&mReadBase), sizeof(uint64_t));

        // Read kmerHash, kmerCov, and kmerFre from the file
        while (inputFile.read(reinterpret_cast<char*>(&kmerHash), sizeof(uint64_t))) {

            // Read kmerCov and kmerFre
            inputFile.read(reinterpret_cast<char*>(&c), sizeof(uint8_t));
            inputFile.read(reinterpret_cast<char*>(&f), sizeof(uint8_t));

            // Read BitVec length
            inputFile.read(reinterpret_cast<char*>(&BitVecLen), sizeof(uint64_t));

            // Read BitVec
            vector<int8_t> BitVec(BitVecLen);
            for (auto& Bit : BitVec) {
                inputFile.read(reinterpret_cast<char*>(&Bit), sizeof(int8_t));
            }

            // Store data in the map
            kmerCovFreBitVecStr.c = c;
            kmerCovFreBitVecStr.f = f;
            kmerCovFreBitVecStr.BitVec = BitVec;
            GraphKmerHashHapStrMap_[kmerHash] = kmerCovFreBitVecStr;
        }

        inputFile.close();
    } else {
        // log
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << inputFileName << "': No such file or directory." << endl;
        exit(1);
    }
}



/**
 * @author zezhen du
 * @date 2023/08/01
 * @version v1.0.1
 * @brief Component index multithreaded functions to reduce thread locks
 * 
 * @param sequenceVec                sequence vector
 * @param kmerLen                    The length of kmer
 * @param GraphKmerHashHapStrMap     Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFreBitVec>
 * 
 * @return hashVec               vector<uint64_t>
**/
vector<uint64_t> fastq_kmer::fastq_file_open_run(
    vector<string> sequenceVec, 
    uint32_t kmerLen, 
    const unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap
) {
    // A temporary list of hash values
    vector<uint64_t> hashVec;

    // Iterate through the list of sequences
    for (const auto& seq : sequenceVec) {
        // save the result
        vector<uint64_t> hashVecTmp = kmerBit::kmer_sketch_fastq(seq, kmerLen, ref(GraphKmerHashHapStrMap));

        // Save to the total vector
        hashVec.insert(hashVec.end(), hashVecTmp.begin(), hashVecTmp.end());
    }

    return hashVec;
}