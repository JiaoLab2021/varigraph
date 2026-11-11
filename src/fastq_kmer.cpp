// g++ -c fastq_kmer.cpp -o fastq_kmer.o -std=c++17 -O3 -march=native

#include "../include/fastq_kmer.hpp"
#include "../include/logger.hpp"
#include "../include/GzChunkReader.hpp"
#include "../include/progress_tracker.hpp"
#include <chrono>
#include <deque>
#include <mutex>

using namespace std;


// kseq.h
KSEQ_INIT(gzFile, read_gzip_checked)


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
        error_stream() << "Parameter error: -f\n";
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
void FastqKmer::fastq_file_open(const string& fastqFileName) {
    stage_stream() << "Collecting read k-mers on CPU: " << fastqFileName << '\n';

    gzFile gzfp = gzopen(fastqFileName.c_str(), "rb");
    std::unique_ptr<gzFile_s, decltype(&gzclose)> input(gzfp, &gzclose);
    if (!gzfp) {
        error_stream() << "'" << fastqFileName << "': No such file or directory.\n";
        exit(1);
    }

    deque<future<void>> results;
    auto progress = ProgressTracker::Every(1000000);
    ThreadPool pool(threads_);
    const size_t batch_size = static_cast<size_t>(threads_) * 100;
    vector<string> sequences;
    auto submit_batch = [&] {
        results.push_back(pool.submit([this, &progress, batch = std::move(sequences)]() mutable {
            const auto read_count = batch.size();
            auto hashes = fastq_kmer::fastq_file_open_run(std::move(batch), kmerLen_, GraphKmerHashHapStrMap_);
            // Sketching only reads keys; serialize coverage writes without changing map structure.
            {
                std::lock_guard<std::mutex> lock(mtx);
                for (const auto hash : hashes) {
                    auto& coverage = GraphKmerHashHapStrMap_.at(hash).c;
                    if (coverage < UINT8_MAX) ++coverage;
                }
            }
            progress.hit(read_count);
        }));
        sequences.clear();
        // The bounded pool supplies backpressure; a slow first batch must not stop input.
        for (auto result = results.begin(); result != results.end();) {
            if (result->wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
                result->get();
                result = results.erase(result);
            } else {
                ++result;
            }
        }
    };

    kseq_t* ks = kseq_init(gzfp);
    std::unique_ptr<kseq_t, decltype(&kseq_destroy)> reader(ks, &kseq_destroy);
    int read_status;
    while ((read_status = kseq_read(ks)) >= 0) {
        string sequence = ks->seq.s;
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        mReadBase += ks->seq.l;
        sequences.push_back(std::move(sequence));
        if (sequences.size() >= batch_size) submit_batch();
    }
    if (read_status < -1) throw std::runtime_error("Invalid sequence or quality record in '" + fastqFileName + "'.");
    reader.reset();
    input.reset();

    if (!sequences.empty()) submit_batch();
    for (auto& result : results) result.get();
    results.clear();
    pool.stop();
    progress.finish();
    malloc_trim(0);
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
    stage_stream() << "Writing read index: " << outputFileName << endl;

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
        error_stream()
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
    stage_stream() << "Reading read index: " << inputFileName << endl;  // print log

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
        error_stream()
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
