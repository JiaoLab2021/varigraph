// g++ -c fastq_kmer.cpp -o fastq_kmer.o -std=c++17 -O3 -march=native

#include "../include/fastq_kmer.hpp"

using namespace std;


// kseq.h 打开文件
KSEQ_INIT(gzFile, gzread)


/**
 * @author zezhen du
 * @date 2023/07/17
 * @version v1.0
 * @brief building the kmer index of sequencing read
 * 
 * @param GraphKmerCovFreMap     Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFre>
 * @param fastqFileNameVec       the vector of sequencing read
 * @param kmerLen                the length of k-mer
 * @param threads                threads
 * 
 * @return void
**/
FastqKmer::FastqKmer(
    unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap, 
    const vector<string>& fastqFileNameVec, 
    const uint32_t& kmerLen, 
    const uint32_t & threads
) : GraphKmerCovFreMap_(GraphKmerCovFreMap), fastqFileNameVec_(fastqFileNameVec), kmerLen_(kmerLen), threads_(threads) {}


/**
 * @author zezhen du
 * @date 2023/06/27
 * @version v1.0.1
 * @brief building the kmer index of sequencing read
 * 
 * @return void
**/
void FastqKmer::build_fastq_index()
{
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
 * @date 2023/06/27
 * @version v1.0
 * @brief building the kmer index of sequencing read
 * 
 * @param fastqFileName       sequencing read
 * 
 * @return void
**/
void FastqKmer::fastq_file_open(
    const string & fastqFileName
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Collecting kmers from reads: " << fastqFileName << endl;  // print log

    // 保存多线程的结果
    vector<future<vector<uint64_t> > > hashVecVec;  // vector<kmerVec>

    // open base file
    gzFile gzfp = gzopen(fastqFileName.c_str(), "rb");

    // 进程池
    ThreadPool pool(threads_);
    const int MAX_THREADS_NUM = threads_*100;  // maximum number of queues

    // 初始化线程池
    pool.init();

    // 保存序列的临时列表，用于多线程提交任务
    vector<string> sequenceVecTmp;

    // 打开文件
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
            // ks->name.s 记录的是名字
            // ks->seq.s 记录的是序列
            // 转大写
            string sequence = ks->seq.s;
            transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

            // 存储序列，用于多线程提交
            sequenceVecTmp.push_back(sequence);

            // 如果大于阈值则提交任务
            if (sequenceVecTmp.size() >= MAX_THREADS_NUM) {
                // submit
                hashVecVec.push_back(
                    pool.submit(
                        fastq_kmer::fastq_file_open_run, 
                        sequenceVecTmp, 
                        kmerLen_, 
                        ref(GraphKmerCovFreMap_)
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
                        if (GraphKmerCovFreMap_[hash].c < UINT8_MAX - 1) {  // Prevent variable out of bounds
                            ++GraphKmerCovFreMap_[hash].c;
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
                ref(GraphKmerCovFreMap_)
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
                if (GraphKmerCovFreMap_[hash].c < UINT8_MAX - 1) {  // Prevent variable out of bounds
                    ++GraphKmerCovFreMap_[hash].c;
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
 * @brief Save the kmer statistical results of fastq
 * 
 * @param outputFileName
 * 
 * @return void
**/
void FastqKmer::save_index(
    const string & outputFileName
)
{
    std::ofstream outputFile(outputFileName, std::ios::binary);
    if (outputFile) {
        for (const auto& kvp : GraphKmerCovFreMap_) {
            const uint64_t kmerHash = kvp.first;
            const kmerCovFre& covFre = kvp.second;

            // Write kmerHash
            outputFile.write(reinterpret_cast<const char*>(&kmerHash), sizeof(uint64_t));

            // Write kmerCov and kmerFre
            outputFile.write(reinterpret_cast<const char*>(&covFre.c), sizeof(uint8_t));
            outputFile.write(reinterpret_cast<const char*>(&covFre.f), sizeof(uint8_t));
        }

        outputFile.close();
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "Read index saved to file: " << outputFileName << endl;
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
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Read index loaded from file: " << inputFileName << endl;  // print log

    std::ifstream inputFile(inputFileName, std::ios::binary);
    if (inputFile) {
        // Clear the existing data in GraphKmerCovFreMap_
        GraphKmerCovFreMap_.clear();

        uint64_t kmerHash;
        kmerCovFre covFre;

        // Read kmerHash, kmerCov, and kmerFre from the file
        while (inputFile.read(reinterpret_cast<char*>(&kmerHash), sizeof(uint64_t))) {
            inputFile.read(reinterpret_cast<char*>(&covFre.c), sizeof(uint8_t));
            inputFile.read(reinterpret_cast<char*>(&covFre.f), sizeof(uint8_t));

            // Store the values in GraphKmerCovFreMap_
            GraphKmerCovFreMap_[kmerHash] = covFre;
        }

        inputFile.close();
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "Read index loaded from file: " << inputFileName << endl;
    } else {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << inputFileName << "': No such file or directory." << endl;
        exit(1);
    }
}



/**
 * @author zezhen du
 * @date 2023/07/17
 * @version v1.0.1
 * @brief 构件索引多线程函数，减少线程锁
 * 
 * @param sequenceVec            sequence vector
 * @param kmerLen                The length of kmer
 * @param GraphKmerCovFreMap     Record the coverage and frequency of all k-mers in the graph: map<kmerHash, kmerCovFre>
 * 
 * @return hashVec               vector<uint64_t>
**/
vector<uint64_t> fastq_kmer::fastq_file_open_run(
    vector<string> sequenceVec, 
    uint32_t kmerLen, 
    const unordered_map<uint64_t, kmerCovFre> & GraphKmerCovFreMap
)
{
    // 存储kmer的hash值
    vector<uint64_t> hashVec;

    // 遍历序列提交任务
    for (const auto& seq : sequenceVec)
    {
        // save the result
        vector<uint64_t> hashVecTmp = kmerBit::kmer_sketch_fastq(seq, kmerLen, ref(GraphKmerCovFreMap));

        // Save to the total vector
        hashVec.insert(hashVec.end(), hashVecTmp.begin(), hashVecTmp.end());
    }

    return hashVec;
}