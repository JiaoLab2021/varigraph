#include "../include/fastq_kmer.cuh"

using namespace std;

__device__ char toupper_device(char c) {
    return (c >= 'a' && c <= 'z') ? (c - 'a' + 'A') : c;
}

// kseq.h
KSEQ_INIT(gzFile, gzread)

/**
 * @author zezhen du
 * @date 2024/04/30
 * @version v1.0.1
 * @brief building the kmer index of sequencing read
 * 
 * @return void
**/
void FastqKmerKernel::build_fastq_index_kernel() {
    if (fastqFileNameVec_.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -f\n";
        exit(1);
    }

    for (auto fastqFileName : fastqFileNameVec_) {
        fastq_file_open_kernel(fastqFileName);
    }

    malloc_trim(0); // 0 is for heap memory
}

/**
 * @author zezhen du
 * @date 2024/04/30
 * @version v1.0
 * @brief building the kmer index of sequencing read
 * 
 * @param fastqFileName       sequencing read
 * 
 * @return void
**/
void FastqKmerKernel::fastq_file_open_kernel(const string & fastqFileName) {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Collecting kmers from read on GPU: " << fastqFileName << endl;  // print log

    // *************************************************** GPU parameters *************************************************** //
    const uint32_t MAX_CHUNK_SIZE = buffer_ * 1024 * 1024;  // Maximum chunk size
    const uint32_t MAX_KMER_NUM = MAX_CHUNK_SIZE - kmerLen_ + 1;  // sequence.size() - kmerLen_ + 1
    uint32_t bufferSize = 0;  // Buffer size

    // Calculate the size of the grid and thread blocks
    int blockSize = 64;
    uint32_t numBlocks = (MAX_KMER_NUM / blockSize) + 1;

    // sequence
    char* strH = new char[MAX_CHUNK_SIZE];  // host memory
    char* strD;  // device memory
    gpuErrchk(cudaMallocManaged((void**)&strD, sizeof(char) * MAX_CHUNK_SIZE));  // Allocate unified memory to store the sequence

    // k-mers
    uint64_t * kmersD;  // device memory
    gpuErrchk(cudaMallocManaged(&kmersD, sizeof(uint64_t) * MAX_KMER_NUM));  // Allocate unified memory to store the hash results
    thrust::device_ptr<uint64_t> kmersD_ptr(kmersD);  // Create device pointers

    // Allocate space to store the unique keys and their counts on the device
    thrust::device_vector<uint64_t> uniqueKmersD(MAX_KMER_NUM);
    thrust::device_vector<uint8_t> uniqueKmersCountsD(MAX_KMER_NUM);

    // Create host vectors
    vector<uint64_t> uniqueKmersH(MAX_KMER_NUM);
    vector<uint8_t> uniqueKmersCountsH(MAX_KMER_NUM);

    // *************************************************** Thread pool *************************************************** //
    ThreadPool pool(threads_);
    pool.init();  // Initializes the thread pool

    vector<future<void> > futureVec;  // Save the result of multithreading

    // *************************************************** Open file *************************************************** //
    gzFile gzfp = gzopen(fastqFileName.c_str(), "rb");

    // open file
    if(!gzfp) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << fastqFileName << "': No such file or directory." << endl;
        exit(1);
    }

    kseq_t *ks;
    ks = kseq_init(gzfp);

    while( kseq_read(ks) >= 0 ) {
        // ks->name.s records the name
        // ks->seq.s records the sequence
        uint32_t seqLen = ks->seq.l;

        if (bufferSize + seqLen > MAX_CHUNK_SIZE) {
            // Copy the sequence to the end of the existing data in strD
            gpuErrchk(cudaMemcpy(strD, strH, sizeof(char) * bufferSize, cudaMemcpyHostToDevice));

            // Convert to upper case on the GPU
            thrust::device_ptr<char> strD_ptr(strD);
            thrust::transform(strD_ptr, strD_ptr + MAX_CHUNK_SIZE, strD_ptr, [] __device__ (char c) {
                return toupper_device(c);
            });

            // ************************************************************* k-mer sketch ************************************************************* //
            // Launch the GPU kernel
            gpuErrchk(cudaMemset(kmersD, -1, sizeof(uint64_t) * MAX_KMER_NUM));  // Initialize the memory to UINT64_MAX
            kmer_sketch_kernel<<<numBlocks, blockSize>>>(strD, bufferSize, kmerLen_, kmersD);
            gpuErrchk(cudaDeviceSynchronize());  // Wait for the GPU to finish operations

            // Fill the vectors with initial values
            thrust::fill(uniqueKmersD.begin(), uniqueKmersD.end(), 0);
            thrust::fill(uniqueKmersCountsD.begin(), uniqueKmersCountsD.end(), 0);

            // Perform parallel quicksort on GPU
            thrust::device_ptr<uint64_t> kmersD_end = kmersD_ptr + bufferSize;
            thrust::sort(kmersD_ptr, kmersD_end);

            // Use reduce_by_key to compress and calculate counts
            auto end = thrust::reduce_by_key(kmersD_ptr, kmersD_end, thrust::make_constant_iterator(1),
                                             uniqueKmersD.begin(), uniqueKmersCountsD.begin(), thrust::equal_to<uint64_t>(), capped_add());
            gpuErrchk(cudaDeviceSynchronize());  // Wait for the GPU to finish operations

            // Calculate the number of unique hash
            uint32_t numUniqueKmers = end.first - uniqueKmersD.begin();

            // Copy k-mer coverage from device to host
            thrust::copy(uniqueKmersD.begin(), uniqueKmersD.end(), uniqueKmersH.begin());
            thrust::copy(uniqueKmersCountsD.begin(), uniqueKmersCountsD.end(), uniqueKmersCountsH.begin());

            // Calculate the size of each chunk
            uint32_t chunkSize = numUniqueKmers / threads_;

            // Add k-mer coverage to the map
            for (uint32_t i = 0; i < threads_; i++) {
                uint32_t startIndex = i * chunkSize;
                uint32_t endIndex = (i == threads_ - 1) ? numUniqueKmers : startIndex + chunkSize;

                // Break the loop if startIndex or endIndex is out of bounds
                if (startIndex >= numUniqueKmers) {
                    break;
                } else if (endIndex > numUniqueKmers) {
                    endIndex = numUniqueKmers;
                }

                futureVec.push_back(
                    pool.submit(
                        fastq_kmer_kernel::add_cov_to_map, 
                        ref(uniqueKmersH), 
                        ref(uniqueKmersCountsH), 
                        ref(GraphKmerHashHapStrMap_),
                        startIndex,
                        endIndex
                    )
                );
            }
            // Wait
            for (auto& future : futureVec) {
                future.get();
            }
            // Clear
            vector<future<void> >().swap(futureVec);

            // Reset the buffer size
            bufferSize = 0;
        }

        // Copy the sequence to the end of the existing data in strH
        memcpy(strH + bufferSize, ks->seq.s, sizeof(char) * seqLen);
        bufferSize += seqLen;  // Update the buffer size

        // Add the 'N' character to the end of the sequence in strH
        strH[bufferSize] = 'N';
        bufferSize += 1;  // Update the buffer sizedelete[] strH;

        // read size
        mReadBase += ks->seq.l;
    }

    // free memory, close file
    kseq_destroy(ks);
    gzclose(gzfp);

    // Last batch
    if (bufferSize > 0) {
        // Copy the sequence to the end of the existing data in strD
        gpuErrchk(cudaMemcpy(strD, strH, sizeof(char) * bufferSize, cudaMemcpyHostToDevice));

        // Convert to upper case on the GPU
        thrust::device_ptr<char> strD_ptr(strD);
        thrust::transform(strD_ptr, strD_ptr + MAX_CHUNK_SIZE, strD_ptr, [] __device__ (char c) {
            return toupper_device(c);
        });

        // ************************************************************* k-mer sketch ************************************************************* //
        // Launch the GPU kernel
        gpuErrchk(cudaMemset(kmersD, -1, sizeof(uint64_t) * MAX_KMER_NUM));  // Initialize the memory to UINT64_MAX
        kmer_sketch_kernel<<<numBlocks, blockSize>>>(strD, bufferSize, kmerLen_, kmersD);
        gpuErrchk(cudaDeviceSynchronize());  // Wait for the GPU to finish operations

        // Fill the vectors with initial values
        thrust::fill(uniqueKmersD.begin(), uniqueKmersD.end(), 0);
        thrust::fill(uniqueKmersCountsD.begin(), uniqueKmersCountsD.end(), 0);

        // Perform parallel quicksort on GPU
        thrust::device_ptr<uint64_t> kmersD_end = kmersD_ptr + bufferSize;
        thrust::sort(kmersD_ptr, kmersD_end);

        // Use reduce_by_key to compress and calculate counts
        auto end = thrust::reduce_by_key(kmersD_ptr, kmersD_end, thrust::make_constant_iterator(1),
                                         uniqueKmersD.begin(), uniqueKmersCountsD.begin(), thrust::equal_to<uint64_t>(), capped_add());
        gpuErrchk(cudaDeviceSynchronize());  // Wait for the GPU to finish operations

        // Calculate the number of unique hash
        uint32_t numUniqueKmers = end.first - uniqueKmersD.begin();

        // Copy k-mer coverage from device to host
        thrust::copy(uniqueKmersD.begin(), uniqueKmersD.end(), uniqueKmersH.begin());
        thrust::copy(uniqueKmersCountsD.begin(), uniqueKmersCountsD.end(), uniqueKmersCountsH.begin());

        // Calculate the size of each chunk
        uint32_t chunkSize = numUniqueKmers / threads_;

        // Add k-mer coverage to the map
        for (uint32_t i = 0; i < threads_; i++) {
            uint32_t startIndex = i * chunkSize;
            uint32_t endIndex = (i == threads_ - 1) ? numUniqueKmers : startIndex + chunkSize;

            // Break the loop if startIndex or endIndex is out of bounds
            if (startIndex >= numUniqueKmers) {
                break;
            } else if (endIndex > numUniqueKmers) {
                endIndex = numUniqueKmers;
            }

            // submit
            futureVec.push_back(
                pool.submit(
                    fastq_kmer_kernel::add_cov_to_map, 
                    ref(uniqueKmersH), 
                    ref(uniqueKmersCountsH), 
                    ref(GraphKmerHashHapStrMap_),
                    startIndex,
                    endIndex
                )
            );
        }

        // Wait
        for (auto& future : futureVec) {
            future.get();
        }
        // Clear
        vector<future<void> >().swap(futureVec);

        // Reset the buffer size
        bufferSize = 0;
    }

    // Close the thread pool
    pool.shutdown();

    // Free
    delete[] strH;

    // Free unified memory
    gpuErrchk(cudaFree(strD));
    gpuErrchk(cudaFree(kmersD));

    // free memory
    malloc_trim(0);	 // 0 is for heap memory
}


/**
 * This function adds coverage to a map of k-mers.
 * 
 * @param kmerHashVec The hash vector of the k-mer.
 * @param kmerCovVec The coverage vector of the k-mer.
 * @param GraphKmerHashHapStrMap The map of k-mers to their coverage and frequency bit vectors.
 * @param startIndex The starting index of the k-mer.
 * @param endIndex The ending index of the k-mer.
*/
void fastq_kmer_kernel::add_cov_to_map(
    const vector<uint64_t>& kmerHashVec, 
    const vector<uint8_t>& kmerCovVec,
    unordered_map<uint64_t, kmerCovFreBitVec>& GraphKmerHashHapStrMap,
    uint32_t startIndex,
    uint32_t endIndex
) {
    for (uint32_t i = startIndex; i < endIndex; i++) {
        uint64_t kmerHash = kmerHashVec[i];
        uint8_t kmerCov = kmerCovVec[i];

        auto it = GraphKmerHashHapStrMap.find(kmerHash);
        if (it == GraphKmerHashHapStrMap.end()) {
            continue;
        }

        auto& c = it->second.c;
        if (static_cast<uint32_t>(c) + static_cast<uint32_t>(kmerCov) <= UINT8_MAX) {  // Prevent variable out of bounds
            c += kmerCov;
        } else {
            c = UINT8_MAX;
        }
    }
}