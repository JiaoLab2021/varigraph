#include "../include/counting_bloom_filter.cuh"


// This function performs an atomic addition operation on a uint8_t value, ensuring that the result does not exceed 255.
__device__ uint8_t atomicAddUint8(uint8_t* address, uint8_t val) {
    unsigned int* base_address = (unsigned int*) ((char*)address - ((size_t)address & 3));
    unsigned int long_val = ((unsigned int)val) << 8 * ((size_t)address & 3);
    unsigned int long_old, long_new_old;
    do {
        long_old = *base_address;
        if ((long_old >> 8 * ((size_t)address & 3) & 0xff) + val > 255) {
            long_val = ((unsigned int)(255 - (long_old >> 8 * ((size_t)address & 3) & 0xff))) << 8 * ((size_t)address & 3);
        }
        long_new_old = atomicCAS(base_address, long_old, long_old + long_val);
    } while (long_old != long_new_old);
    return (long_new_old >> 8 * ((size_t)address & 3)) & 0xff;
}

// This kernel function performs counting for a bloom filter. It takes an array of unique hashes, their counts, and a filter array.
// Each thread handles one unique hash. It calculates the position in the filter array by taking the modulo of the hash value with the size of the filter.
// Then it atomically adds the count of the hash to the value at the calculated position in the filter array.
__global__ void add_filter_kernel(uint64_t* uniqueHashsD, uint8_t* uniqueHashCountsD, uint32_t numUniqueHashs, uint8_t* filter, uint64_t size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < numUniqueHashs) {
        uint64_t position = uniqueHashsD[idx] % size;
        atomicAddUint8(&filter[position], uniqueHashCountsD[idx]);
    }
}

// This kernel function checks if a k-mer exists in the bloom filter. It takes an array of hash values, a filter array, and an output array.
// Each thread handles one k-mer. It calculates the position in the filter array by taking the modulo of the hash value with the size of the filter.
// Then it checks if the value at the calculated position in the filter array is zero. If it is, it sets the corresponding value in the output array to false.
// This is used to query the bloom filter on the GPU.
__global__ void check_filter_kernel(uint64_t* hashValues, uint32_t numKmers, uint8_t* filter, uint64_t size, bool* findBools) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < numKmers) {
        uint64_t hashValue = hashValues[idx];
        if (filter[hashValue % size] == 0) {
            findBools[idx] = false;
        }
    }
}

// This kernel function checks if a k-mer exists in the bloom filter and updates the counts array accordingly.
// It takes an array of hash values, a filter array, the size of the filter, and an output counts array.
// Each thread handles one k-mer. It calculates the position in the filter array by taking the modulo of the hash value with the size of the filter.
// Then it checks if the value at the calculated position in the filter array is less than the current count for the k-mer.
// If it is, it updates the count for the k-mer to the value in the filter array.
// This is used to query the bloom filter on the GPU and update the counts of k-mers.
__global__ void count_filter_kernel(uint64_t* hashValues, uint32_t numKmers, uint8_t* filter, uint64_t size, uint8_t* counts) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < numKmers) {
        uint64_t hashValue = hashValues[idx];
        uint64_t position = hashValue % size;
        if (filter[position] < counts[idx]) {
            counts[idx] = filter[position];
        }
    }
}

/*
 * Variable: kmersD
 * -----------------
 * This variable is declared on the device (GPU) and is only accessible 
 * from the device code (i.e., CUDA kernels). It cannot be directly accessed 
 * from the host (CPU) code.
 */
void BloomFilterKernel::add_kernel(uint64_t * kmersD, uint32_t numKmers) {
    // Allocate unified memory to store the hash results
    thrust::device_vector<uint64_t> hashValues(numKmers);

    // Allocate space to store the unique keys and their counts on the device
    thrust::device_vector<uint64_t> uniqueHashsD(numKmers);
    thrust::device_vector<uint8_t> uniqueHashsCountsD(numKmers);

    // Calculate the size of the grid and thread blocks
    int blockSize = 32;
    uint64_t numBlocks = (numKmers / blockSize) + 1;

    for (int i = 0; i < _numHashes; ++i) {
        // Fill the vectors with initial values
        thrust::fill(uniqueHashsD.begin(), uniqueHashsD.end(), 0);
        thrust::fill(uniqueHashsCountsD.begin(), uniqueHashsCountsD.end(), 0);

        // Calculate the hash value of the k-mer on GPU
        MurmurHash3_x64_128_kernel<<<numBlocks, blockSize>>>(kmersD, numKmers, sizeof(uint64_t), _seeds[i], thrust::raw_pointer_cast(hashValues.data()));
        gpuErrchk(cudaPeekAtLastError());

        // Perform parallel quicksort on GPU
        thrust::sort(hashValues.begin(), hashValues.end());

        // Use reduce_by_key to compress and calculate counts
        auto end = thrust::reduce_by_key(hashValues.begin(), hashValues.end(), thrust::make_constant_iterator(1),
                                         uniqueHashsD.begin(), uniqueHashsCountsD.begin(), thrust::equal_to<uint64_t>(), capped_add());

        // Wait for the GPU to finish operations
        gpuErrchk(cudaDeviceSynchronize());

        // Calculate the number of unique hash
        uint32_t numUniqueHashs = end.first - uniqueHashsD.begin();

        // Use a GPU kernel to perform the counting
        add_filter_kernel<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(uniqueHashsD.data()), thrust::raw_pointer_cast(uniqueHashsCountsD.data()), numUniqueHashs, _filterD, _size);
        gpuErrchk(cudaDeviceSynchronize());
    }
}

/*
 * Variable: kmersD, findBoolsD
 * -----------------
 * This variable is declared on the device (GPU) and is only accessible 
 * from the device code (i.e., CUDA kernels). It cannot be directly accessed 
 * from the host (CPU) code.
*/
void BloomFilterKernel::find_kernel(uint64_t * kmersD, uint32_t numKmers, bool * findBoolsD) {
    // Initialize findBoolsD to all true on the GPU
    gpuErrchk(cudaMemset(findBoolsD, true, numKmers * sizeof(bool)));

    // Allocate unified memory to store the hash results
    uint64_t * hashValues;
    gpuErrchk(cudaMallocManaged(&hashValues, sizeof(uint64_t) * numKmers));

    // Calculate the size of the grid and thread blocks
    int blockSize = 32;
    uint64_t numBlocks = (numKmers / blockSize) + 1;

    for (int i = 0; i < _numHashes; ++i) {
        // Calculate the hash value of the k-mer on GPU
        MurmurHash3_x64_128_kernel<<<numBlocks, blockSize>>>(kmersD, numKmers, sizeof(uint64_t), _seeds[i], hashValues);
        gpuErrchk(cudaPeekAtLastError());

        // Wait for the GPU to finish operations
        gpuErrchk(cudaDeviceSynchronize());

        // Determine whether the k-mer exists in the bloom filter
        check_filter_kernel<<<numBlocks, blockSize>>>(hashValues, numKmers, _filterD, _size, findBoolsD);
        gpuErrchk(cudaPeekAtLastError());

        // Wait for the GPU to finish operations
        gpuErrchk(cudaDeviceSynchronize());
    }

    // Free unified memory
    gpuErrchk(cudaFree(hashValues));
}

/*
 * Variable: kmersD, countsD
 * -----------------
 * This variable is declared on the device (GPU) and is only accessible 
 * from the device code (i.e., CUDA kernels). It cannot be directly accessed 
 * from the host (CPU) code.
*/
void BloomFilterKernel::count_kernel(uint64_t * kmersD, uint32_t numKmers, uint8_t * countsD) {
    // Initialize countsD to all UINT8_MAX on the GPU
    gpuErrchk(cudaMemset(countsD, 255, numKmers * sizeof(uint8_t)));

    // Allocate unified memory to store the hash results
    uint64_t * hashValues;
    gpuErrchk(cudaMallocManaged(&hashValues, sizeof(uint64_t) * numKmers));

    // Calculate the size of the grid and thread blocks
    int blockSize = 32;
    uint64_t numBlocks = (numKmers / blockSize) + 1;

    for (int i = 0; i < _numHashes; ++i) {
        // Calculate the hash value of the k-mer on GPU
        MurmurHash3_x64_128_kernel<<<numBlocks, blockSize>>>(kmersD, numKmers, sizeof(uint64_t), _seeds[i], hashValues);
        gpuErrchk(cudaPeekAtLastError());

        // Wait for the GPU to finish operations
        gpuErrchk(cudaDeviceSynchronize());

        // Determine whether the k-mer exists in the bloom filter
        count_filter_kernel<<<numBlocks, blockSize>>>(hashValues, numKmers, _filterD, _size, countsD);
        gpuErrchk(cudaPeekAtLastError());

        // Wait for the GPU to finish operations
        gpuErrchk(cudaDeviceSynchronize());
    }

    // Free unified memory
    gpuErrchk(cudaFree(hashValues));
}